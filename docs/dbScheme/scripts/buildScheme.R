#!/usr/bin/env Rscript
# Prism dbScheme: build the scheme tables and the workbook.
#
#   Rscript buildScheme.R --datamodel DIR --out-dir OUT [--scheme-dir DIR]
#
# --datamodel   ORANGE data model sources, from fetchOrangeDatamodel.R.
# --out-dir     A finished oncoanalyser output folder (one folder per
#               group). Used to see which fields are filled in matched and
#               tumor-only samples. Only all/some/none per field is
#               written; no sample names and no values.
# --scheme-dir  The folder with curated/; tables/ and dbScheme.xlsx are
#               written there. Default: the parent of this script's folder.
#
# Inputs in curated/ (edited by hand):
#   orange_descriptions.csv  description and source of ORANGE fields,
#                            keyed by data model type and JSON key
#   group_a_fields.csv       every field of the Group A (one value set per
#                            sample) tool files, with its ORANGE path
# Outputs:
#   tables/orange_fields.csv  every ORANGE JSON path
#   tables/orange_enums.csv   allowed values of the enum fields
#   tables/group_a_fields.csv group_a_fields.csv plus checks against OUT
#   dbScheme.xlsx             the three tables and a readme sheet

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
  library(jsonlite)
  library(openxlsx)
})

# ---- arguments ---------------------------------------------------------

parse_args <- function(args) {
  opts <- list(datamodel = NULL, out_dir = NULL, scheme_dir = NULL)
  i <- 1
  while (i <= length(args)) {
    key <- args[i] |> str_remove("^--") |> str_replace_all("-", "_")
    if (!key %in% names(opts) || i == length(args)) {
      stop(str_glue("Bad argument: {args[i]}"), call. = FALSE)
    }
    opts[[key]] <- args[i + 1]
    i <- i + 2
  }
  if (is.null(opts$datamodel) || is.null(opts$out_dir)) {
    stop("--datamodel and --out-dir are required", call. = FALSE)
  }
  opts
}

script_dir <- function() {
  file_arg <- commandArgs(trailingOnly = FALSE) |> str_subset("^--file=")
  if (length(file_arg) == 0) {
    return(path_wd())
  }
  file_arg[1] |> str_remove("^--file=") |> path_abs() |> path_dir()
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
model_dir <- path_abs(opts$datamodel)
out_dir <- path_abs(opts$out_dir)
scheme_dir <- if (is.null(opts$scheme_dir)) {
  path_dir(script_dir())
} else {
  path_abs(opts$scheme_dir)
}
curated_dir <- path(scheme_dir, "curated")
tables_dir <- path(scheme_dir, "tables")
dir_create(tables_dir)

log_msg <- function(...) {
  message(str_glue(..., .envir = parent.frame()))
}

warn_state <- new.env()
warn_state$n <- 0L
warn_msg <- function(...) {
  warn_state$n <- warn_state$n + 1L
  message("WARNING: ", str_glue(..., .envir = parent.frame()))
}

# ---- data model --------------------------------------------------------

NUMBER_TYPES <- c("double", "Double", "float", "Float")
INTEGER_TYPES <- c("int", "Integer", "long", "Long")
BOOLEAN_TYPES <- c("boolean", "Boolean")
STRING_TYPES <- c("String")
DATE_TYPES <- c("LocalDate")

# An abstract accessor: annotations, type, name, "();".
METHOD_RE <- regex(str_c(
  "((?:@[\\w.]+(?:\\([^)]*\\))?\\s+)*)",
  "([\\w.]+(?:<[\\w.<>, ?]+>)?)\\s+(\\w+)\\s*\\(\\s*\\)\\s*;"
))
# A derived attribute: computed by the class, may still be serialized.
DERIVED_RE <- regex(
  "@Value\\.Derived\\s+default\\s+([\\w.]+(?:<[\\w.<>, ?]+>)?)\\s+(\\w+)\\s*\\(\\s*\\)"
)

strip_comments <- function(text) {
  text |>
    str_remove_all(regex("/\\*.*?\\*/", dotall = TRUE)) |>
    str_remove_all(regex("//[^\n]*", multiline = TRUE))
}

parse_enum <- function(text) {
  m <- str_match(text, regex("\\benum\\s+(\\w+)\\s*\\{(.*?)(?:;|\\})", dotall = TRUE))
  values <- m[1, 3] |>
    str_remove_all("\\([^)]*\\)") |>
    str_split_1(",") |>
    str_squish()
  list(name = m[1, 2], values = values[values != ""])
}

parse_interface <- function(text) {
  m <- str_match(
    text,
    regex("\\binterface\\s+(\\w+)(?:\\s+extends\\s+([\\w\\s,]+?))?\\s*\\{(.*)\\}", dotall = TRUE)
  )
  body <- m[1, 4]
  methods <- str_match_all(body, METHOD_RE)[[1]]
  fields <- tibble(
    annot = methods[, 2], java_type = methods[, 3], name = methods[, 4]
  ) |>
    filter(!java_type %in% c("return", "new", "throw")) |>
    mutate(
      nullable = str_detect(annot, "@Nullable\\b"),
      json_key = coalesce(str_match(annot, "@SerializedName\\(\"([^\"]+)\"\\)")[, 2], name),
      derived = FALSE
    )
  derived <- str_match_all(body, DERIVED_RE)[[1]]
  if (nrow(derived) > 0) {
    fields <- bind_rows(fields, tibble(
      annot = "", java_type = derived[, 2], name = derived[, 3],
      nullable = FALSE, json_key = derived[, 3], derived = TRUE
    ))
  }
  extends <- if (is.na(m[1, 3])) character() else m[1, 3] |> str_split_1(",") |> str_squish()
  list(name = m[1, 2], extends = extends, fields = select(fields, -annot))
}

java_files <- dir_ls(model_dir, glob = "*.java")
if (length(java_files) == 0) {
  stop(str_glue("No .java files in {model_dir}"), call. = FALSE)
}
sources <- java_files |> map_chr(\(f) read_file(f) |> strip_comments())
enums <- sources |>
  keep(\(s) str_detect(s, "\\bpublic\\s+enum\\s")) |>
  map(parse_enum) |>
  unname()
enums <- set_names(enums, map_chr(enums, "name"))
types <- sources |>
  keep(\(s) str_detect(s, "\\bpublic\\s+interface\\s")) |>
  map(parse_interface) |>
  unname()
types <- set_names(types, map_chr(types, "name"))
model_tag <- if (file_exists(path(model_dir, "SOURCE.txt"))) {
  read_lines(path(model_dir, "SOURCE.txt")) |> str_subset("^tag: ") |> str_remove("^tag: ")
} else {
  "unknown"
}
log_msg("Data model {model_tag}: {length(types)} types, {length(enums)} enums")

all_fields <- function(type) {
  t <- types[[type]]
  inherited <- t$extends |> keep(\(p) p %in% names(types)) |> map(all_fields) |> list_rbind()
  bind_rows(inherited, t$fields) |> distinct(json_key, .keep_all = TRUE)
}

type_info <- function(java_type) {
  inner <- str_match(java_type, "^(?:List|Set|Collection)<(.+)>$")[, 2]
  is_list <- !is.na(inner)
  base <- (if (is_list) inner else java_type) |> str_remove("^.*\\.")
  kind <- case_when(
    base %in% NUMBER_TYPES ~ "number",
    base %in% INTEGER_TYPES ~ "integer",
    base %in% BOOLEAN_TYPES ~ "boolean",
    base %in% STRING_TYPES ~ "string",
    base %in% DATE_TYPES ~ "date",
    base %in% names(enums) ~ "enum",
    base %in% names(types) ~ "object",
    .default = "unknown"
  )
  list(is_list = is_list, base = base, kind = kind)
}

# One row per JSON path, depth first in declaration order. LEVEL is the
# list a path belongs to ("sample" if none).
expand <- function(type, prefix = "", level = "sample", depth = 0) {
  if (depth > 12) {
    stop(str_glue("Data model nesting too deep at {prefix}"), call. = FALSE)
  }
  fields <- all_fields(type)
  pmap(fields, \(java_type, json_key, nullable, derived, ...) {
    info <- type_info(java_type)
    path <- if (prefix == "") json_key else str_c(prefix, ".", json_key)
    row <- tibble(
      path = path, parent_type = type, field = json_key, java_type = java_type,
      json_type = if (info$is_list) "array" else info$kind,
      item_type = case_when(
        !info$is_list ~ NA_character_,
        info$kind == "object" ~ info$base,
        .default = info$kind
      ),
      enum = if (info$kind == "enum") info$base else NA_character_,
      nullable = nullable, derived = derived, level = level
    )
    if (info$kind != "object") {
      return(row)
    }
    child_prefix <- if (info$is_list) str_c(path, "[]") else path
    child_level <- if (info$is_list) path else level
    bind_rows(row, expand(info$base, child_prefix, child_level, depth + 1))
  }) |>
    list_rbind()
}

declared <- expand("OrangeRecord") |> mutate(in_datamodel = TRUE)
unknown_types <- declared |> filter(json_type == "unknown" | item_type %in% "unknown")
if (nrow(unknown_types) > 0) {
  warn_msg("Unresolved data model types: {str_c(unique(unknown_types$java_type), collapse = ', ')}")
}

# ---- observed ORANGE JSON ----------------------------------------------

# Every key path in one JSON object, with whether its value is filled
# (not null, not an empty list). Arrays of objects add "[]".
observe_json <- function(node, prefix = "") {
  if (!is.list(node) || is.null(names(node))) {
    return(tibble(path = character(), filled = logical()))
  }
  imap(node, \(value, key) {
    path <- if (prefix == "") key else str_c(prefix, ".", key)
    filled <- !is.null(value) && !(is.list(value) && length(value) == 0)
    children <- if (is.list(value) && is.null(names(value))) {
      value |>
        keep(\(v) is.list(v) && !is.null(names(v))) |>
        map(\(v) observe_json(v, str_c(path, "[]"))) |>
        list_rbind()
    } else if (is.list(value)) {
      observe_json(value, path)
    }
    bind_rows(tibble(path = path, filled = filled), children)
  }) |>
    list_rbind()
}

orange_files <- dir_ls(out_dir, recurse = 2, glob = "*/orange/*.orange.json")
if (length(orange_files) == 0) {
  stop(str_glue("No ORANGE JSON files under {out_dir}"), call. = FALSE)
}
samples <- map(orange_files, \(f) {
  json <- read_json(f, simplifyVector = FALSE)
  list(
    group_dir = path_dir(path_dir(f)),
    tumor = json$sampleId,
    normal = json$referenceId %||% NA_character_,
    mode = if (is.null(json$referenceId)) "tumor_only" else "matched",
    observed = observe_json(json) |> summarize(filled = any(filled), .by = path)
  )
})
n_mode <- c(
  matched = sum(map_chr(samples, "mode") == "matched"),
  tumor_only = sum(map_chr(samples, "mode") == "tumor_only")
)
log_msg(
  "{length(samples)} ORANGE files: {n_mode['matched']} matched, ",
  "{n_mode['tumor_only']} tumor-only"
)

fill_status <- function(n_filled, n_total) {
  case_when(
    n_total == 0 ~ "no samples",
    n_filled == n_total ~ "all",
    n_filled == 0 ~ "none",
    .default = "some"
  )
}

observed <- samples |>
  map(\(s) mutate(s$observed, mode = s$mode)) |>
  list_rbind() |>
  summarize(n_filled = sum(filled), .by = c(path, mode)) |>
  pivot_wider(names_from = mode, values_from = n_filled, names_prefix = "n_")
for (col in c("n_matched", "n_tumor_only")) {
  if (!col %in% names(observed)) observed[[col]] <- 0L
}
observed <- observed |>
  mutate(
    in_json = TRUE,
    filled_matched = fill_status(coalesce(n_matched, 0L), n_mode[["matched"]]),
    filled_tumor_only = fill_status(coalesce(n_tumor_only, 0L), n_mode[["tumor_only"]])
  ) |>
  select(path, in_json, filled_matched, filled_tumor_only)

# ---- ORANGE fields table -----------------------------------------------

# Group of a path, as in the per-sample readout catalog: A one value set
# per sample, B somatic events, C germline. First matching prefix wins.
GROUP_RULES <- tribble(
  ~prefix, ~group,
  "plots", "display",
  "isofox", "RNA",
  "purple.germlineDrivers", "C",
  "purple.germlineVariants", "C",
  "purple.germlineGainsDels", "C",
  "linx.germlineBreakends", "C",
  "purple.somaticDrivers", "B",
  "purple.somaticVariants", "B",
  "purple.somaticGainsDels", "B",
  "purple.armCopyNumberAbberations", "B",
  "linx", "B"
)

assign_group <- function(paths) {
  map_chr(paths, \(p) {
    hit <- GROUP_RULES |>
      filter(
        p == prefix |
          str_starts(p, fixed(str_c(prefix, "."))) |
          str_starts(p, fixed(str_c(prefix, "[]")))
      )
    if (nrow(hit) == 0) "A" else hit$group[1]
  })
}

descriptions <- read_csv(
  path(curated_dir, "orange_descriptions.csv"),
  col_types = cols(.default = col_character()), show_col_types = FALSE
)
dup_desc <- descriptions |> count(type, field) |> filter(n > 1)
if (nrow(dup_desc) > 0) {
  warn_msg("Duplicate descriptions: {str_c(dup_desc$type, '.', dup_desc$field, collapse = ', ')}")
}

orange_fields <- declared |>
  full_join(observed, by = "path") |>
  mutate(
    order = row_number(),
    in_datamodel = coalesce(in_datamodel, FALSE),
    in_json = coalesce(in_json, FALSE),
    filled_matched = coalesce(filled_matched, fill_status(0L, n_mode[["matched"]])),
    filled_tumor_only = coalesce(filled_tumor_only, fill_status(0L, n_mode[["tumor_only"]])),
    group = assign_group(path)
  ) |>
  left_join(descriptions, by = c(parent_type = "type", field = "field")) |>
  select(
    order, path, group, level, json_type, item_type, enum, nullable, derived,
    in_datamodel, in_json, filled_matched, filled_tumor_only,
    description, source_file, source_column, parent_type, field, java_type
  )

unused_desc <- descriptions |>
  anti_join(orange_fields, by = c(type = "parent_type", field = "field"))
if (nrow(unused_desc) > 0) {
  keys <- str_c(unused_desc$type, ".", unused_desc$field, collapse = ", ")
  warn_msg("Descriptions with no data model field: {keys}")
}
not_in_model <- orange_fields |> filter(!in_datamodel)
if (nrow(not_in_model) > 0) {
  warn_msg("JSON paths not in the data model: {str_c(not_in_model$path, collapse = ', ')}")
}
undescribed_a <- orange_fields |> filter(group == "A", is.na(description))
if (nrow(undescribed_a) > 0) {
  warn_msg("Group A paths without a description: {str_c(undescribed_a$path, collapse = ', ')}")
}

orange_enums <- orange_fields |>
  filter(!is.na(enum)) |>
  summarize(used_by = str_c(path, collapse = "; "), .by = enum) |>
  mutate(value = map(enum, \(e) enums[[e]]$values)) |>
  unnest_longer(value) |>
  select(enum, value, used_by)

# ---- Group A tool files ------------------------------------------------

group_a <- read_csv(
  path(curated_dir, "group_a_fields.csv"),
  col_types = cols(.default = col_character()), show_col_types = FALSE
)

# The files PATTERN names in one group folder: <T> tumor, <N> normal,
# <S> each of the two.
resolve_files <- function(sample, pattern) {
  ids <- c(sample$tumor, sample$normal) |> discard(is.na)
  candidates <- if (str_detect(pattern, "<S>")) {
    map_chr(ids, \(id) str_replace_all(pattern, "<S>", id))
  } else {
    pattern
  }
  candidates <- candidates |> str_replace_all("<T>", sample$tumor)
  if (is.na(sample$normal)) {
    candidates <- candidates |> discard(\(x) str_detect(x, "<N>"))
  } else {
    candidates <- candidates |> str_replace_all("<N>", sample$normal)
  }
  files <- path(sample$group_dir, candidates)
  files[file_exists(files)]
}

# The field names a file holds under LAYOUT.
file_fields <- function(file, layout) {
  lines <- read_lines(file, lazy = FALSE)
  lines <- lines[lines != ""]
  if (length(lines) == 0) {
    return(character())
  }
  if (layout %in% c("wide", "long")) {
    return(str_split_1(lines[1], "\t"))
  }
  if (layout == "key-value") {
    return(str_extract(lines, "^[^\t]+"))
  }
  if (str_starts(layout, "row:")) {
    col <- str_remove(layout, "^row:")
    header <- str_split_1(lines[1], "\t")
    idx <- match(col, header)
    if (is.na(idx) || length(lines) < 2) {
      return(character())
    }
    return(lines[-1] |> str_split("\t") |> map_chr(\(x) x[idx]) |> unique())
  }
  if (layout == "text") {
    return("(samtools flagstat lines)")
  }
  stop(str_glue("Unknown layout {layout}"), call. = FALSE)
}

file_layouts <- group_a |> distinct(file, layout)
found <- samples |>
  map(\(s) {
    pmap(file_layouts, \(file, layout) {
      files <- resolve_files(s, file)
      tibble(
        file = file, layout = layout, mode = s$mode,
        field = list(unique(unlist(map(files, \(f) file_fields(f, layout)))))
      )
    }) |>
      list_rbind()
  }) |>
  list_rbind() |>
  mutate(sample = rep(seq_along(samples), each = nrow(file_layouts)))

presence <- group_a |>
  select(file, layout, field) |>
  inner_join(
    found |> rename(fields = field),
    by = c("file", "layout"), relationship = "many-to-many"
  ) |>
  mutate(present = map2_lgl(field, fields, \(f, fs) f %in% fs)) |>
  summarize(n_present = sum(present), .by = c(file, layout, field, mode)) |>
  pivot_wider(names_from = mode, values_from = n_present, names_prefix = "n_")
for (col in c("n_matched", "n_tumor_only")) {
  if (!col %in% names(presence)) presence[[col]] <- 0L
}

undocumented <- found |>
  filter(layout != "text") |>
  unnest_longer(field) |>
  filter(!is.na(field)) |>
  distinct(file, layout, field) |>
  anti_join(group_a, by = c("file", "layout", "field"))
if (nrow(undocumented) > 0) {
  walk(seq_len(nrow(undocumented)), \(i) {
    warn_msg("Field not in group_a_fields.csv: {undocumented$file[i]} {undocumented$field[i]}")
  })
}

bad_paths <- group_a |>
  filter(!is.na(orange_path), !orange_path %in% orange_fields$path)
if (nrow(bad_paths) > 0) {
  warn_msg("orange_path not in ORANGE: {str_c(unique(bad_paths$orange_path), collapse = ', ')}")
}

group_a_fields <- group_a |>
  left_join(presence, by = c("file", "layout", "field")) |>
  mutate(
    in_orange = if_else(is.na(orange_path), "no", "yes"),
    present_matched = fill_status(coalesce(n_matched, 0L), n_mode[["matched"]]),
    present_tumor_only = fill_status(coalesce(n_tumor_only, 0L), n_mode[["tumor_only"]])
  ) |>
  select(
    readout, file, layout, field, description, in_orange, orange_path,
    present_matched, present_tumor_only, notes
  )

# ---- write -------------------------------------------------------------

write_csv(orange_fields, path(tables_dir, "orange_fields.csv"), na = "")
write_csv(orange_enums, path(tables_dir, "orange_enums.csv"), na = "")
write_csv(group_a_fields, path(tables_dir, "group_a_fields.csv"), na = "")

readme <- tribble(
  ~sheet, ~column, ~description,
  "(all)", "", str_glue(
    "ORANGE data model {model_tag}. Filled/present columns from ",
    "{n_mode['matched']} matched and {n_mode['tumor_only']} tumor-only samples: ",
    "all, some or none of them."
  ),
  "group_a_fields", "", "Every field of the per-sample summary (Group A) tool files and its ORANGE path.",
  "group_a_fields", "readout", "Short name of the readout.",
  "group_a_fields", "file", "File under out/<group>/. <T> tumor ID, <N> normal ID, <S> either (one file per BAM).",
  "group_a_fields", "layout", "wide: one row, one column per field. long: one row per item. key-value: two columns, no header. row:COL: the field is a value of column COL. text: not a table.",
  "group_a_fields", "field", "Column name, key or row value.",
  "group_a_fields", "description", "What the field holds (from the tool README where it has one).",
  "group_a_fields", "in_orange", "yes if ORANGE holds the same value.",
  "group_a_fields", "orange_path", "Path of that value in the ORANGE JSON (see orange_fields).",
  "group_a_fields", "present_matched", "Whether the field is in the files of the matched samples.",
  "group_a_fields", "present_tumor_only", "Whether the field is in the files of the tumor-only samples.",
  "group_a_fields", "notes", "Differences from the tool README and other remarks.",
  "orange_fields", "", "Every path of the ORANGE JSON: from the data model, checked against the JSON files.",
  "orange_fields", "order", "Order in the data model (depth first).",
  "orange_fields", "path", "JSON path. [] marks a list; fields after it belong to each list item.",
  "orange_fields", "group", "A: one value set per sample. B: somatic events. C: germline. RNA: needs RNA data. display: plot file paths.",
  "orange_fields", "level", "sample: one value per sample. Otherwise the list the path belongs to.",
  "orange_fields", "json_type", "string, number, integer, boolean, date, enum, object or array.",
  "orange_fields", "item_type", "For arrays: type of the items.",
  "orange_fields", "enum", "For enum fields: the enum name (values in orange_enums).",
  "orange_fields", "nullable", "TRUE if the data model allows null.",
  "orange_fields", "derived", "TRUE if computed from other fields of the same record.",
  "orange_fields", "in_datamodel", "FALSE if the path is in the JSON files but not in the data model.",
  "orange_fields", "in_json", "FALSE if the path is in no JSON file (its parent is null or empty everywhere).",
  "orange_fields", "filled_matched", "Whether the value is filled (not null, not empty) in the matched samples.",
  "orange_fields", "filled_tumor_only", "Whether the value is filled in the tumor-only samples.",
  "orange_fields", "description", "What the field holds. Empty for groups B and C (not yet documented).",
  "orange_fields", "source_file", "Tool file the value comes from; (ORANGE) if ORANGE computes it.",
  "orange_fields", "source_column", "Column or key in that file.",
  "orange_fields", "parent_type", "Data model type that holds the field.",
  "orange_fields", "field", "JSON key.",
  "orange_fields", "java_type", "Type in the data model.",
  "orange_enums", "", "Allowed values of each enum in the data model.",
  "orange_enums", "enum", "Enum name.",
  "orange_enums", "value", "Allowed value.",
  "orange_enums", "used_by", "ORANGE paths with this enum."
)

TEXT_COLS <- c("description", "notes", "used_by")
TEXT_WIDTH <- 60
MIN_COL_WIDTH <- 6
MAX_COL_WIDTH <- 50

auto_widths <- function(df) {
  content <- map_dbl(df, \(x) {
    n <- nchar(as.character(x))
    n <- n[!is.na(n) & n > 0]
    if (length(n) == 0) 0 else ceiling(quantile(n, 0.95, names = FALSE)) + 2
  })
  header <- nchar(names(df)) + 4
  pmin(pmax(content, header, MIN_COL_WIDTH), MAX_COL_WIDTH) |> unname()
}

write_workbook <- function(sheets, file) {
  wb <- createWorkbook()
  modifyBaseFont(wb, fontSize = 10, fontName = "Arial")
  header_style <- createStyle(textDecoration = "bold", fgFill = "#DDDDDD")
  wrap_style <- createStyle(wrapText = TRUE, valign = "top")
  iwalk(sheets, \(df, name) {
    addWorksheet(wb, name)
    writeData(wb, name, df, headerStyle = header_style, keepNA = FALSE)
    freezePane(wb, name, firstRow = TRUE)
    if (nrow(df) > 0) {
      addFilter(wb, name, rows = 1, cols = seq_len(ncol(df)))
    }
    widths <- auto_widths(df)
    text <- which(names(df) %in% TEXT_COLS | name == "readme" & names(df) == "description")
    widths[text] <- TEXT_WIDTH
    setColWidths(wb, name, cols = seq_len(ncol(df)), widths = widths)
    if (length(text) > 0 && nrow(df) > 0) {
      addStyle(wb, name, wrap_style, rows = seq_len(nrow(df)) + 1, cols = text, gridExpand = TRUE)
    }
  })
  # Explicit window size, so Excel does not open the workbook full screen.
  wb$workbook$bookViews <- str_c(
    '<bookViews><workbookView xWindow="4000" yWindow="3000" ',
    'windowWidth="32000" windowHeight="20000" firstSheet="0" ',
    'activeTab="0"/></bookViews>'
  )
  saveWorkbook(wb, file, overwrite = TRUE)
}

workbook <- path(scheme_dir, "dbScheme.xlsx")
write_workbook(
  list(
    readme = readme, group_a_fields = group_a_fields,
    orange_fields = orange_fields, orange_enums = orange_enums
  ),
  workbook
)

log_msg("orange_fields: {nrow(orange_fields)} paths ({sum(orange_fields$group == 'A')} in group A)")
log_msg("orange_enums: {nrow(orange_enums)} values")
n_in_orange <- sum(group_a_fields$in_orange == "yes")
log_msg("group_a_fields: {nrow(group_a_fields)} fields, {n_in_orange} in ORANGE")
log_msg("Wrote {tables_dir} and {workbook}")
log_msg("{warn_state$n} warnings")
