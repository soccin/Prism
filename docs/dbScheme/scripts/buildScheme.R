#!/usr/bin/env Rscript
# Prism dbScheme: build the scheme tables and the workbook.
#
#   Rscript buildScheme.R --datamodel DIR --out-dir OUT [--scheme-dir DIR]
#
# --datamodel   ORANGE data model sources, from fetchOrangeDatamodel.R.
# --out-dir     A finished oncoanalyser output folder (one folder per
#               group). Used to see which fields are present in matched and
#               tumor-only samples and to check the copies. Only all/some/
#               none and counts are written; no sample names, no values.
# --scheme-dir  The folder with curated/; tables/ and dbScheme.xlsx are
#               written there. Default: the parent of this script's folder.
#
# Inputs in curated/ (edited by hand):
#   files.csv       one row per per-sample file: tool, role (primary or
#                   secondary), what it holds
#   fields.csv      every field of the files with detail "fields": role
#                   (primary, copy, derived, annotation, key), description,
#                   source of copies and derived values
#   orange_map.csv  where each part of the ORANGE JSON comes from, and the
#                   values ORANGE computes itself
# Outputs:
#   tables/files.csv   files.csv plus field counts by role
#   tables/fields.csv  fields.csv plus presence and the copy checks
#   tables/orange.csv  orange_map.csv plus JSON type, fill status and the
#                      number of JSON paths each row covers
#   dbScheme.xlsx      the three tables and a readme sheet

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

read_curated <- function(name) {
  read_csv(
    path(curated_dir, name),
    col_types = cols(.default = col_character()), show_col_types = FALSE
  )
}

FIELD_ROLES <- c("primary", "copy", "derived", "annotation", "key")
FILE_ROLES <- c("primary", "secondary")
MAP_ROLES <- c("copy", "derived", "annotation", "metadata", "plot", "rna")

# ---- curated tables ----------------------------------------------------

files <- read_curated("files.csv")
fields <- read_curated("fields.csv")
orange_map <- read_curated("orange_map.csv")

bad <- files |> filter(!role %in% FILE_ROLES)
if (nrow(bad) > 0) warn_msg("files.csv: bad role for {str_c(bad$file, collapse = ', ')}")
bad <- fields |> filter(!role %in% FIELD_ROLES)
if (nrow(bad) > 0) warn_msg("fields.csv: bad role for {str_c(bad$field, collapse = ', ')}")
bad <- orange_map |> filter(!role %in% MAP_ROLES)
if (nrow(bad) > 0) warn_msg("orange_map.csv: bad role for {str_c(bad$path, collapse = ', ')}")

no_source <- fields |> filter(role %in% c("copy", "derived", "annotation"), is.na(source))
if (nrow(no_source) > 0) {
  warn_msg("fields.csv: no source for {str_c(no_source$field, collapse = ', ')}")
}
dup <- fields |> count(file, layout, field) |> filter(n > 1)
if (nrow(dup) > 0) warn_msg("fields.csv: duplicate rows for {str_c(dup$field, collapse = ', ')}")

detailed <- files |> filter(detail == "fields") |> pull(file)
orphan <- setdiff(unique(fields$file), detailed)
if (length(orphan) > 0) {
  warn_msg("fields.csv files not marked detail = fields in files.csv: {str_c(orphan, collapse = ', ')}")
}
missing <- setdiff(detailed, fields$file)
if (length(missing) > 0) warn_msg("No fields for: {str_c(missing, collapse = ', ')}")

role_check <- fields |>
  summarize(n_primary = sum(role == "primary"), .by = file) |>
  inner_join(select(files, file, file_role = role), by = "file") |>
  filter((file_role == "primary") != (n_primary > 0))
if (nrow(role_check) > 0) {
  warn_msg("File role does not match its field roles: {str_c(role_check$file, collapse = ', ')}")
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

parse_enum_name <- function(text) {
  str_match(text, "\\benum\\s+(\\w+)")[1, 2]
}

parse_interface <- function(text) {
  m <- str_match(
    text,
    regex("\\binterface\\s+(\\w+)(?:\\s+extends\\s+([\\w\\s,]+?))?\\s*\\{(.*)\\}", dotall = TRUE)
  )
  body <- m[1, 4]
  methods <- str_match_all(body, METHOD_RE)[[1]]
  fields <- tibble(annot = methods[, 2], java_type = methods[, 3], name = methods[, 4]) |>
    filter(!java_type %in% c("return", "new", "throw")) |>
    mutate(json_key = coalesce(str_match(annot, "@SerializedName\\(\"([^\"]+)\"\\)")[, 2], name))
  derived <- str_match_all(body, DERIVED_RE)[[1]]
  if (nrow(derived) > 0) {
    fields <- bind_rows(fields, tibble(
      annot = "", java_type = derived[, 2], name = derived[, 3], json_key = derived[, 3]
    ))
  }
  extends <- if (is.na(m[1, 3])) character() else m[1, 3] |> str_split_1(",") |> str_squish()
  list(name = m[1, 2], extends = extends, fields = select(fields, java_type, json_key))
}

java_files <- dir_ls(model_dir, glob = "*.java")
if (length(java_files) == 0) {
  stop(str_glue("No .java files in {model_dir}"), call. = FALSE)
}
sources <- java_files |> map_chr(\(f) read_file(f) |> strip_comments())
enum_names <- sources |>
  keep(\(s) str_detect(s, "\\bpublic\\s+enum\\s")) |>
  map_chr(parse_enum_name) |>
  unname()
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
log_msg("Data model {model_tag}: {length(types)} types, {length(enum_names)} enums")

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
    base %in% enum_names ~ "enum",
    base %in% names(types) ~ "object",
    .default = "unknown"
  )
  list(is_list = is_list, base = base, kind = kind)
}

# One row per JSON path, depth first in declaration order. "[]" marks a
# list; paths after it belong to each item.
expand <- function(type, prefix = "", depth = 0) {
  if (depth > 12) {
    stop(str_glue("Data model nesting too deep at {prefix}"), call. = FALSE)
  }
  pmap(all_fields(type), \(java_type, json_key) {
    info <- type_info(java_type)
    path <- if (prefix == "") json_key else str_c(prefix, ".", json_key)
    row <- tibble(path = path, json_type = if (info$is_list) "array" else info$kind)
    if (info$kind != "object") {
      return(row)
    }
    child_prefix <- if (info$is_list) str_c(path, "[]") else path
    bind_rows(row, expand(info$base, child_prefix, depth + 1))
  }) |>
    list_rbind()
}

declared <- expand("OrangeRecord") |> mutate(order = row_number())
if (any(declared$json_type == "unknown")) {
  warn_msg("Unresolved data model types at: {str_c(declared$path[declared$json_type == 'unknown'], collapse = ', ')}")
}

# ---- samples -----------------------------------------------------------

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
}) |>
  unname()
n_mode <- c(
  matched = sum(map_chr(samples, "mode") == "matched"),
  tumor_only = sum(map_chr(samples, "mode") == "tumor_only")
)
log_msg("{length(samples)} samples: {n_mode['matched']} matched, {n_mode['tumor_only']} tumor-only")

fill_status <- function(n_filled, n_total) {
  case_when(
    n_total == 0 ~ "no samples",
    n_filled == n_total ~ "all",
    n_filled == 0 ~ "none",
    .default = "some"
  )
}

# ---- reading the tool files --------------------------------------------

# The files PATTERN names in one group folder, named by role: T for <T>
# and for the tumor's <S> file, N for <N> and the normal's <S> file.
resolve_files <- function(sample, pattern) {
  ids <- c(T = sample$tumor, N = sample$normal) |> discard(is.na)
  candidates <- if (str_detect(pattern, "<S>")) {
    map_chr(ids, \(id) str_replace_all(pattern, "<S>", id))
  } else if (str_detect(pattern, "<N>")) {
    if (is.na(sample$normal)) character() else c(N = str_replace_all(pattern, "<N>", sample$normal))
  } else {
    c(T = pattern)
  }
  candidates <- candidates |> str_replace_all("<T>", sample$tumor) |> set_names(names(candidates))
  files <- path(sample$group_dir, candidates) |> set_names(names(candidates))
  files[file_exists(files)]
}

read_tab <- function(file, col_names = TRUE) {
  read_tsv(file, col_names = col_names, col_types = cols(.default = col_character()),
           show_col_types = FALSE, na = character())
}

# Line label and QC-passed count of each samtools flagstat line.
read_flagstat <- function(file) {
  lines <- read_lines(file) |> discard(\(x) x == "")
  tibble(
    field = lines |>
      str_remove("^\\d+ \\+ \\d+ ") |>
      str_remove(" \\((QC-passed|[0-9.]+%|N/A)[^)]*\\)$"),
    value = str_extract(lines, "^\\d+")
  )
}

# Field names in one file under LAYOUT.
file_fields <- function(file, layout) {
  if (layout == "flagstat") {
    return(read_flagstat(file)$field)
  }
  if (layout == "key-value") {
    return(read_tab(file, col_names = FALSE)[[1]])
  }
  tab <- read_tab(file)
  if (layout %in% c("wide", "long")) {
    return(names(tab))
  }
  if (str_starts(layout, "row:")) {
    key_col <- str_split_1(layout, ":")[2]
    return(if (key_col %in% names(tab)) unique(tab[[key_col]]) else character())
  }
  stop(str_glue("Unknown layout {layout}"), call. = FALSE)
}

# Values of one field in one file, by role (T or N). NULL for layouts
# with more than one value per field (long).
field_values <- function(file, layout, field, file_role) {
  if (layout == "flagstat") {
    v <- read_flagstat(file) |> filter(field == !!field)
    return(tibble(role = file_role, value = v$value))
  }
  if (layout == "key-value") {
    v <- read_tab(file, col_names = FALSE) |> filter(X1 == field)
    return(tibble(role = file_role, value = v$X2))
  }
  tab <- read_tab(file)
  if (layout == "wide") {
    if (!field %in% names(tab) || nrow(tab) == 0) return(tibble(role = character(), value = character()))
    return(tibble(role = file_role, value = tab[[field]][1]))
  }
  if (str_starts(layout, "row:")) {
    cols <- str_split_1(layout, ":")
    rows <- tab |> filter(.data[[cols[2]]] == field)
    role <- if ("SampleType" %in% names(rows)) {
      if_else(rows$SampleType == "NORMAL", "N", "T")
    } else {
      rep(file_role, nrow(rows))
    }
    return(tibble(role = role, value = rows[[cols[3]]]))
  }
  NULL
}

# ---- presence of each field --------------------------------------------

file_layouts <- fields |> distinct(file, layout)
found <- samples |>
  imap(\(s, i) {
    pmap(file_layouts, \(file, layout) {
      paths <- resolve_files(s, file)
      tibble(
        sample = i, file = file, layout = layout, mode = s$mode,
        field = list(unique(unlist(map(paths, \(p) file_fields(p, layout)))))
      )
    }) |>
      list_rbind()
  }) |>
  list_rbind()

presence <- fields |>
  select(file, layout, field) |>
  inner_join(found |> rename(fields_found = field), by = c("file", "layout"),
             relationship = "many-to-many") |>
  mutate(present = map2_lgl(field, fields_found, \(f, fs) f %in% fs)) |>
  summarize(n_present = sum(present), .by = c(file, layout, field, mode)) |>
  pivot_wider(names_from = mode, values_from = n_present, names_prefix = "n_")
for (col in c("n_matched", "n_tumor_only")) {
  if (!col %in% names(presence)) presence[[col]] <- 0L
}

undocumented <- found |>
  unnest_longer(field) |>
  filter(!is.na(field)) |>
  distinct(file, layout, field) |>
  anti_join(fields, by = c("file", "layout", "field"))
walk(seq_len(nrow(undocumented)), \(i) {
  warn_msg("Field not in fields.csv: {undocumented$file[i]} {undocumented$field[i]}")
})

# ---- copy checks -------------------------------------------------------

decimals <- function(x) {
  d <- str_extract(x, "(?<=\\.)\\d+$")
  if_else(is.na(d), 0L, nchar(d))
}

# Equal as numbers (to the coarser precision; an integer only if exactly
# equal to the other value) or, if not numbers, as upper-case text.
same_value <- function(a, b) {
  na <- suppressWarnings(as.numeric(a))
  nb <- suppressWarnings(as.numeric(b))
  da <- decimals(a)
  db <- decimals(b)
  d <- pmin(da, db, 4L)
  tol <- if_else(d == 0L & da != db, 1e-9, 0.5 * 10^(-d) + 1e-9)
  if_else(!is.na(na) & !is.na(nb), abs(na - nb) <= tol, str_to_upper(a) == str_to_upper(b))
}

layout_of <- fields |> distinct(file, layout, field)

check_copy <- function(file, layout, field, source) {
  src <- str_match(source, "^(.+)::(.+)$")
  if (is.na(src[1, 1])) {
    return("not checked")
  }
  src_layout <- layout_of |> filter(file == src[1, 2], field == src[1, 3]) |> pull(layout)
  if (length(src_layout) == 0) {
    return("not checked (source file not described here)")
  }
  if (layout == "long" || src_layout[1] == "long") {
    return("not checked (multi-row file)")
  }
  pairs <- map(samples, \(s) {
    get <- function(pattern, lay, fld) {
      paths <- resolve_files(s, pattern)
      imap(paths, \(p, r) field_values(p, lay, fld, r)) |> list_rbind()
    }
    a <- get(file, layout, field)
    b <- get(src[1, 2], src_layout[1], src[1, 3])
    if (is.null(a) || is.null(b) || nrow(a) == 0 || nrow(b) == 0) return(NULL)
    inner_join(a, b, by = "role", suffix = c("_a", "_b"))
  }) |>
    list_rbind() |>
    filter(value_a != "", value_b != "")
  if (nrow(pairs) == 0) {
    return("not checked (no values)")
  }
  n_same <- sum(same_value(pairs$value_a, pairs$value_b))
  if (n_same < nrow(pairs)) {
    warn_msg("Copy differs: {file} {field} vs {source} ({nrow(pairs) - n_same} of {nrow(pairs)})")
    return(str_glue("DIFFERS in {nrow(pairs) - n_same} of {nrow(pairs)}"))
  }
  str_glue("equal in {n_same} of {nrow(pairs)}")
}

fields_out <- fields |>
  left_join(presence, by = c("file", "layout", "field")) |>
  mutate(
    check = pmap_chr(list(file, layout, field, source, role), \(f, l, fl, s, r) {
      if (r == "copy") check_copy(f, l, fl, s) else ""
    }),
    present_matched = fill_status(coalesce(n_matched, 0L), n_mode[["matched"]]),
    present_tumor_only = fill_status(coalesce(n_tumor_only, 0L), n_mode[["tumor_only"]])
  ) |>
  left_join(select(files, file, tool), by = "file") |>
  select(tool, file, layout, field, role, description, source, check,
         present_matched, present_tumor_only, notes)

# ---- files table -------------------------------------------------------

role_counts <- fields |>
  count(file, role) |>
  pivot_wider(names_from = role, values_from = n, names_prefix = "n_", values_fill = 0L)
for (r in FIELD_ROLES) {
  col <- str_c("n_", r)
  if (!col %in% names(role_counts)) role_counts[[col]] <- 0L
}
files_present <- samples |>
  map(\(s) {
    tibble(file = files$file, mode = s$mode,
           found = map_lgl(files$file, \(f) length(resolve_files(s, f)) > 0))
  }) |>
  list_rbind() |>
  summarize(n = sum(found), .by = c(file, mode)) |>
  pivot_wider(names_from = mode, values_from = n, names_prefix = "n_")
for (col in c("n_matched", "n_tumor_only")) {
  if (!col %in% names(files_present)) files_present[[col]] <- 0L
}
files_out <- files |>
  left_join(role_counts, by = "file") |>
  left_join(files_present, by = "file") |>
  mutate(
    present_matched = fill_status(coalesce(n_matched, 0L), n_mode[["matched"]]),
    present_tumor_only = fill_status(coalesce(n_tumor_only, 0L), n_mode[["tumor_only"]])
  ) |>
  select(tool, file, role, holds, layout, detail, n_primary, n_copy, n_derived,
         n_annotation, n_key, matched_only, present_matched, present_tumor_only, notes)
missing_files <- files_out |> filter(present_matched == "none", present_tumor_only == "none")
if (nrow(missing_files) > 0) {
  warn_msg("Files in no sample: {str_c(missing_files$file, collapse = ', ')}")
}

# ---- ORANGE map --------------------------------------------------------

observed <- samples |>
  map(\(s) mutate(s$observed, mode = s$mode)) |>
  list_rbind() |>
  summarize(n_filled = sum(filled), .by = c(path, mode)) |>
  pivot_wider(names_from = mode, values_from = n_filled, names_prefix = "n_")
for (col in c("n_matched", "n_tumor_only")) {
  if (!col %in% names(observed)) observed[[col]] <- 0L
}

not_in_model <- setdiff(observed$path, declared$path)
if (length(not_in_model) > 0) {
  warn_msg("JSON paths not in the data model: {str_c(not_in_model, collapse = ', ')}")
}
bad_map <- setdiff(orange_map$path, declared$path)
if (length(bad_map) > 0) warn_msg("orange_map.csv paths not in the data model: {str_c(bad_map, collapse = ', ')}")

# The map row of each data model path: the longest map path that is the
# path itself or one of its ancestors.
is_under <- function(p, prefix) {
  p == prefix | str_starts(p, fixed(str_c(prefix, "."))) | str_starts(p, fixed(str_c(prefix, "[]")))
}
declared <- declared |>
  mutate(map_path = map_chr(path, \(p) {
    hits <- orange_map$path[map_lgl(orange_map$path, \(m) is_under(p, m))]
    if (length(hits) == 0) NA_character_ else hits[which.max(nchar(hits))]
  }))
# Top-level containers (purple, linx) are covered by their children.
has_mapped_child <- function(p) {
  others <- declared |> filter(path != p, !is.na(map_path))
  any(is_under(others$path, p))
}
unmapped <- declared |>
  filter(is.na(map_path)) |>
  filter(!map_lgl(path, has_mapped_child))
if (nrow(unmapped) > 0) {
  warn_msg("ORANGE paths not covered by orange_map.csv: {str_c(unmapped$path, collapse = ', ')}")
}

path_types <- declared |> select(path, json_type, order)
path_counts <- declared |> count(map_path, name = "n_paths")
orange_out <- orange_map |>
  left_join(path_types, by = "path") |>
  left_join(path_counts, by = c(path = "map_path")) |>
  left_join(observed, by = "path") |>
  mutate(
    n_paths = coalesce(n_paths, 0L),
    filled_matched = fill_status(coalesce(n_matched, 0L), n_mode[["matched"]]),
    filled_tumor_only = fill_status(coalesce(n_tumor_only, 0L), n_mode[["tumor_only"]])
  ) |>
  arrange(order) |>
  select(path, role, source, rule, json_type, n_paths, filled_matched, filled_tumor_only, notes)

# ---- write -------------------------------------------------------------

write_csv(files_out, path(tables_dir, "files.csv"), na = "")
write_csv(fields_out, path(tables_dir, "fields.csv"), na = "")
write_csv(orange_out, path(tables_dir, "orange.csv"), na = "")

readme <- tribble(
  ~sheet, ~column, ~description,
  "(all)", "", str_glue(
    "Which file first holds each per-sample value. Present/filled columns from ",
    "{n_mode['matched']} matched and {n_mode['tumor_only']} tumor-only samples: all, some or none. ",
    "ORANGE data model {model_tag}."
  ),
  "(all)", "role", str_c(
    "primary: first written here, by the tool that computes it. ",
    "copy: the same value as a field of another file (renamed, rounded or recoded). ",
    "derived: computed from values that another file holds, by the rule in source. ",
    "annotation: looked up from a fixed reference table (drug links, names, thresholds), not measured. ",
    "key: identifies the row (sample, gene, classifier), holds no result."
  ),
  "files", "", "One row per per-sample file. primary: holds values no other file has. secondary: holds no primary values; only copies, derived values and annotations.",
  "files", "tool", "Tool that writes the file.",
  "files", "file", "Path under out/<group>/. <T> tumor ID, <N> normal ID, <S> either (one file per BAM).",
  "files", "role", "primary or secondary (see above).",
  "files", "holds", "What the file holds.",
  "files", "layout", "wide: one row, one column per field. long: one row per item. key-value: two columns, no header. flagstat: samtools flagstat text. json, other.",
  "files", "detail", "fields: described field by field in the fields sheet. orange: in the orange sheet. none: not described field by field.",
  "files", "n_primary", "Number of fields by role (for files with detail fields). Also n_copy, n_derived, n_annotation, n_key.",
  "files", "matched_only", "yes if the tool runs only with a matched normal.",
  "files", "present_matched", "Whether the file exists for the matched samples.",
  "files", "present_tumor_only", "Whether the file exists for the tumor-only samples.",
  "files", "notes", "Remarks.",
  "fields", "", "Every field of the files with detail fields.",
  "fields", "layout", "As in files. row:K:V: the field is a value of column K, its value is in column V.",
  "fields", "field", "Column name, key, row value or flagstat line.",
  "fields", "role", "primary, copy, derived, annotation or key (see above).",
  "fields", "description", "What the field holds (from the tool README where it has one).",
  "fields", "source", "copy: FILE::FIELD the value is copied from. derived: the rule. annotation: the reference table.",
  "fields", "check", "For copies: whether the values equal the source in every sample (tumor and normal values separately). Multi-row files are not checked by this script.",
  "fields", "present_matched", "Whether the field is in the files of the matched samples.",
  "fields", "present_tumor_only", "Whether the field is in the files of the tumor-only samples.",
  "fields", "notes", "Differences from the tool README and other remarks.",
  "orange", "", "Where each part of the ORANGE JSON comes from. A row covers its path and every path under it that has no row of its own. ORANGE adds the rows whose role is not copy.",
  "orange", "path", "JSON path. [] marks a list; paths after it belong to each item.",
  "orange", "role", "copy, derived, annotation (see above); metadata: run setting; plot: image path; rna: null without RNA data.",
  "orange", "source", "Tool file (copy) or input values (derived, annotation).",
  "orange", "rule", "For lists: which records ORANGE keeps. For derived values: how ORANGE computes them.",
  "orange", "json_type", "string, number, integer, boolean, date, enum, object or array.",
  "orange", "n_paths", "Number of JSON paths the row covers.",
  "orange", "filled_matched", "Whether the value is filled (not null, not empty) in the matched samples.",
  "orange", "filled_tumor_only", "Whether the value is filled in the tumor-only samples.",
  "orange", "notes", "Remarks, including errors in ORANGE."
)

TEXT_COLS <- c("description", "notes", "holds", "rule", "source")
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
    text <- which(names(df) %in% TEXT_COLS)
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
  list(readme = readme, files = files_out, fields = fields_out, orange = orange_out),
  workbook
)

log_msg("files: {nrow(files_out)} ({sum(files_out$role == 'primary')} primary)")
log_msg(
  "fields: {nrow(fields_out)} ({sum(fields_out$role == 'primary')} primary, ",
  "{sum(fields_out$role == 'copy')} copy, {sum(fields_out$role == 'derived')} derived, ",
  "{sum(fields_out$role == 'annotation')} annotation, {sum(fields_out$role == 'key')} key)"
)
log_msg("copies checked equal: {sum(str_starts(fields_out$check, 'equal'))}")
log_msg("orange: {nrow(orange_out)} rows covering {sum(orange_out$n_paths)} of {nrow(declared)} paths")
log_msg("Wrote {tables_dir} and {workbook}")
log_msg("{warn_state$n} warnings")
