#!/usr/bin/env Rscript
# Prism: cohort aggregation of a finished oncoanalyser run. Called by
# postOA.sh; see that script for the usage and the output layout.
#
#   Rscript postOA.R --out-dir OUT --res-dir RES
#
# Each per-sample table is stacked over all groups into one file in
# RES/data, with group_id and sample_id in front. The lines are copied as
# the tool wrote them: no row is left out and no value is changed or
# reformatted. A table is written only if every file it needs is there,
# has the same header as the others and has as many fields on each line
# as in its header; otherwise the script stops. So far: LINX.

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
})

# ---- arguments ---------------------------------------------------------

parse_args <- function(args) {
  opts <- list(out_dir = NULL, res_dir = NULL)
  i <- 1
  while (i <= length(args)) {
    key <- args[i] |> str_remove("^--") |> str_replace_all("-", "_")
    if (!key %in% names(opts) || i == length(args)) {
      stop(str_glue("Bad argument: {args[i]}"), call. = FALSE)
    }
    opts[[key]] <- args[i + 1]
    i <- i + 2
  }
  if (is.null(opts$out_dir) || is.null(opts$res_dir)) {
    stop("--out-dir and --res-dir are required", call. = FALSE)
  }
  opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
out_dir <- path_abs(opts$out_dir)
res_dir <- path_abs(opts$res_dir)
data_dir <- path(res_dir, "data")
dir_create(data_dir)

# ---- helpers -----------------------------------------------------------

log_msg <- function(...) {
  message(format(Sys.time(), "%H:%M:%S"), " ", str_glue(..., .envir = parent.frame()))
}

first_file <- function(dir, glob) {
  if (!dir_exists(dir)) {
    return(NA_character_)
  }
  found <- dir_ls(dir, glob = glob, type = "file")
  if (length(found) == 0) NA_character_ else as.character(found[1])
}

# Stack the per-sample files of one table into OUT_FILE: the header of the
# files with group_id and sample_id in front, then the lines of each file
# with its group_id and sample_id in front. FILES has the columns group_id,
# sample_id and file. Returns FILES with the number of rows each file
# added.
stack_table <- function(files, out_file) {
  lines <- map(files$file, read_lines)
  empty <- files$file[lengths(lines) == 0]
  if (length(empty) > 0) {
    stop("Empty file (no header): ", str_c(empty, collapse = ", "), call. = FALSE)
  }
  header <- lines[[1]][1]
  other_header <- files$file[map_chr(lines, 1) != header]
  if (length(other_header) > 0) {
    stop(
      "Header differs from ", files$file[1], ": ",
      str_c(other_header, collapse = ", "), call. = FALSE
    )
  }
  n_tabs <- str_count(header, "\t")
  body <- pmap(list(files$group_id, files$sample_id, files$file, lines), \(g, s, f, l) {
    rows <- l[-1]
    bad <- which(str_count(rows, "\t") != n_tabs)
    if (length(bad) > 0) {
      stop(
        str_glue("{f}: line {bad[1] + 1} has a different number of fields than the header"),
        call. = FALSE
      )
    }
    if (length(rows) == 0) character() else str_c(g, "\t", s, "\t", rows)
  })
  out_lines <- c(str_c("group_id\tsample_id\t", header), list_c(body, ptype = character()))
  write_lines(out_lines, out_file)
  files |> mutate(rows = lengths(body))
}

# ---- samples -----------------------------------------------------------

# One row per group with a PURPLE purity file. sample_id is the tumor
# sample, named by that file. normal_id is the other sample folder in
# bamtools/, if any. run_mode is PURPLE's (TUMOR_GERMLINE: with a matched
# normal; TUMOR: tumor-only).
group_dirs <- dir_ls(out_dir, type = "directory")
samples <- tibble(group_dir = as.character(group_dirs)) |>
  mutate(
    group_id = path_file(group_dir),
    purity_file = map_chr(
      group_dir, \(d) first_file(path(d, "purple"), "*.purple.purity.tsv")
    )
  ) |>
  filter(!is.na(purity_file)) |>
  arrange(group_id)
if (nrow(samples) == 0) {
  stop(str_glue("No groups with a PURPLE purity file in {out_dir}"), call. = FALSE)
}

other_samples <- function(group_dir, sample_id) {
  dir <- path(group_dir, "bamtools")
  if (!dir_exists(dir)) {
    return(NA_character_)
  }
  found <- dir_ls(dir, type = "directory") |> path_file() |> setdiff(sample_id)
  if (length(found) == 0) NA_character_ else str_c(found, collapse = ";")
}

read_version <- function(file) {
  if (!file_exists(file)) {
    return(NA_character_)
  }
  read_lines(file) |> str_subset("^version=") |> str_remove("^version=") |> first()
}

samples <- samples |>
  mutate(
    sample_id = purity_file |> path_file() |> str_remove("\\.purple\\.purity\\.tsv$"),
    normal_id = map2_chr(group_dir, sample_id, other_samples),
    run_mode = map_chr(purity_file, \(f) {
      read_tsv(f, col_types = cols(.default = col_character()), show_col_types = FALSE)$runMode
    }),
    matched_normal = run_mode == "TUMOR_GERMLINE",
    linx_version = map_chr(
      path(group_dir, "linx", "somatic_annotations", "linx.version"), read_version
    )
  )
log_msg("{nrow(samples)} groups in {out_dir}, {sum(samples$matched_normal)} with a matched normal")

mismatch <- samples |> filter(matched_normal == is.na(normal_id))
if (nrow(mismatch) > 0) {
  warning(
    "run_mode and the bamtools/ folders disagree on the normal for: ",
    str_c(mismatch$group_id, collapse = ", "), call. = FALSE
  )
}

# ---- LINX --------------------------------------------------------------

# The LINX tables. somatic: every group, from
# linx/somatic_annotations/<sample_id>.linx.<table>.tsv. germline: the
# groups with a matched normal, from
# linx/germline_annotations/<sample_id>.linx.germline.<table>.tsv.
# Not aggregated: neoepitope (input to NEO) and vis_* (input to the LINX
# plots). The descriptions were checked against LINX 2.3.1 (source code
# at tag linx-v2.3.1, task logs, files); see docs/ACTUAL_METHODS.md.
linx_tables <- tribble(
  ~mode, ~table, ~description,
  "somatic", "svs", str_c(
    "One row per SV LINX analyzed: the SVs of purple.sv.vcf.gz with FILTER ",
    "PASS or INFERRED; SVs with FILTER PON are left out. type INF: a ",
    "breakend PURPLE inferred from copy number (FILTER INFERRED). ",
    "vcfIdStart and vcfIdEnd: the VCF IDs of the two breakends; for some ",
    "SVs whose two breakends have the same position both hold the start ",
    "ID. Key svId; clusterId joins clusters."
  ),
  "somatic", "clusters", str_c(
    "One row per SV cluster: category, resolvedType, clusterCount (number ",
    "of SVs). The SVs of clusters with category ARTIFACT are kept in svs. ",
    "Key clusterId."
  ),
  "somatic", "links", str_c(
    "One row per link between two breakends of a chain. lowerSvId and ",
    "upperSvId join svs.svId; clusterId joins clusters."
  ),
  "somatic", "breakend", str_c(
    "One row per SV breakend and transcript LINX writes for it: the ",
    "canonical transcript of every gene whose span contains the breakend ",
    "or whose 5' end is up to 100 kb after it, plus driver panel and ",
    "fusion transcripts (regionType: where relative to the transcript, ",
    "most rows UPSTREAM; canonical: true for most rows), with disruptive ",
    "and reportedStatus. Key id; svId joins svs."
  ),
  "somatic", "fusion", str_c(
    "One row per fusion LINX found; reported: true for those LINX ",
    "reports. fivePrimeBreakendId and threePrimeBreakendId join breakend.id."
  ),
  "somatic", "drivers", str_c(
    "One row per driver gene event and the SV cluster LINX assigns it to; ",
    "clusterId -1: no cluster."
  ),
  "somatic", "driver.catalog", str_c(
    "The driver catalog: the rows of PURPLE's somatic driver catalog with ",
    "reportedStatus REPORTED, unchanged, plus LINX's gene disruption rows ",
    "(driver DISRUPTION)."
  ),
  "germline", "svs", str_c(
    "One row per SV of purple.sv.germline.vcf.gz with FILTER PASS or PON ",
    "(LINX keeps both in germline mode). Key svId; clusterId joins ",
    "germline clusters."
  ),
  "germline", "breakend", str_c(
    "One row per germline SV breakend and gene transcript LINX writes for ",
    "it (LINX loads canonical transcripts only in germline mode), with ",
    "disruptive and reportedStatus. Key id; svId joins germline svs."
  ),
  "germline", "disruption", str_c(
    "LINX's germline disruption list: one row per germline SV and gene, ",
    "with the SV's positions, type, FILTER, read support and cluster. ",
    "svId joins germline svs."
  ),
  "germline", "clusters", "One row per germline SV cluster. Key clusterId.",
  "germline", "links", str_c(
    "One row per link between two breakends of a germline chain. ",
    "lowerSvId and upperSvId join germline svs.svId."
  ),
  "germline", "driver.catalog", str_c(
    "LINX's germline driver catalog (gene disruptions by germline SVs), ",
    "with the columns of the somatic driver catalog."
  )
) |>
  mutate(
    name = str_c("linx_", mode, "_", str_replace_all(table, fixed("."), "_")),
    file = str_c("data/", name, ".tsv.gz"),
    infix = if_else(mode == "somatic", "linx", "linx.germline"),
    source = str_c(
      "out/<group_id>/linx/", mode, "_annotations/<sample_id>.", infix, ".", table, ".tsv"
    )
  )

sample_dirs <- samples |> select(group_id, sample_id, group_dir, matched_normal)
linx_files <- linx_tables |>
  select(mode, table, infix) |>
  cross_join(sample_dirs) |>
  filter(mode == "somatic" | matched_normal) |>
  mutate(
    file = path(
      group_dir, "linx", str_c(mode, "_annotations"),
      str_c(sample_id, ".", infix, ".", table, ".tsv")
    )
  )
missing <- linx_files |> filter(!file_exists(file))
if (nrow(missing) > 0) {
  message("Missing LINX files:\n", str_c(missing$file, collapse = "\n"))
  stop(nrow(missing), " LINX files missing (listed above)", call. = FALSE)
}

row_counts <- linx_tables |>
  pmap(\(mode, table, name, file, ...) {
    files <- linx_files |> filter(.data$mode == .env$mode, .data$table == .env$table)
    counts <- stack_table(files, path(res_dir, file))
    log_msg("LINX {mode} {table}: {nrow(files)} groups, {sum(counts$rows)} rows -> {file}")
    counts |> transmute(group_id, name = .env$name, rows)
  }) |>
  list_rbind()

# ---- write -------------------------------------------------------------

samples |>
  select(group_id, sample_id, normal_id, run_mode, linx_version) |>
  write_tsv(path(res_dir, "samples.tsv"), na = "")

samples |>
  select(group_id, sample_id, run_mode) |>
  left_join(
    row_counts |> pivot_wider(names_from = name, values_from = rows),
    by = "group_id"
  ) |>
  write_tsv(path(res_dir, "rows_per_sample.tsv"), na = "")

table_totals <- row_counts |>
  summarize(groups = n(), rows = sum(rows), .by = name)

manifest <- bind_rows(
  tibble(
    file = c("manifest.tsv", "samples.tsv", "rows_per_sample.tsv"),
    description = c(
      "This file: what each file holds and where it comes from.",
      str_c(
        "One row per group: group_id, sample_id (tumor), normal_id (empty ",
        "if tumor-only), run_mode (PURPLE: TUMOR_GERMLINE with a matched ",
        "normal, TUMOR tumor-only), linx_version."
      ),
      str_c(
        "The number of rows each group adds to each data table; empty where ",
        "the table does not apply (germline tables of tumor-only groups)."
      )
    )
  ),
  linx_tables |>
    left_join(table_totals, by = "name") |>
    transmute(
      file, tool = "LINX", mode, groups, rows, source,
      description = str_c(
        description, " Columns: group_id, sample_id, then those of the ",
        "LINX file (see the LINX README)."
      )
    )
) |>
  relocate(file, tool, mode, groups, rows, source, description)
write_tsv(manifest, path(res_dir, "manifest.tsv"), na = "")

log_msg("Done: {nrow(linx_tables)} LINX tables in {data_dir}")
