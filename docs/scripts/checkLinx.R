#!/usr/bin/env Rscript
# Check what LINX holds against its input, for every group of an
# oncoanalyser output folder. Evidence for the LINX rows of
# docs/ACTUAL_METHODS.md; rerun after a pipeline upgrade.
#
#   Rscript --vanilla Prism/docs/scripts/checkLinx.R [OUT_DIR]
#
# OUT_DIR defaults to out. Prints one table per check, one row per group;
# a check passes where its column is TRUE. It reads the PURPLE SV VCFs
# and the LINX tables only (small files; seconds).
#
# Somatic (linx/somatic_annotations):
#   svs_eq_log      svs.tsv rows = "procesed N SVs" in the task log
#   ids_eq_kept     the VCF IDs in svs.tsv = the IDs of the PURPLE SV VCF
#                   records with FILTER PASS or INFERRED (or empty)
#   records_ok      2 x two-breakend SVs + single-breakend SVs = those
#                   records
#   inf_eq_inferred svs.tsv type INF rows = FILTER INFERRED records
#   same_pos_svs, same_pos_id_repeated
#                   SVs whose two breakends have the same position, and
#                   those of them with vcfIdEnd = vcfIdStart (the reason
#                   ids_eq_kept can be FALSE while records_ok is TRUE)
#   keys_ok         breakend and links svIds in svs; svs and clusters
#                   clusterIds equal; clusterCount sums to the svs rows;
#                   fusion breakend IDs in breakend; drivers clusterIds in
#                   clusters or -1
#   catalog_ok      driver.catalog.tsv = the REPORTED rows of PURPLE's
#                   somatic driver catalog, unchanged, plus LINX rows
#                   (driver values in catalog_linx_drivers)
# Germline (linx/germline_annotations, groups with a matched normal):
#   g_svs_eq_log    germline svs.tsv rows = "loaded" = "procesed" in the
#                   task log
#   g_keys_ok       germline breakend svIds in svs; svs and clusters
#                   clusterIds equal

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
})
options(width = 250)

# SV VCF lines (assembled insert sequences) are longer than readr's
# default buffer for compressed files.
Sys.setenv(VROOM_CONNECTION_SIZE = 2^24)

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) > 0) args[1] else "out"

read_chr <- function(file, ...) {
  read_tsv(
    file, col_types = cols(.default = col_character()), quote = "",
    na = character(), show_col_types = FALSE, lazy = FALSE, trim_ws = FALSE, ...
  )
}

read_vcf_body <- function(file) {
  n_skip <- sum(str_starts(read_lines(file, n_max = 20000), "##"))
  read_chr(file, skip = n_skip, comment = "")
}

log_count <- function(file, pattern) {
  hit <- str_subset(read_lines(file), pattern)
  if (length(hit) == 0) NA_integer_ else as.integer(str_match(hit[1], pattern)[, 2])
}

check_somatic <- function(group, tumor, group_dir) {
  lx <- path(group_dir, "linx", "somatic_annotations", str_c(tumor, ".linx"))
  vcf <- read_vcf_body(path(group_dir, "purple", str_c(tumor, ".purple.sv.vcf.gz")))
  kept <- vcf |> filter(FILTER %in% c("PASS", "INFERRED", ".", ""))
  svs <- read_chr(str_c(lx, ".svs.tsv"))
  clusters <- read_chr(str_c(lx, ".clusters.tsv"))
  links <- read_chr(str_c(lx, ".links.tsv"))
  breakend <- read_chr(str_c(lx, ".breakend.tsv"))
  fusion <- read_chr(str_c(lx, ".fusion.tsv"))
  drivers <- read_chr(str_c(lx, ".drivers.tsv"))
  catalog <- read_chr(str_c(lx, ".driver.catalog.tsv"))
  purple_catalog <- read_chr(
    path(group_dir, "purple", str_c(tumor, ".purple.driver.catalog.somatic.tsv"))
  )
  log_file <- path(out_dir, "logs", group, str_c("linx_somatic.", group, ".command.log"))

  linx_ids <- c(svs$vcfIdStart, svs$vcfIdEnd) |> discard(\(x) x == "")
  two_sided <- svs$coordsEnd != ""
  same_pos <- two_sided &
    str_remove(svs$coordsStart, ":-?1$") == str_remove(svs$coordsEnd, ":-?1$")

  reported <- purple_catalog |> filter(reportedStatus == "REPORTED")
  key <- c("gene", "transcript", "driver")
  from_purple <- catalog |> semi_join(purple_catalog, by = key)
  linx_rows <- catalog |> anti_join(purple_catalog, by = key)
  catalog_ok <- nrow(from_purple) == nrow(reported) &&
    nrow(anti_join(from_purple, reported, by = names(reported))) == 0 &&
    nrow(anti_join(reported, from_purple, by = names(reported))) == 0

  keys_ok <- all(breakend$svId %in% svs$svId) &&
    all(c(links$lowerSvId, links$upperSvId) %in% svs$svId) &&
    setequal(svs$clusterId, clusters$clusterId) &&
    sum(as.integer(clusters$clusterCount)) == nrow(svs) &&
    all(c(fusion$fivePrimeBreakendId, fusion$threePrimeBreakendId) %in% breakend$id) &&
    all(drivers$clusterId %in% c(clusters$clusterId, "-1"))

  artifact <- clusters$clusterId[clusters$category == "ARTIFACT"]
  tibble(
    vcf_records = nrow(vcf),
    vcf_pass = sum(vcf$FILTER == "PASS"),
    vcf_inferred = sum(vcf$FILTER == "INFERRED"),
    vcf_pon = sum(vcf$FILTER == "PON"),
    vcf_other = nrow(vcf) - vcf_pass - vcf_inferred - vcf_pon,
    log_loaded = log_count(log_file, "loaded (\\d+) SV data records"),
    log_processed = log_count(log_file, "procesed (\\d+) SVs"),
    svs_rows = nrow(svs),
    svs_eq_log = svs_rows == log_processed,
    ids_eq_kept = setequal(linx_ids, kept$ID),
    records_ok = 2 * sum(two_sided) + sum(!two_sided) == nrow(kept),
    inf_eq_inferred = sum(svs$type == "INF") == vcf_inferred,
    same_pos_svs = sum(same_pos),
    same_pos_id_repeated = sum(same_pos & svs$vcfIdStart == svs$vcfIdEnd),
    artifact_svs = sum(svs$clusterId %in% artifact),
    keys_ok = keys_ok,
    catalog_ok = catalog_ok,
    catalog_linx_rows = nrow(linx_rows),
    catalog_linx_drivers = str_c(sort(unique(linx_rows$driver)), collapse = ";")
  )
}

check_germline <- function(group, tumor, group_dir) {
  gl_dir <- path(group_dir, "linx", "germline_annotations")
  if (!dir_exists(gl_dir)) {
    return(tibble(g_svs_rows = NA_integer_))
  }
  gl <- path(gl_dir, str_c(tumor, ".linx.germline"))
  svs <- read_chr(str_c(gl, ".svs.tsv"))
  clusters <- read_chr(str_c(gl, ".clusters.tsv"))
  breakend <- read_chr(str_c(gl, ".breakend.tsv"))
  log_file <- path(out_dir, "logs", group, str_c("linx_germline.", group, ".command.log"))
  loaded <- log_count(log_file, "loaded (\\d+) germline SV data records")
  processed <- log_count(log_file, "procesed (\\d+) SVs")
  tibble(
    g_svs_rows = nrow(svs),
    g_svs_eq_log = nrow(svs) == loaded && nrow(svs) == processed,
    g_keys_ok = all(breakend$svId %in% svs$svId) &&
      setequal(svs$clusterId, clusters$clusterId)
  )
}

groups <- tibble(group_dir = as.character(dir_ls(out_dir, type = "directory"))) |>
  mutate(group = path_file(group_dir)) |>
  filter(dir_exists(path(group_dir, "linx"))) |>
  mutate(
    tumor = map_chr(group_dir, \(d) {
      dir_ls(path(d, "purple"), glob = "*.purple.purity.tsv") |>
        path_file() |>
        str_remove("\\.purple\\.purity\\.tsv$") |>
        first()
    }),
    run_mode = map_chr(group_dir, \(d) {
      read_chr(dir_ls(path(d, "purple"), glob = "*.purple.purity.tsv")[1])$runMode
    })
  )

res <- groups |>
  mutate(
    somatic = pmap(list(group, tumor, group_dir), check_somatic),
    germline = pmap(list(group, tumor, group_dir), check_germline)
  ) |>
  select(-group_dir, -tumor) |>
  unnest(c(somatic, germline))

print(res |> select(group, run_mode, vcf_records:svs_rows), n = Inf)
print(res |> select(group, svs_eq_log:catalog_linx_drivers), n = Inf)
print(res |> select(group, starts_with("g_")), n = Inf)
