#!/usr/bin/env Rscript
# Prism: one table of the somatic SV events of all groups, with every field
# the pipeline wrote for them and technical evidence columns. Called by
# postOA.sh after postOA.R; see that script for the usage and the output.
#
#   Rscript svCohort.R --out-dir OUT --db-dir DB --track-dir TRACKS --hmf-dir HMF
#
# Rows: every record of PURPLE's purple.sv.vcf.gz of every group (FILTER
# PASS, PON, INFERRED; nothing is left out), one row per event. The two
# breakends of an SV (linked by MATEID) form one row; a single breakend or
# an inferred breakend forms one row. The breakend that comes first in the
# VCF is the start.
#
# Columns, each prefixed by its source (see the column dictionary written
# next to the table, which also gives each column's type):
#   no prefix         per event: VCF fields with the same value on both
#                     breakends of an SV (qual, filter, info_*, tumor_*,
#                     normal_*, esvee_bnd_*), checked
#   start_* / end_*   per breakend: the other VCF fields (info_*: INFO,
#                     tumor_* and normal_*: FORMAT of each sample), the
#                     ESVEE breakend.tsv columns that no VCF field repeats
#                     (esvee_bnd_*), LINX's breakend.tsv rows for that end
#                     (linx_bnd_*, one list element per transcript), region
#                     tracks (region_*), and the same breakend in other
#                     groups (xs_*) and in the group's own germline file
#                     (own_germline_*)
#   linx_*            the LINX svs.tsv row, its cluster, and lists of its
#                     links, fusions and cluster drivers
#   xs_pair_*         the same two-breakend event in other groups
# A cell holds one value, except where the source holds a list of any
# length; lists are joined with "|" (VCF lists with ",") in source order.
# Fixed pairs (AD, CIPOS, IHOMPOS, TPR) are split into two columns; columns
# that repeat another column are dropped; each of these steps is checked on
# every row (see "unpack" below). Derived columns are documented in the
# dictionary. No value is interpreted and no row is removed.

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
})

# GenomicRanges, IRanges and S4Vectors are called through their namespaces:
# attaching them would mask fs::path and dplyr verbs (rename, first, slice).

# Positions or intervals (1-based, both ends included) as GRanges.
as_ranges <- function(chrom, start, end = start) {
  GenomicRanges::GRanges(chrom, IRanges::IRanges(as.integer(start), as.integer(end)))
}

# SV VCF lines (assembled insert sequences) are longer than readr's
# default buffer for compressed files.
Sys.setenv(VROOM_CONNECTION_SIZE = 2^26)

SEP <- "|"
# Distances (bases) for "the same breakend" in other files: exact position,
# and within XS_WINDOW bases. The orientation must be the same.
XS_WINDOW <- 50

# ---- arguments ---------------------------------------------------------

parse_args <- function(args) {
  opts <- list(out_dir = NULL, db_dir = NULL, track_dir = NULL, hmf_dir = NULL)
  i <- 1
  while (i <= length(args)) {
    key <- args[i] |> str_remove("^--") |> str_replace_all("-", "_")
    if (!key %in% names(opts) || i == length(args)) {
      stop(str_glue("Bad argument: {args[i]}"), call. = FALSE)
    }
    opts[[key]] <- args[i + 1]
    i <- i + 2
  }
  missing <- names(opts)[map_lgl(opts, is.null)]
  if (length(missing) > 0) {
    stop("Missing arguments: ", str_c("--", missing, collapse = " "), call. = FALSE)
  }
  opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
out_dir <- path_abs(opts$out_dir)
db_dir <- path_abs(opts$db_dir)
track_dir <- path_abs(opts$track_dir)
hmf_dir <- path_abs(opts$hmf_dir)
out_file <- path(db_dir, "sv_cohort_filtered.tsv.gz")
dict_file <- path(db_dir, "sv_cohort_filtered.columns.tsv")

# ---- helpers -----------------------------------------------------------

log_msg <- function(...) {
  message(format(Sys.time(), "%H:%M:%S"), " ", str_glue(..., .envir = parent.frame()))
}

read_chr <- function(file, ...) {
  read_tsv(
    file, col_types = cols(.default = col_character()), quote = "",
    na = character(), show_col_types = FALSE, lazy = FALSE, trim_ws = FALSE, ...
  )
}

# A VCF as text: the header lines and the body (all columns character).
read_vcf <- function(file) {
  header <- character()
  con <- gzfile(file, "r")
  repeat {
    line <- readLines(con, n = 1)
    if (length(line) == 0 || !str_starts(line, "##")) break
    header <- c(header, line)
  }
  close(con)
  body <- read_chr(file, skip = length(header)) |> rename(CHROM = `#CHROM`)
  list(header = header, body = body)
}

# INFO or FORMAT definitions (Number, Type, Description) from VCF header lines.
header_fields <- function(header, kind) {
  m <- str_match(
    header, str_glue('^##{kind}=<ID=([^,]+),Number=([^,]+),Type=([^,]+),.*Description="([^"]*)"')
  )
  tibble(key = m[, 2], number = m[, 3], type = m[, 4], description = m[, 5]) |> filter(!is.na(key))
}

# Breakend orientation from the VCF ALT (breakend notation): 1 if the ALT
# starts with the reference base (sequence kept on the left), else -1. The
# same as ESVEE's and LINX's orientation (checked against both).
alt_orientation <- function(alt) {
  if_else(str_detect(alt, "^[ACGTNacgtn]"), 1L, -1L)
}

# INFO strings to one column per key (info_<KEY>); flags get "true".
info_wide <- function(info) {
  tibble(row = seq_along(info), kv = str_split(info, ";")) |>
    unnest(kv) |>
    filter(kv != "" & kv != ".") |>
    mutate(
      key = str_remove(kv, "=.*$"),
      value = if_else(str_detect(kv, "="), str_remove(kv, "^[^=]*="), "true")
    ) |>
    select(row, key, value) |>
    pivot_wider(names_from = key, values_from = value, names_prefix = "info_") |>
    right_join(tibble(row = seq_along(info)), by = "row") |>
    arrange(row) |>
    select(-row)
}

# FORMAT keys and one sample's values to one column per key (<prefix><KEY>).
format_wide <- function(format, values, prefix) {
  keys <- str_split(format, ":")
  vals <- str_split(values, ":")
  vals <- map2(keys, vals, \(k, v) c(v, rep(NA_character_, max(0, length(k) - length(v))))[seq_along(k)])
  tibble(row = seq_along(format), key = keys, value = vals) |>
    unnest(c(key, value)) |>
    pivot_wider(names_from = key, values_from = value, names_prefix = prefix) |>
    right_join(tibble(row = seq_along(format)), by = "row") |>
    arrange(row) |>
    select(-row)
}

# Join the rows of DF per BY into one row: each column's values joined with
# SEP in source order, and <prefix>n, the number of rows.
collapse_rows <- function(df, by, prefix) {
  value_cols <- setdiff(names(df), by)
  hits <- map_lgl(value_cols, \(col) any(str_detect(df[[col]], fixed(SEP)), na.rm = TRUE))
  if (any(hits)) {
    stop(
      "Separator '", SEP, "' found in values of: ", str_c(value_cols[hits], collapse = ", "),
      call. = FALSE
    )
  }
  df |>
    summarize(
      n = as.character(n()),
      across(all_of(value_cols), \(x) str_c(replace_na(x, ""), collapse = SEP)),
      .by = all_of(by)
    ) |>
    rename_with(\(x) str_c(prefix, x), -all_of(by))
}

coords_key <- function(chrom, pos, orient) str_c(chrom, ":", pos, ":", orient)

# ---- samples -----------------------------------------------------------

samples <- read_chr(path(db_dir, "samples.tsv"))
log_msg("{nrow(samples)} groups from {path(db_dir, 'samples.tsv')}")

linx <- c("svs", "clusters", "links", "breakend", "fusion", "drivers") |>
  set_names() |>
  map(\(t) read_chr(path(db_dir, "data", str_c("linx_somatic_", t, ".tsv.gz"))))

# ---- per-group events --------------------------------------------------

# Returns the group's events and the INFO and FORMAT descriptions of its
# VCF header.
build_group <- function(group_id, sample_id, normal_id, run_mode, ...) {
  g_dir <- path(out_dir, group_id)
  vcf <- read_vcf(path(g_dir, "purple", str_c(sample_id, ".purple.sv.vcf.gz")))

  body <- vcf$body
  sample_cols <- setdiff(names(body), c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"))
  if (!sample_id %in% sample_cols) {
    stop(str_glue("{group_id}: no column {sample_id} in the PURPLE SV VCF"), call. = FALSE)
  }
  other_cols <- setdiff(sample_cols, sample_id)
  if (length(other_cols) > 1) {
    stop(str_glue("{group_id}: more than two samples in the PURPLE SV VCF"), call. = FALSE)
  }

  bnd <- body |>
    transmute(
      id = ID, chrom = CHROM, pos = POS, orient = as.character(alt_orientation(ALT)),
      ref = REF, alt = ALT, qual = QUAL, filter = FILTER, vcf_order = as.character(row_number())
    ) |>
    bind_cols(
      info_wide(body$INFO),
      format_wide(body$FORMAT, body[[sample_id]], "tumor_"),
      if (length(other_cols) == 1) format_wide(body$FORMAT, body[[other_cols]], "normal_")
    )
  if (anyDuplicated(bnd$id) > 0) stop(str_glue("{group_id}: duplicate VCF IDs"), call. = FALSE)
  # "." is the VCF's missing value: stored as empty.
  bnd <- bnd |>
    mutate(across(c(qual, filter, starts_with(c("info_", "tumor_", "normal_"))), \(x) na_if(x, ".")))

  # ESVEE breakend.tsv: every column, for the records in the PURPLE file
  # (inferred breakends are PURPLE's and have none).
  esvee <- read_chr(path(g_dir, "esvee", str_c(sample_id, ".esvee.breakend.tsv"))) |>
    filter(Id %in% bnd$id) |>
    rename_with(\(x) str_c("esvee_bnd_", x), -Id)
  bnd <- bnd |> left_join(esvee, by = c(id = "Id"))
  esvee_bad <- bnd |>
    filter(filter != "INFERRED") |>
    filter(
      is.na(esvee_bnd_Chromosome) | esvee_bnd_Chromosome != chrom | esvee_bnd_Position != pos |
        esvee_bnd_Orientation != orient
    )
  if (nrow(esvee_bad) > 0) {
    stop(
      str_glue("{group_id}: {nrow(esvee_bad)} records without a matching ESVEE breakend.tsv row"),
      call. = FALSE
    )
  }

  # Pair the breakends: MATEID both ways; the first in the VCF is the start.
  mate <- bnd$info_MATEID
  if (is.null(mate)) mate <- rep(NA_character_, nrow(bnd))
  has_mate <- !is.na(mate) & mate %in% bnd$id
  mate_of_mate <- mate[match(mate, bnd$id)]
  if (any(has_mate & mate_of_mate != bnd$id)) {
    stop(str_glue("{group_id}: MATEID links that are not reciprocal"), call. = FALSE)
  }
  order <- as.integer(bnd$vcf_order)
  is_start <- !has_mate | order < order[match(mate, bnd$id)]
  events <- tibble(start_id = bnd$id[is_start], end_id = if_else(has_mate, mate, NA_character_)[is_start])
  if (2 * sum(!is.na(events$end_id)) + sum(is.na(events$end_id)) != nrow(bnd)) {
    stop(str_glue("{group_id}: breakends not accounted for by the events"), call. = FALSE)
  }

  start_cols <- bnd |> rename_with(\(x) str_c("start_", x))
  end_cols <- bnd |> rename_with(\(x) str_c("end_", x))
  events <- events |>
    left_join(start_cols, by = "start_id") |>
    left_join(end_cols, by = "end_id") |>
    mutate(
      group_id = group_id, sample_id = sample_id, normal_id = normal_id, run_mode = run_mode,
      vcf_normal_sample = if (length(other_cols) == 1) other_cols else NA_character_,
      .before = 1
    )
  log_msg("{group_id}: {nrow(bnd)} VCF records, {nrow(events)} events")
  list(
    events = events,
    info = header_fields(vcf$header, "INFO"),
    format = header_fields(vcf$header, "FORMAT")
  )
}

groups <- samples |> pmap(build_group)
events <- groups |> map("events") |> list_rbind()
header_defs <- function(kind) {
  defs <- groups |> map(kind) |> list_rbind() |> distinct(key, number, type, .keep_all = TRUE)
  differ <- defs$key[duplicated(defs$key)]
  if (length(differ) > 0) {
    stop(str_glue("{kind} keys with Number or Type differing between groups: {str_c(differ, collapse = ', ')}"),
         call. = FALSE)
  }
  defs
}
info_desc <- header_defs("info")
format_desc <- header_defs("format")
rm(groups)

events <- events |>
  mutate(
    event_key = str_c(sample_id, ":", start_id),
    n_breakends = if_else(is.na(end_id), 1L, 2L),
    same_chrom = if_else(is.na(end_id), NA, start_chrom == end_chrom),
    length = if_else(
      same_chrom %in% TRUE, abs(as.integer(end_pos) - as.integer(start_pos)), NA_integer_
    ),
    start_coords = coords_key(start_chrom, start_pos, start_orient),
    end_coords = if_else(is.na(end_id), NA_character_, coords_key(end_chrom, end_pos, end_orient)),
    .after = run_mode
  )

# ---- LINX --------------------------------------------------------------

lx_svs <- linx$svs |>
  rename_with(\(x) str_c("linx_", x), -c(group_id, sample_id)) |>
  mutate(start_id = linx_vcfIdStart)
matched <- events |>
  select(group_id, sample_id, start_id, end_id) |>
  inner_join(lx_svs, by = c("group_id", "sample_id", "start_id"))
bad_end <- matched |>
  filter(!(linx_vcfIdEnd == coalesce(end_id, "") | linx_vcfIdEnd == start_id))
if (nrow(bad_end) > 0) stop("LINX vcfIdEnd does not match the event's end", call. = FALSE)
if (nrow(matched) != nrow(linx$svs)) {
  stop(
    str_glue("{nrow(linx$svs) - nrow(matched)} LINX SVs do not match a PURPLE event"),
    call. = FALSE
  )
}

clusters <- linx$clusters |> rename_with(\(x) str_c("linx_cluster_", x), -c(group_id, sample_id))
events <- events |>
  left_join(
    matched |> select(-end_id),
    by = c("group_id", "sample_id", "start_id")
  ) |>
  left_join(clusters, by = c("group_id", "sample_id", linx_clusterId = "linx_cluster_clusterId")) |>
  mutate(in_linx = !is.na(linx_svId), .after = n_breakends)

# LINX start and end are the VCF's start and end except for some SVs whose
# two breakends share a position (see docs/ACTUAL_METHODS.md); the LINX
# columns are kept as LINX wrote them.
events <- events |>
  mutate(
    linx_start_is_vcf_start = if_else(in_linx, linx_coordsStart == start_coords, NA),
    .after = in_linx
  )
flipped <- events |> filter(linx_start_is_vcf_start %in% FALSE)
if (any(flipped$start_pos != coalesce(flipped$end_pos, ""))) {
  stop("LINX start differs from the VCF start for an SV whose breakends differ in position",
       call. = FALSE)
}
log_msg("LINX: {sum(events$in_linx)} of {nrow(events)} events; {nrow(flipped)} with LINX start = VCF end")

# LINX breakend rows (one per transcript) to an event end by svId and
# isStart: LINX's start is the VCF start unless linx_start_is_vcf_start is
# FALSE. Not by coordinates: an SV can have both ends at the same
# coordinates, and LINX writes isStart=false rows for single breakends
# (at a location that is not a VCF record, no vcfId), which go to the end.
lx_bnd <- linx$breakend
bnd_cols <- setdiff(names(lx_bnd), c("group_id", "sample_id", "svId"))
bnd_side <- lx_bnd |>
  inner_join(
    events |>
      filter(in_linx) |>
      select(group_id, sample_id, svId = linx_svId, linx_start_is_vcf_start),
    by = c("group_id", "sample_id", "svId")
  ) |>
  mutate(side = if_else((isStart == "true") == linx_start_is_vcf_start, "start", "end"))
if (nrow(bnd_side) != nrow(lx_bnd) || any(!lx_bnd$isStart %in% c("true", "false"))) {
  stop("LINX breakend rows without an event, or with isStart not true/false", call. = FALSE)
}
for (s in c("start", "end")) {
  side_bnd <- bnd_side |>
    filter(side == s) |>
    select(group_id, sample_id, svId, all_of(bnd_cols)) |>
    collapse_rows(c("group_id", "sample_id", "svId"), str_c(s, "_linx_bnd_"))
  events <- events |>
    left_join(side_bnd, by = c("group_id", "sample_id", linx_svId = "svId"))
}
n_bnd <- sum(as.integer(events$start_linx_bnd_n), as.integer(events$end_linx_bnd_n), na.rm = TRUE)
if (n_bnd != nrow(lx_bnd)) {
  stop(str_glue("{nrow(lx_bnd) - n_bnd} LINX breakend rows not placed on an event end"), call. = FALSE)
}
n_sgl_end <- sum(is.na(events$end_id) & !is.na(events$end_linx_bnd_n))
log_msg("LINX breakend: {nrow(lx_bnd)} rows placed; {n_sgl_end} single-breakend events with an end row")

# Fusions, per SV: each fusion is listed for the SV(s) of its 5' and 3'
# breakends, with the side.
be_sv <- lx_bnd |> select(group_id, sample_id, id, svId)
be_5p <- be_sv |> rename(sv_5p = svId)
be_3p <- be_sv |> rename(sv_3p = svId)
fusions <- linx$fusion |>
  left_join(be_5p, by = c("group_id", "sample_id", fivePrimeBreakendId = "id")) |>
  left_join(be_3p, by = c("group_id", "sample_id", threePrimeBreakendId = "id"))
fusion_cols <- setdiff(names(linx$fusion), c("group_id", "sample_id"))
fusion_by_sv <- bind_rows(
  fusions |> transmute(group_id, sample_id, svId = sv_5p, side = "5p", across(all_of(fusion_cols))),
  fusions |> transmute(group_id, sample_id, svId = sv_3p, side = "3p", across(all_of(fusion_cols)))
) |>
  collapse_rows(c("group_id", "sample_id", "svId"), "linx_fusion_")
links_by_sv <- bind_rows(
  linx$links |> mutate(svId = lowerSvId),
  linx$links |> filter(upperSvId != lowerSvId) |> mutate(svId = upperSvId)
) |>
  collapse_rows(c("group_id", "sample_id", "svId"), "linx_link_")
drivers_by_cluster <- linx$drivers |>
  collapse_rows(c("group_id", "sample_id", "clusterId"), "linx_cluster_driver_")
events <- events |>
  left_join(fusion_by_sv, by = c("group_id", "sample_id", linx_svId = "svId")) |>
  left_join(links_by_sv, by = c("group_id", "sample_id", linx_svId = "svId")) |>
  left_join(drivers_by_cluster, by = c("group_id", "sample_id", linx_clusterId = "clusterId"))

# ---- breakend points ---------------------------------------------------

start_points <- events |>
  transmute(
    row = row_number(), side = "start", group_id, chrom = start_chrom,
    pos = as.integer(start_pos), orient = start_orient
  )
end_points <- events |>
  mutate(row = row_number()) |>
  filter(!is.na(end_id)) |>
  transmute(row, side = "end", group_id, chrom = end_chrom, pos = as.integer(end_pos), orient = end_orient)
points <- bind_rows(start_points, end_points)
points_gr <- as_ranges(points$chrom, points$pos)

# Add per-breakend columns (tibble with columns row, side and values) to
# events as start_<col> and end_<col>.
add_point_cols <- function(events, vals) {
  for (s in c("start", "end")) {
    v <- vals |> filter(side == s) |> select(-side)
    names(v)[-1] <- str_c(s, "_", names(v)[-1])
    events <- events |>
      mutate(row = row_number()) |>
      left_join(v, by = "row") |>
      select(-row)
  }
  events
}

# ---- region tracks -----------------------------------------------------

log_msg("region tracks from {track_dir}")
tracks <- read_chr(path(track_dir, "tracks.tsv"))
bed_tracks <- tracks |> filter(str_starts(name, "giab_") | name == "encode_blacklist")

# Labels of the intervals of GR that contain each point, unique, joined
# with SEP; NA where none does.
point_labels <- function(gr, labels) {
  hits <- IRanges::findOverlaps(points_gr, gr)
  tibble(q = S4Vectors::queryHits(hits), label = labels[S4Vectors::subjectHits(hits)]) |>
    summarize(label = str_c(unique(label), collapse = SEP), .by = q) |>
    right_join(tibble(q = seq_along(points_gr)), by = "q") |>
    arrange(q) |>
    pull(label)
}

# BED: 0-based start, end exclusive; a breakend at POS p is in it when
# start < p <= end, so the interval is [start + 1, end] 1-based.
read_bed <- function(file) {
  read_tsv(file, col_names = FALSE, col_types = cols(.default = "c"), show_col_types = FALSE)
}

region <- points |> select(row, side)
for (i in seq_len(nrow(bed_tracks))) {
  name <- bed_tracks$name[i]
  b <- read_bed(path(track_dir, bed_tracks$file[i]))
  gr <- as_ranges(b$X1, as.integer(b$X2) + 1L, b$X3)
  region[[str_c("region_", name)]] <- IRanges::overlapsAny(points_gr, gr)
  if (name == "encode_blacklist") {
    region$region_encode_blacklist_reason <- point_labels(gr, b$X4)
  }
}

# Hartwig resource regions: 1-based, both ends included.
# Read as Hartwig's own loader does (hmf-common ChrBaseRegion
# .loadChrBaseRegions): a .bed file is 0-based (start + 1), any other
# file is 1-based as written.
read_hmf_regions <- function(file) {
  b <- read_chr(file)
  shift <- if (str_detect(file, "\\.bed(\\.gz)?$")) 1L else 0L
  as_ranges(b$Chromosome, as.integer(b$PosStart) + shift, b$PosEnd)
}
region$region_hmf_sv_prep_blacklist <- IRanges::overlapsAny(
  points_gr, read_hmf_regions(path(hmf_dir, "dna", "sv", "sv_prep_blacklist.38.bed"))
)
region$region_hmf_unmap <- IRanges::overlapsAny(
  points_gr, read_hmf_regions(path(hmf_dir, "common", "unmap_regions.38.tsv"))
)

# RepeatMasker: name, class and family of the repeats at the breakend.
rmsk <- read_tsv(
  path(track_dir, tracks$file[tracks$name == "ucsc_rmsk"]), col_names = FALSE,
  col_select = c(6, 7, 8, 11, 12, 13), col_types = cols(.default = "c"), show_col_types = FALSE
)
rmsk_gr <- as_ranges(rmsk$X6, as.integer(rmsk$X7) + 1L, rmsk$X8)
region$region_rmsk_name <- point_labels(rmsk_gr, rmsk$X11)
region$region_rmsk_class <- point_labels(rmsk_gr, rmsk$X12)
region$region_rmsk_family <- point_labels(rmsk_gr, rmsk$X13)
rm(rmsk, rmsk_gr)

# Distance to the nearest centromere model on the same chromosome: bases
# strictly between (0 inside or at the next base; NA with no model).
cen <- read_bed(path(track_dir, tracks$file[tracks$name == "ucsc_centromeres"]))
cen_gr <- as_ranges(cen$X2, as.integer(cen$X3) + 1L, cen$X4)
near <- IRanges::distanceToNearest(points_gr, cen_gr)
region$region_centromere_distance <- NA_integer_
region$region_centromere_distance[S4Vectors::queryHits(near)] <- S4Vectors::mcols(near)$distance

events <- add_point_cols(events, region)

# ---- the same breakends in other files ---------------------------------

log_msg("cross-sample sets")
read_breakends <- function(file, group_id) {
  b <- read_vcf(file)$body
  tibble(
    group_id = group_id, id = b$ID, chrom = b$CHROM, pos = as.integer(b$POS),
    orient = as.character(alt_orientation(b$ALT)), filter = b$FILTER,
    mate = str_match(b$INFO, "(?:^|;)MATEID=([^;]+)")[, 2]
  )
}
sets <- list(
  # PURPLE somatic records, without the breakends inferred from copy number.
  somatic = samples |>
    pmap(\(group_id, sample_id, ...) {
      read_breakends(path(out_dir, group_id, "purple", str_c(sample_id, ".purple.sv.vcf.gz")), group_id)
    }) |>
    list_rbind() |>
    filter(filter != "INFERRED"),
  # PURPLE germline SV records (groups with a matched normal).
  germline = samples |>
    filter(!is.na(normal_id) & normal_id != "") |>
    pmap(\(group_id, sample_id, ...) {
      read_breakends(path(out_dir, group_id, "purple", str_c(sample_id, ".purple.sv.germline.vcf.gz")), group_id)
    }) |>
    list_rbind(),
  # Every breakend ESVEE wrote, any FILTER.
  unfiltered = samples |>
    pmap(\(group_id, sample_id, ...) {
      read_breakends(path(out_dir, group_id, "esvee", str_c(sample_id, ".esvee.unfiltered.vcf.gz")), group_id)
    }) |>
    list_rbind()
)
walk2(names(sets), sets, \(n, s) log_msg("  {n}: {nrow(s)} breakends in {n_distinct(s$group_id)} groups"))

# Hits of each point against a set: same chromosome, same orientation,
# position within XS_WINDOW.
point_hits <- function(set) {
  set_gr <- as_ranges(set$chrom, set$pos)
  h <- IRanges::findOverlaps(points_gr, set_gr, maxgap = XS_WINDOW)
  tibble(p = S4Vectors::queryHits(h), s = S4Vectors::subjectHits(h)) |>
    mutate(dist = abs(points$pos[p] - set$pos[s])) |>
    filter(dist <= XS_WINDOW, points$orient[p] == set$orient[s]) |>
    mutate(point_group = points$group_id[p], set_group = set$group_id[s])
}

# Per point: the number of distinct groups among hits D, and their list.
count_groups <- function(d) {
  d |>
    summarize(n = n_distinct(set_group), list = str_c(sort(unique(set_group)), collapse = SEP), .by = p) |>
    right_join(tibble(p = seq_len(nrow(points))), by = "p") |>
    arrange(p) |>
    mutate(n = replace_na(n, 0L))
}

hits <- map(sets, \(set) point_hits(set) |> mutate(set_id = set$id[s], set_mate = set$mate[s]))
xs <- points |> select(row, side)
for (n in names(sets)) {
  other <- hits[[n]] |> filter(set_group != point_group)
  exact <- count_groups(other |> filter(dist == 0))
  near <- count_groups(other)
  xs[[str_c("xs_", n, "_groups_exact")]] <- exact$n
  xs[[str_c("xs_", n, "_groups_w", XS_WINDOW)]] <- near$n
  xs[[str_c("xs_", n, "_group_list_w", XS_WINDOW)]] <- near$list
}
# The group's own germline file (matched groups; empty for tumor-only).
own <- hits$germline |>
  filter(set_group == point_group) |>
  summarize(exact = sum(dist == 0), near = n(), .by = p) |>
  right_join(tibble(p = seq_len(nrow(points))), by = "p") |>
  arrange(p)
has_germline <- points$group_id %in% (samples |> filter(!is.na(normal_id) & normal_id != "") |> pull(group_id))
xs$own_germline_exact <- if_else(has_germline, replace_na(own$exact, 0L), NA_integer_)
xs[[str_c("own_germline_w", XS_WINDOW)]] <- if_else(has_germline, replace_na(own$near, 0L), NA_integer_)
events <- add_point_cols(events, xs)

# Two-breakend events: the same pair (both breakends within XS_WINDOW, same
# orientations) in another group.
pair_rows <- which(!is.na(events$end_id))
for (n in names(sets)) {
  h <- hits[[n]] |> filter(set_group != point_group)
  hs <- h |>
    filter(points$side[p] == "start") |>
    transmute(row = points$row[p], set_group, id1 = set_id, mate1 = set_mate)
  he <- h |>
    filter(points$side[p] == "end") |>
    transmute(row = points$row[p], set_group, id2 = set_id)
  pairs <- inner_join(hs, he, by = c("row", "set_group"), relationship = "many-to-many") |>
    filter(!is.na(mate1), mate1 == id2) |>
    summarize(n = n_distinct(set_group), list = str_c(sort(unique(set_group)), collapse = SEP), .by = row)
  col_n <- str_c("xs_pair_", n, "_groups_w", XS_WINDOW)
  col_l <- str_c("xs_pair_", n, "_group_list_w", XS_WINDOW)
  events[[col_n]] <- NA_integer_
  events[[col_n]][pair_rows] <- 0L
  events[[col_n]][pairs$row] <- pairs$n
  events[[col_l]] <- NA_character_
  events[[col_l]][pairs$row] <- pairs$list
}

# ---- unpack, deduplicate, event-level columns ----------------------------
#
# The table is for analysis: a cell holds one value, except where the source
# holds a list of any length (LINX transcripts, fusions, links and drivers,
# ASMLNKS, RepeatMasker labels, group lists). Each step below checks its
# assumption on every row and stops the run, naming the column, if a row
# breaks it. docs/scripts/checkSvCohort.R rebuilds every source value from
# the result.

log_msg("unpack, deduplicate, event-level columns")
OTHER_SIDE <- c(start = "end", end = "start")
n_values <- function(x, sep) str_count(x, fixed(sep)) + 1L
nth_value <- function(x, sep, i) str_split_i(x, fixed(sep), i)
blank_na <- function(x) if_else(x == "", NA_character_, as.character(x))
has_col <- function(col) col %in% names(events)
stop_unless <- function(ok, what, col, values) {
  bad <- which(!coalesce(ok, FALSE))
  if (length(bad) > 0) {
    stop(str_glue("{what}: {col}, {length(bad)} rows, e.g. '{values[bad[1]]}'"), call. = FALSE)
  }
}
pair_rows <- events$n_breakends == 2L

# 1. Columns that repeat another column: dropped after checking that every
# value they hold equals the kept column ({o}: the other breakend). Lists
# are compared after mapping SEP_FROM to ",".
DUPLICATES <- tribble(
  ~drop, ~keep, ~sep_from,
  "esvee_bnd_AssemblyId", "info_ASMID", NA,
  "esvee_bnd_AssemblyInfo", "info_ASMINFO", NA,
  "esvee_bnd_Type", "info_SVTYPE", NA,
  "esvee_bnd_Chromosome", "chrom", NA,
  "esvee_bnd_Position", "pos", NA,
  "esvee_bnd_Orientation", "orient", NA,
  "info_BEOR", "orient", NA,
  "esvee_bnd_MateId", "{o}_id", NA,
  "info_MATEID", "{o}_id", NA,
  "esvee_bnd_MateChr", "{o}_chrom", NA,
  "esvee_bnd_MatePos", "{o}_pos", NA,
  "esvee_bnd_MateOrient", "{o}_orient", NA,
  "esvee_bnd_ConfidenceInterval", "info_CIPOS", NA,
  "esvee_bnd_InexactOffset", "info_IHOMPOS", NA,
  "esvee_bnd_Homology", "info_HOMSEQ", NA,
  "esvee_bnd_Qual", "qual", NA,
  "esvee_bnd_SplitFragments", "info_SF", NA,
  "esvee_bnd_DiscFragments", "info_DF", NA,
  "esvee_bnd_SequenceLength", "info_ASMLEN", NA,
  "esvee_bnd_SegmentIndex", "info_ASMSEG", NA,
  "esvee_bnd_SequenceIndex", "info_BEAPOS", NA,
  "esvee_bnd_AlignedBases", "info_SEGALEN", NA,
  "esvee_bnd_MapQual", "info_SEGMAPQ", NA,
  "esvee_bnd_Score", "info_SEGSCO", NA,
  "esvee_bnd_AdjAlignedBases", "info_SEGRL", NA,
  "esvee_bnd_FacingBreakendIds", "info_ASMLNKS", ";",
  "esvee_bnd_AltAlignments", "info_INSALN", NA,
  "esvee_bnd_UniqueFragPos", "info_UFP", NA,
  "esvee_bnd_ThreePrimeRange", "info_TPR", ";"
)
dropped <- character()
for (s in c("start", "end")) {
  for (i in seq_len(nrow(DUPLICATES))) {
    d_col <- str_c(s, "_", DUPLICATES$drop[i])
    if (!has_col(d_col)) next
    k_col <- if (str_starts(DUPLICATES$keep[i], fixed("{o}"))) {
      str_replace(DUPLICATES$keep[i], fixed("{o}"), OTHER_SIDE[[s]])
    } else {
      str_c(s, "_", DUPLICATES$keep[i])
    }
    d <- blank_na(events[[d_col]])
    if (!is.na(DUPLICATES$sep_from[i])) d <- str_replace_all(d, fixed(DUPLICATES$sep_from[i]), ",")
    k <- if (has_col(k_col)) blank_na(events[[k_col]]) else rep(NA_character_, nrow(events))
    stop_unless(is.na(d) | d == k, str_glue("values differ from {k_col}"), d_col, d)
    events[[d_col]] <- NULL
    dropped <- c(dropped, d_col)
  }
  # ESVEE StrandBias: one value per sample (tumor, then normal), FORMAT SB
  # rounded to two decimals (format "%.2f" in ESVEE's BreakendWriter).
  sb_col <- str_c(s, "_esvee_bnd_StrandBias")
  if (has_col(sb_col)) {
    sb <- blank_na(events[[sb_col]])
    no_values <- rep(NA_character_, nrow(events))
    t_sb <- if (has_col(str_c(s, "_tumor_SB"))) blank_na(events[[str_c(s, "_tumor_SB")]]) else no_values
    n_sb <- if (has_col(str_c(s, "_normal_SB"))) blank_na(events[[str_c(s, "_normal_SB")]]) else no_values
    expect <- if_else(is.na(n_sb), t_sb, str_c(t_sb, ";", n_sb))
    stop_unless(
      map2_lgl(sb, expect, \(a, b) {
        if (is.na(a)) return(TRUE)
        x <- as.numeric(str_split_1(a, ";"))
        y <- as.numeric(str_split_1(coalesce(b, ""), ";"))
        length(x) == length(y) && !anyNA(c(x, y)) && all(abs(x - y) <= 0.005 + 1e-9)
      }),
      "values are not tumor_SB;normal_SB to two decimals", sb_col, sb
    )
    events[[sb_col]] <- NULL
    dropped <- c(dropped, sb_col)
  }
}
vns <- coalesce(events$vcf_normal_sample, "")
stop_unless(vns == coalesce(events$normal_id, ""), "values differ from normal_id", "vcf_normal_sample", vns)
events$vcf_normal_sample <- NULL
log_msg("  {length(dropped) + 1} duplicate columns dropped")

# 2. PURPLE_AF, PURPLE_CN, PURPLE_CN_CHANGE: PURPLE writes [start, end] on the
# start record and [end, start] on the end record (SomaticSvCache at
# purple-v4.4); keep each record's first value, its own.
for (k in c("PURPLE_AF", "PURPLE_CN", "PURPLE_CN_CHANGE")) {
  s_col <- str_c("start_info_", k)
  e_col <- str_c("end_info_", k)
  if (!has_col(s_col)) next
  sv <- events[[s_col]]
  ev <- if (has_col(e_col)) events[[e_col]] else rep(NA_character_, nrow(events))
  ok <- if_else(
    pair_rows,
    (is.na(sv) & is.na(ev)) |
      (n_values(sv, ",") %in% 2L & n_values(ev, ",") %in% 2L &
         nth_value(sv, ",", 2) == nth_value(ev, ",", 1) & nth_value(ev, ",", 2) == nth_value(sv, ",", 1)),
    is.na(ev) & (is.na(sv) | n_values(sv, ",") %in% 1L)
  )
  stop_unless(ok, "not [own, mate] on the two records", str_c("info_", k), coalesce(sv, ""))
  events[[s_col]] <- nth_value(sv, ",", 1)
  if (has_col(e_col)) events[[e_col]] <- nth_value(ev, ",", 1)
}

# 3. Fixed pairs: split into two columns.
PAIRS <- tribble(
  ~field, ~part1, ~part2,
  "info_CIPOS", "lower", "upper",
  "info_IHOMPOS", "lower", "upper",
  "info_TPR", "five_prime_span", "three_prime_span",
  "tumor_AD", "ref", "alt",
  "normal_AD", "ref", "alt"
)
for (s in c("start", "end")) {
  for (i in seq_len(nrow(PAIRS))) {
    col <- str_c(s, "_", PAIRS$field[i])
    if (!has_col(col)) next
    x <- events[[col]]
    stop_unless(is.na(x) | n_values(x, ",") == 2L, "not exactly two values", col, x)
    events <- events |>
      mutate(
        "{col}_{PAIRS$part1[i]}" := nth_value(x, ",", 1),
        "{col}_{PAIRS$part2[i]}" := nth_value(x, ",", 2),
        .after = all_of(col)
      ) |>
      select(-all_of(col))
  }
}

# 4. Fields defined as lists that hold one value in every row: kept as one
# value; a row with more stops the run.
SINGLE <- tribble(
  ~field, ~sep,
  "filter", ";",
  "info_ASMID", ",", "info_ASMSEG", ",", "info_BEAOR", ",", "info_BEAPOS", ",", "info_SEGID", ",",
  "region_encode_blacklist_reason", SEP
)
for (s in c("start", "end")) {
  for (i in seq_len(nrow(SINGLE))) {
    col <- str_c(s, "_", SINGLE$field[i])
    if (!has_col(col)) next
    x <- events[[col]]
    stop_unless(is.na(x) | !str_detect(x, fixed(SINGLE$sep[i])), "more than one value", col, x)
  }
}

# 5. Fields with the same value on both breakends of every SV: one column
# without the start_/end_ prefix. A two-breakend event whose records differ
# stops the run.
EVENT_FIELDS <- c(
  "qual", "filter", "info_SVTYPE", "info_PURPLE_JCN", "info_SVID", "info_ASMID", "info_ASMINFO",
  "info_ASMLEN", "info_AVGLEN", "info_SF", "info_DF", "info_VF", "info_MLR", "info_INSALN",
  "info_INSRMP", "info_INSRMRC", "info_INSRMRT", "info_PON_COUNT", "info_TINC_RECOVERED",
  "tumor_GT", "tumor_SF", "tumor_DF", "tumor_VF", "normal_GT", "normal_SF", "normal_DF", "normal_VF",
  "esvee_bnd_PhaseGroupId", "esvee_bnd_PhaseSetId", "esvee_bnd_Length", "esvee_bnd_RefSplitFragments",
  "esvee_bnd_RefDiscFragments", "esvee_bnd_SegmentCount", "esvee_bnd_IncompleteFragments",
  "esvee_bnd_AvgFragmentLength",
  "esvee_bnd_NonPrimaryFragments"
)
event_cols <- character()
for (f in EVENT_FIELDS) {
  s_col <- str_c("start_", f)
  e_col <- str_c("end_", f)
  if (!has_col(s_col) && !has_col(e_col)) next
  sv <- if (has_col(s_col)) blank_na(events[[s_col]]) else rep(NA_character_, nrow(events))
  ev <- if (has_col(e_col)) blank_na(events[[e_col]]) else rep(NA_character_, nrow(events))
  stop_unless(!pair_rows | (is.na(sv) & is.na(ev)) | sv == ev, "differs between the two breakends", f, sv)
  stop_unless(pair_rows | is.na(ev), "end value on a one-breakend event", f, ev)
  events[[f]] <- sv
  events[[s_col]] <- NULL
  events[[e_col]] <- NULL
  event_cols <- c(event_cols, f)
}
log_msg("  {length(event_cols)} event-level columns")

# ---- column order and dictionary ---------------------------------------

key_cols <- c(
  "group_id", "sample_id", "normal_id", "run_mode", "event_key",
  "n_breakends", "in_linx", "linx_start_is_vcf_start", "same_chrom", "length",
  "start_coords", "end_coords", "start_id", "end_id"
)
side_cols <- function(s) str_subset(names(events), str_glue("^{s}_")) |> setdiff(key_cols)
linx_cols <- str_subset(names(events), "^linx_") |> setdiff(key_cols)
pair_cols <- str_subset(names(events), "^xs_pair_")
n_cols <- ncol(events)
events <- events |>
  select(all_of(key_cols), all_of(event_cols), all_of(side_cols("start")), all_of(side_cols("end")),
         all_of(linx_cols), all_of(pair_cols), everything())
stopifnot(ncol(events) == n_cols)

# The two columns a fixed pair was split into: what each part is.
PART_DESC <- c(
  lower = "first value, the lower offset", upper = "second value, the upper offset",
  ref = "first value, reads of the reference allele", alt = "second value, reads of the alt allele",
  five_prime_span = str_c(
    "first value: spread (max - min) of the 5' positions of the unpaired reads ",
    "(ESVEE Breakend.readOrientationRange; -1: no reads)"
  ),
  three_prime_span = str_c(
    "second value: spread (max - min) of the 3' positions of the unpaired reads ",
    "(ESVEE Breakend.readOrientationRange; -1: no reads)"
  )
)
PART_RE <- str_c("^(.+)_(", str_c(names(PART_DESC), collapse = "|"), ")$")
PART_FIELDS <- str_c(PAIRS$field, collapse = "|")

describe <- function(col) {
  base <- str_remove(col, "^(start|end)_")
  side <- str_match(col, "^(start|end)_")[, 2]
  where <- if_else(is.na(side), "", str_c(side, " breakend: "))
  # A part of a split pair: describe the field, then the part.
  part <- str_match(base, PART_RE)
  is_part <- !is.na(part[, 2]) & str_detect(coalesce(part[, 2], ""), str_c("^(", PART_FIELDS, ")$"))
  part_note <- if_else(is_part, str_c(" [", PART_DESC[coalesce(part[, 3], "lower")], "]"), "")
  base <- if_else(is_part, part[, 2], base)
  purple_note <- if_else(
    str_detect(base, "^info_PURPLE_(AF|CN|CN_CHANGE)$"),
    str_c(
      " [this breakend's value: PURPLE writes [start, end] on the start record and [end, start] on the ",
      "end record; the first value of the record is kept (checked against the mate's record)]"
    ),
    ""
  )
  fixed_desc <- c(
    group_id = "oncoanalyser group", sample_id = "tumor sample",
    normal_id = "matched normal sample (empty: tumor-only)",
    run_mode = "PURPLE run mode: TUMOR_GERMLINE (matched normal) or TUMOR (tumor-only)",
    event_key = "sample_id:start_id, unique over the cohort",
    n_breakends = "derived: 2 for an SV with two breakends (MATEID), 1 for a single or inferred breakend",
    in_linx = "derived: TRUE if LINX analyzed the event (LINX keeps FILTER PASS and INFERRED)",
    linx_start_is_vcf_start = str_c(
      "derived: TRUE if LINX's start breakend is the VCF start; FALSE for some SVs whose ",
      "breakends share a position, where LINX's start and end columns refer to the other end"
    ),
    same_chrom = "derived: both breakends on one chromosome",
    length = "derived: end position - start position, both breakends on one chromosome",
    start_coords = "derived: start chromosome:position:orientation",
    end_coords = "derived: end chromosome:position:orientation",
    start_id = "VCF ID of the start breakend (the first of the two in the VCF)",
    end_id = "VCF ID of the end breakend (MATEID of the start)",
    id = "VCF ID", chrom = "VCF CHROM", pos = "VCF POS", ref = "VCF REF", alt = "VCF ALT",
    qual = "VCF QUAL", filter = "VCF FILTER (PURPLE SV VCF)",
    orient = "derived from ALT: 1 if ALT starts with the reference base, else -1 (as ESVEE and LINX)",
    vcf_order = "record number in the PURPLE SV VCF"
  )
  out <- case_when(
    col %in% names(fixed_desc) ~ fixed_desc[col],
    str_starts(col, "xs_pair_") & str_detect(col, "_groups_") ~ str_c(
      "number of other groups with the same two-breakend event (both breakends within ", XS_WINDOW,
      " bases, same orientations, linked by MATEID) in their ", str_match(col, "^xs_pair_(\\w+?)_group")[, 2], " set"
    ),
    str_starts(col, "xs_pair_") ~ "those groups, joined with |",
    base %in% names(fixed_desc) ~ str_c(where, fixed_desc[base]),
    str_starts(base, "info_") ~ str_c(
      where, "PURPLE SV VCF INFO ", str_remove(base, "^info_"), ": ",
      info_desc$description[match(str_remove(base, "^info_"), info_desc$key)]
    ),
    str_starts(base, "tumor_") ~ str_c(
      where, "PURPLE SV VCF FORMAT ", str_remove(base, "^tumor_"), ", tumor: ",
      format_desc$description[match(str_remove(base, "^tumor_"), format_desc$key)]
    ),
    str_starts(base, "normal_") ~ str_c(
      where, "PURPLE SV VCF FORMAT ", str_remove(base, "^normal_"), ", normal: ",
      format_desc$description[match(str_remove(base, "^normal_"), format_desc$key)]
    ),
    str_starts(base, "esvee_bnd_") ~ str_c(
      where, "ESVEE esvee/<sample_id>.esvee.breakend.tsv column ", str_remove(base, "^esvee_bnd_"),
      " (empty for breakends inferred by PURPLE)"
    ),
    base == "linx_bnd_n" ~ str_c(where, "number of LINX breakend.tsv rows (transcripts)"),
    str_starts(base, "linx_bnd_") ~ str_c(
      where, "LINX breakend.tsv column ", str_remove(base, "^linx_bnd_"),
      ", one value per transcript, joined with |"
    ),
    base == "region_encode_blacklist_reason" ~ str_c(
      where, "ENCODE blacklist v2 label (column 4) at the position (one label; checked)"
    ),
    str_starts(base, "region_rmsk_") ~ str_c(
      where, "UCSC RepeatMasker ", str_remove(base, "^region_rmsk_"), " of the repeat(s) at the position, joined with |"
    ),
    base == "region_centromere_distance" ~ str_c(
      where, "bases strictly between the position and the nearest UCSC centromere model interval on ",
      "the chromosome (0: inside an interval or at the base next to one; empty: no model on the chromosome)"
    ),
    base == "region_hmf_sv_prep_blacklist" ~ str_c(
      where, "in Hartwig dna/sv/sv_prep_blacklist.38.bed (0-based BED, read as Hartwig's loader does; ",
      "a Hartwig resource file that no task of this run reads)"
    ),
    base == "region_hmf_unmap" ~ str_c(where, "in Hartwig common/unmap_regions.38.tsv"),
    str_starts(base, "region_") ~ str_c(
      where, "position in track ", str_remove(base, "^region_"), ": ",
      tracks$url[match(str_remove(base, "^region_"), tracks$name)]
    ),
    str_starts(base, "xs_") & str_detect(base, "_groups_exact$") ~ str_c(
      where, "number of other groups with a breakend at this position and orientation in their ",
      str_match(base, "^xs_(\\w+?)_groups")[, 2], " set"
    ),
    str_starts(base, "xs_") & str_detect(base, "_groups_w\\d+$") ~ str_c(
      where, "number of other groups with a breakend within ", XS_WINDOW,
      " bases, same orientation, in their ", str_match(base, "^xs_(\\w+?)_groups")[, 2], " set"
    ),
    str_starts(base, "xs_") & str_detect(base, "_group_list_") ~ str_c(
      where, "those groups, joined with |"
    ),
    str_starts(base, "own_germline_exact") ~ str_c(
      where, "breakends at this position and orientation in the group's own PURPLE germline SV file"
    ),
    str_starts(base, "own_germline_w") ~ str_c(
      where, "breakends within ", XS_WINDOW, " bases, same orientation, in the group's own germline file"
    ),
    str_starts(col, "linx_cluster_driver_") ~ str_c(
      "LINX drivers.tsv column ", str_remove(col, "^linx_cluster_driver_"), " for the event's cluster, joined with |"
    ),
    str_starts(col, "linx_cluster_") ~ str_c("LINX clusters.tsv column ", str_remove(col, "^linx_cluster_")),
    str_starts(col, "linx_fusion_") ~ str_c(
      "LINX fusion.tsv column ", str_remove(col, "^linx_fusion_"),
      " for the fusions this SV takes part in (side: 5p or 3p), joined with |"
    ),
    str_starts(col, "linx_link_") ~ str_c(
      "LINX links.tsv column ", str_remove(col, "^linx_link_"), " for the links of this SV, joined with |"
    ),
    str_starts(col, "linx_") ~ str_c("LINX svs.tsv column ", str_remove(col, "^linx_"), " (as LINX wrote it)"),
    .default = NA_character_
  )
  unname(str_c(out, part_note, purple_note))
}

# Storage type of each column, for the Parquet copy: string, integer, double,
# boolean, or list<...> for lists of any length (LIST_SEP in the TSV;
# LIST_COUNT names the column holding the number of elements, where one
# exists). VCF fields follow the header's Number and Type, as unpacked
# above; other columns follow the source, then all of their values.
RE_INT <- "^-?[0-9]+$"
RE_DBL <- "^-?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][-+]?[0-9]+)?$|^-?(NaN|Inf|Infinity|nan|inf)$"
RE_BOOL <- "^(true|false|TRUE|FALSE)$"
value_class <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  case_when(
    length(x) == 0 ~ "string",
    all(str_detect(x, RE_BOOL)) ~ "boolean",
    all(str_detect(x, RE_INT)) ~ "integer",
    all(str_detect(x, RE_DBL) | str_detect(x, RE_INT)) ~ "double",
    .default = "string"
  )
}
VCF_TYPE <- c(Integer = "integer", Float = "double", String = "string", Character = "string", Flag = "boolean")
# VCF header types that the values contradict (any other mismatch stops
# tsvToParquet.py): ESVEE declares SEGID Integer but writes
# <assembly id>_<segment index>, e.g. 32022_0.
VCF_TYPE_OVERRIDE <- c(SEGID = "string")
FIXED_TYPE <- c(
  group_id = "string", sample_id = "string", normal_id = "string", run_mode = "string",
  event_key = "string", n_breakends = "integer", in_linx = "boolean", linx_start_is_vcf_start = "boolean",
  same_chrom = "boolean", length = "integer", coords = "string", id = "string", chrom = "string",
  pos = "integer", orient = "integer", ref = "string", alt = "string", qual = "double", filter = "string",
  vcf_order = "integer", linx_vcfIdStart = "string", linx_vcfIdEnd = "string"
)
# LINX svs.tsv text fields holding ;-joined items (1 to N).
LINX_SEMI_LISTS <- c("linx_clusterReason", "linx_lineTypeStart", "linx_lineTypeEnd", "linx_geneStart", "linx_geneEnd")
LIST_FAMILIES <- c(
  "^(start|end)_linx_bnd_" = "{side}_linx_bnd_n", "^linx_fusion_" = "linx_fusion_n",
  "^linx_link_" = "linx_link_n", "^linx_cluster_driver_" = "linx_cluster_driver_n"
)
column_type <- function(col) {
  base <- str_remove(col, "^(start|end)_")
  side <- str_match(col, "^(start|end)_")[, 2]
  vcf <- str_match(base, "^(info|tumor|normal)_(.+)$")
  family <- names(LIST_FAMILIES)[str_detect(col, names(LIST_FAMILIES))]
  res <- function(type, sep = NA_character_, count = NA_character_) {
    tibble(type = as.character(type), list_sep = sep, list_count_column = as.character(count))
  }
  if (col %in% names(FIXED_TYPE)) return(res(FIXED_TYPE[[col]]))
  if (col %in% LINX_SEMI_LISTS) return(res("list<string>", ";"))
  if (base %in% names(FIXED_TYPE)) return(res(FIXED_TYPE[[base]]))
  if (length(family) == 1) {
    count <- str_glue(LIST_FAMILIES[[family]], side = coalesce(side, ""))
    if (col == count) return(res("integer"))
    elems <- unlist(str_split(events[[col]][!is.na(events[[col]])], fixed(SEP)))
    elem <- if (str_detect(col, "(vcfId|VcfId)$")) "string" else value_class(elems)
    return(res(str_glue("list<{elem}>"), SEP, count))
  }
  if (!is.na(vcf[, 1])) {
    key <- vcf[, 3]
    part <- str_match(key, PART_RE)
    if (!is.na(part[, 2]) && str_c(vcf[, 2], "_", part[, 2]) %in% PAIRS$field) key <- part[, 2]
    defs <- if (vcf[, 2] == "info") info_desc else format_desc
    d <- defs |> filter(key == !!key)
    if (nrow(d) != 1) stop(str_glue("No VCF header definition for {col}"), call. = FALSE)
    elem <- if (key %in% names(VCF_TYPE_OVERRIDE)) VCF_TYPE_OVERRIDE[[key]] else VCF_TYPE[[d$type]]
    one_value <- d$number %in% c("0", "1") || key != vcf[, 3] ||
      key %in% c("PURPLE_AF", "PURPLE_CN", "PURPLE_CN_CHANGE") ||
      str_c(vcf[, 2], "_", key) %in% SINGLE$field
    if (d$number == "0") return(res("boolean"))
    if (one_value) return(res(elem))
    return(res(str_glue("list<{elem}>"), ","))
  }
  if (str_detect(base, "^region_rmsk_") || str_detect(col, "_group_list_w[0-9]+$")) {
    return(res("list<string>", SEP))
  }
  if (base == "region_encode_blacklist_reason") return(res("string"))
  if (base == "region_centromere_distance") return(res("integer"))
  if (str_starts(base, "region_")) return(res("boolean"))
  if (str_starts(base, "xs_") || str_starts(col, "xs_pair_") || str_starts(base, "own_germline_")) {
    return(res("integer"))
  }
  res(value_class(events[[col]]))
}

set_desc <- c(
  somatic = "PURPLE purple.sv.vcf.gz records without FILTER INFERRED",
  germline = "PURPLE purple.sv.germline.vcf.gz records (groups with a matched normal)",
  unfiltered = "ESVEE esvee.unfiltered.vcf.gz records, any FILTER"
)
dict <- tibble(column = names(events)) |>
  mutate(
    description = describe(column),
    description = if_else(
      column %in% event_cols,
      str_c(description, " [one value per event: the same on both breakends of an SV (checked)]"),
      description
    ),
    description = if_else(
      str_detect(column, "xs_"),
      str_c(description, ". Sets: ", str_c(names(set_desc), ": ", set_desc, collapse = "; ")),
      description
    ),
    description = str_replace_all(description, fixed("joined with |"), "a list (in the TSV joined with |)"),
    non_empty = map_int(column, \(c) sum(!is.na(events[[c]]) & events[[c]] != ""))
  ) |>
  bind_cols(map(names(events), column_type) |> list_rbind())
undescribed <- dict$column[is.na(dict$description)]
if (length(undescribed) > 0) {
  stop("Columns without a description: ", str_c(undescribed, collapse = ", "), call. = FALSE)
}
log_msg("  types: {str_c(names(table(dict$type)), ' ', table(dict$type), collapse = ', ')}")

# Written without quoting, so every reader sees the values as they are; a
# value with a tab, newline or double quote would need quoting: stop.
quote_cols <- names(events)[map_lgl(events, \(x) any(str_detect(as.character(x), "[\t\n\"]"), na.rm = TRUE))]
if (length(quote_cols) > 0) {
  stop("Tab, newline or double quote in values of: ", str_c(quote_cols, collapse = ", "), call. = FALSE)
}
write_tsv(events, out_file, na = "", quote = "none")
write_tsv(dict, dict_file, na = "", quote = "none")
log_msg("{nrow(events)} events x {ncol(events)} columns -> {out_file}")
log_msg("column dictionary -> {dict_file}")

# ---- manifest ----------------------------------------------------------

manifest_file <- path(db_dir, "manifest.tsv")
table_desc <- str_c(
  "One row per somatic SV event of PURPLE's SV VCF (all FILTER values; nothing left out), ",
  "with every field PURPLE, ESVEE and LINX wrote for it, region tracks and the same breakends in ",
  "other groups; fixed pairs split, repeated columns dropped. Derived columns are marked in the ",
  "column dictionary."
)
own_rows <- tibble(
  file = c("sv_cohort_filtered.parquet", "sv_cohort_filtered.tsv.gz", "sv_cohort_filtered.columns.tsv"),
  tool = c("PURPLE, ESVEE, LINX, region tracks", "PURPLE, ESVEE, LINX, region tracks", ""),
  mode = c("somatic", "somatic", ""),
  groups = c(rep(as.character(n_distinct(events$group_id)), 2), ""),
  rows = c(rep(as.character(nrow(events)), 2), as.character(nrow(dict))),
  source = c(
    "sv_cohort_filtered.tsv.gz, typed by scripts/tsvToParquet.py",
    "out/<group_id>/purple/<sample_id>.purple.sv.vcf.gz and the files named in the column dictionary", ""
  ),
  description = c(
    str_c(table_desc, " Typed (column dictionary); lists as arrays."),
    str_c(table_desc, " As text; lists joined with | (VCF lists with ,)."),
    "Column dictionary of sv_cohort_filtered: description, type, list separator and count column, non-empty values."
  )
)
if (file_exists(manifest_file)) {
  manifest <- read_chr(manifest_file) |>
    filter(!file %in% own_rows$file) |>
    bind_rows(own_rows)
  write_tsv(manifest, manifest_file, na = "")
}
log_msg("Done")
