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
# next to the table):
#   start_* / end_*   per breakend: the VCF fields (info_*: INFO, tumor_* and
#                     normal_*: FORMAT of each sample), every column of
#                     ESVEE's breakend.tsv (esvee_bnd_*), LINX's breakend.tsv
#                     rows for that end (linx_bnd_*, one value per
#                     transcript), region tracks (region_*), and the same
#                     breakend in other groups (xs_*) and in the group's own
#                     germline file (own_germline_*)
#   linx_*            the LINX svs.tsv row, its cluster, links, fusions and
#                     cluster drivers
#   xs_pair_*         the same two-breakend event in other groups
# Values are copied as the tools wrote them. Several values for one cell
# (one per transcript, fusion, link or driver) are joined with "|" in the
# order of the source file. Derived columns are documented in the
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

# INFO or FORMAT descriptions from VCF header lines.
header_fields <- function(header, kind) {
  m <- str_match(header, str_glue('^##{kind}=<ID=([^,]+),.*Description="([^"]*)"'))
  tibble(key = m[, 2], description = m[, 3]) |> filter(!is.na(key))
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
info_desc <- groups |> map("info") |> list_rbind() |> distinct(key, .keep_all = TRUE)
format_desc <- groups |> map("format") |> list_rbind() |> distinct(key, .keep_all = TRUE)
rm(groups)

events <- events |>
  mutate(
    event_key = str_c(sample_id, ":", start_id),
    n_breakends = if_else(is.na(end_id), 1L, 2L),
    same_chrom = if_else(is.na(end_id), NA, start_chrom == end_chrom),
    length = if_else(
      same_chrom %in% TRUE, abs(as.numeric(end_pos) - as.numeric(start_pos)), NA_real_
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

# LINX breakend rows (one per transcript) to the event end with the same
# coordinates.
lx_bnd <- linx$breakend
bnd_cols <- setdiff(names(lx_bnd), c("group_id", "sample_id", "svId", "coords"))
for (side in c("start", "end")) {
  side_bnd <- lx_bnd |>
    inner_join(
      events |>
        filter(in_linx) |>
        select(group_id, sample_id, svId = linx_svId, coords = all_of(str_c(side, "_coords"))),
      by = c("group_id", "sample_id", "svId", "coords")
    ) |>
    select(group_id, sample_id, svId, all_of(bnd_cols)) |>
    collapse_rows(c("group_id", "sample_id", "svId"), str_c(side, "_linx_bnd_"))
  events <- events |>
    left_join(side_bnd, by = c("group_id", "sample_id", linx_svId = "svId"))
}
n_bnd <- sum(as.integer(events$start_linx_bnd_n), as.integer(events$end_linx_bnd_n), na.rm = TRUE)
if (n_bnd != nrow(lx_bnd)) {
  stop(str_glue("{nrow(lx_bnd) - n_bnd} LINX breakend rows not placed on an event end"), call. = FALSE)
}

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
read_hmf_regions <- function(file) {
  b <- read_chr(file)
  as_ranges(b$Chromosome, b$PosStart, b$PosEnd)
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

# Distance to the nearest centromere model on the same chromosome (0 inside).
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
has_germline <- points$group_id %in% unique(sets$germline$group_id)
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

# ---- column order and dictionary ---------------------------------------

key_cols <- c(
  "group_id", "sample_id", "normal_id", "run_mode", "vcf_normal_sample", "event_key",
  "n_breakends", "in_linx", "linx_start_is_vcf_start", "same_chrom", "length",
  "start_coords", "end_coords", "start_id", "end_id"
)
side_cols <- function(s) str_subset(names(events), str_glue("^{s}_")) |> setdiff(key_cols)
linx_cols <- str_subset(names(events), "^linx_") |> setdiff(key_cols)
pair_cols <- str_subset(names(events), "^xs_pair_")
n_cols <- ncol(events)
events <- events |>
  select(all_of(key_cols), all_of(side_cols("start")), all_of(side_cols("end")),
         all_of(linx_cols), all_of(pair_cols), everything())
stopifnot(ncol(events) == n_cols)

describe <- function(col) {
  base <- str_remove(col, "^(start|end)_")
  side <- str_match(col, "^(start|end)_")[, 2]
  where <- if_else(is.na(side), "", str_c(side, " breakend: "))
  fixed_desc <- c(
    group_id = "oncoanalyser group", sample_id = "tumor sample",
    normal_id = "matched normal sample (empty: tumor-only)",
    run_mode = "PURPLE run mode: TUMOR_GERMLINE (matched normal) or TUMOR (tumor-only)",
    vcf_normal_sample = "the second sample column of the PURPLE SV VCF (the normal), if any",
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
    str_starts(base, "linx_bnd_n") ~ str_c(where, "number of LINX breakend.tsv rows (transcripts)"),
    str_starts(base, "linx_bnd_") ~ str_c(
      where, "LINX breakend.tsv column ", str_remove(base, "^linx_bnd_"),
      ", one value per transcript, joined with |"
    ),
    str_starts(base, "region_encode_blacklist_reason") ~ str_c(where, "ENCODE blacklist v2 label(s)"),
    str_starts(base, "region_rmsk_") ~ str_c(
      where, "UCSC RepeatMasker ", str_remove(base, "^region_rmsk_"), " of the repeat(s) at the position, joined with |"
    ),
    base == "region_centromere_distance" ~ str_c(
      where, "bases to the nearest UCSC centromere model on the chromosome (0: inside)"
    ),
    base == "region_hmf_sv_prep_blacklist" ~ str_c(
      where, "in Hartwig dna/sv/sv_prep_blacklist.38.bed (regions ESVEE prep leaves out)"
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
  unname(out)
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
      str_detect(column, "xs_"),
      str_c(description, ". Sets: ", str_c(names(set_desc), ": ", set_desc, collapse = "; ")),
      description
    ),
    non_empty = map_int(column, \(c) sum(!is.na(events[[c]]) & events[[c]] != ""))
  )
undescribed <- dict$column[is.na(dict$description)]
if (length(undescribed) > 0) {
  warning("Columns without a description: ", str_c(undescribed, collapse = ", "), call. = FALSE)
}

write_tsv(events, out_file, na = "")
write_tsv(dict, dict_file, na = "")
log_msg("{nrow(events)} events x {ncol(events)} columns -> {out_file}")
log_msg("column dictionary -> {dict_file}")

# ---- manifest ----------------------------------------------------------

manifest_file <- path(db_dir, "manifest.tsv")
own_rows <- tibble(
  file = c("sv_cohort_filtered.tsv.gz", "sv_cohort_filtered.columns.tsv"),
  tool = c("PURPLE, ESVEE, LINX, region tracks", ""),
  mode = c("somatic", ""),
  groups = c(as.character(n_distinct(events$group_id)), ""),
  rows = c(as.character(nrow(events)), as.character(nrow(dict))),
  source = c("out/<group_id>/purple/<sample_id>.purple.sv.vcf.gz and the files named in the column dictionary", ""),
  description = c(
    str_c(
      "One row per somatic SV event of PURPLE's SV VCF (all FILTER values; nothing left out), ",
      "with every field PURPLE, ESVEE and LINX wrote for it, region tracks and the same breakends in ",
      "other groups. Derived columns are marked in the column dictionary."
    ),
    "Column dictionary of sv_cohort_filtered.tsv.gz: description and number of non-empty values."
  )
)
if (file_exists(manifest_file)) {
  manifest <- read_chr(manifest_file) |>
    filter(!file %in% own_rows$file) |>
    bind_rows(own_rows)
  write_tsv(manifest, manifest_file, na = "")
}
log_msg("Done")
