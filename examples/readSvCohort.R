# Bootstrap for analysis of the Prism SV cohort table (postOA.sh):
# reads res/prism/db/sv_cohort_filtered.parquet and its column dictionary,
# and defines helpers for the list columns. Source it at the top of an
# analysis script, from the project folder:
#
#   source("Prism/examples/readSvCohort.R")
#   # or, for another location:
#   SV_DB_DIR <- "path/to/res/prism/db"; source("Prism/examples/readSvCohort.R")
#
# The workbook for people (PASS events, LINX genes and fusions, evidence
# flags, high-confidence set) is made from this table by
# Prism/scripts/svReport.R; see Prism/docs/SV_REPORT.md.
#
# Defines:
#   sv        the table: one row per somatic SV event, typed (integer,
#             double, logical, character; list columns hold vectors)
#   sv_dict   the column dictionary: column, description, type, list_sep,
#             list_count_column, non_empty
#   sv_find(pattern)        dictionary rows whose column or description
#                           matches PATTERN (regex, case-insensitive)
#   sv_transcripts(events)  LINX breakend rows: one row per event, breakend
#                           side and transcript, with gene_role (driver
#                           gene panel: ONCO, TSG or NA)
#   sv_fusions(events)      LINX fusion rows: one row per event and fusion
#                           side, with gene_5p and gene_3p (the genes at
#                           the fusion's 5' and 3' breakends). A fusion is
#                           listed for the SV of its 5' breakend (side
#                           "5p") and of its 3' breakend ("3p"), so a
#                           fusion made by one SV has two rows on that
#                           event; filter(side == "5p") gives each fusion
#                           once
#   sv_links(events)        LINX link rows: one row per event and link
#   sv_cluster_drivers(events)
#                           LINX driver rows of the event's cluster
#   sv_flatten(events)      the same rows with every list column joined
#                           into one string with "|" (no value contains
#                           it; checked). A missing element becomes an
#                           empty slot (positions stay aligned); an empty
#                           list, NA
#
# Layout (details in Prism/docs/prism_methods.md and the dictionary):
#   key columns first: the event (group, sample, FILTER, type, positions),
#   the genes LINX attaches to each breakend (containing it or starting
#   up to 100 kb after it; not the nearest gene) with their role and
#   LINX annotation, fusions, cluster drivers and category, copy number
#   and support
#   no prefix         per event (identity, LINX, and VCF fields equal on
#                     both breakends: qual, filter, info_*, tumor_*, ...)
#   start_* / end_*   per breakend; end_* is empty for single and inferred
#                     breakends, except LINX columns LINX fills from the
#                     alignments of a single breakend's inserted sequence
#                     (end_linx_bnd_* and end_gene_role on a few single
#                     breakends in a fusion)
#   linx_*            LINX svs.tsv row, cluster, and list columns
#   {side}_region_*   region tracks at the breakend (GIAB, ENCODE,
#                     RepeatMasker, Hartwig, centromere distance)
#   {side}_xs_*, {side}_own_germline_*, xs_pair_*
#                     the same breakend, or the same pair, in other files
#   {side}_bam_*      reads at the breakend recounted in the BAMs: own
#                     tumor (bam_tumor_*), own normal (bam_normal_*),
#                     other groups (bam_other_normals_*, bam_other_tumors_*)
#   {side}_<source>_w100_*, <source>_bp100_*, <source>_ro50_*
#                     population SV call sets (gnomad_sv, kgp_sv,
#                     dbvar_nstd166, dbvar_nstd200): sites near the
#                     breakend, sites matching the event
# Events: filter == "INFERRED" are PURPLE copy-number breakpoints, not SV
# calls; filter == "PON" are SV calls LINX does not analyze (in_linx FALSE).

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
})

if (!exists("SV_DB_DIR")) SV_DB_DIR <- "res/prism/db"

sv <- nanoparquet::read_parquet(path(SV_DB_DIR, "sv_cohort_filtered.parquet")) |> as_tibble()
sv_dict <- read_tsv(
  path(SV_DB_DIR, "sv_cohort_filtered.columns.tsv"),
  col_types = cols(.default = col_character()), quote = "", show_col_types = FALSE
)
stopifnot(identical(names(sv), sv_dict$column))

sv_find <- function(pattern) {
  sv_dict |>
    filter(str_detect(column, regex(pattern, ignore_case = TRUE)) |
             str_detect(description, regex(pattern, ignore_case = TRUE))) |>
    select(column, type, description)
}

# One row per element of the list columns whose names start with PREFIX
# (they line up by position), with the event's key columns in front.
# PREFIX is stripped from the column names. EXTRA: other list columns
# aligned with the family, as c(new_name = "column").
unnest_family <- function(events, prefix, keep = c("group_id", "sample_id", "event_key"), extra = character()) {
  count_col <- str_c(prefix, "n")
  list_cols <- names(events)[str_starts(names(events), prefix) & map_lgl(events, is.list)]
  events |>
    filter(!is.na(.data[[count_col]]), .data[[count_col]] > 0) |>
    select(all_of(keep), all_of(list_cols), all_of(extra)) |>
    unnest(c(all_of(list_cols), all_of(names(extra)))) |>
    rename_with(\(x) str_remove(x, str_c("^", prefix)), all_of(list_cols))
}

sv_transcripts <- function(events = sv) {
  bind_rows(
    unnest_family(events, "start_linx_bnd_", extra = c(gene_role = "start_gene_role")) |>
      mutate(side = "start", .after = event_key),
    unnest_family(events, "end_linx_bnd_", extra = c(gene_role = "end_gene_role")) |>
      mutate(side = "end", .after = event_key)
  )
}
sv_fusions <- function(events = sv) {
  unnest_family(events, "linx_fusion_", extra = c(gene_5p = "fusion_gene_5p", gene_3p = "fusion_gene_3p"))
}
sv_links <- function(events = sv) unnest_family(events, "linx_link_")
sv_cluster_drivers <- function(events = sv) unnest_family(events, "linx_cluster_driver_")

# One row per event as before; each list column becomes one string.
sv_flatten <- function(events = sv) {
  list_cols <- names(events)[map_lgl(events, is.list)]
  clash <- list_cols[map_lgl(list_cols, \(c) any(str_detect(unlist(events[[c]]), fixed("|")), na.rm = TRUE))]
  if (length(clash) > 0) stop("| occurs inside values of: ", str_c(clash, collapse = ", "), call. = FALSE)
  events |>
    mutate(across(all_of(list_cols), \(col) {
      map_chr(col, \(x) if (length(x) == 0) NA_character_ else str_c(coalesce(as.character(x), ""), collapse = "|"))
    }))
}

message(
  "sv: ", nrow(sv), " events x ", ncol(sv), " columns from ", path(SV_DB_DIR, "sv_cohort_filtered.parquet"),
  "; ", sum(sv$filter != "INFERRED", na.rm = TRUE), " SV calls (filter != INFERRED)"
)

# ---- examples ------------------------------------------------------------
#
# SV calls only (no inferred copy-number breakpoints), PASS only:
#   calls <- sv |> filter(filter != "INFERRED")
#   pass <- sv |> filter(filter == "PASS")
#
# Events per group and SV type:
#   pass |> count(group_id, linx_type)
#
# Columns about split reads:
#   sv_find("split")
#
# Genes with a disruptive breakend in a canonical transcript, per group:
#   sv_transcripts(pass) |>
#     filter(canonical, disruptive) |>
#     distinct(group_id, gene) |>
#     count(gene, sort = TRUE)
#
# Fusions LINX reports:
#   sv_fusions(pass) |> filter(reported) |> distinct(group_id, name)
#
# Events breaking one gene, or any tumor suppressor (TSG), at either end:
#   sv_transcripts(pass) |> filter(gene == "CD274") |> distinct(event_key, side, regionType, disruptive)
#   sv_transcripts(pass) |> filter(gene_role == "TSG", canonical, disruptive) |> count(gene, sort = TRUE)
#
# Fusions with a given 3' gene:
#   sv_fusions(pass) |> filter(gene_3p == "JAK1") |> distinct(group_id, name, phased, reportedType)
#
# Lists as strings, one row per event (for a spreadsheet or a join):
#   sv_flatten(pass) |> select(event_key, start_linx_bnd_gene, start_gene_role, start_linx_bnd_biotype)
#
# Breakends in GIAB low-mappability regions:
#   sv |> count(start_region_giab_lowmappability, filter)
#
# Two-breakend PASS events with split reads in the group's own normal
# (bam_normal_* is empty for tumor-only groups):
#   pass |> filter(start_bam_normal_split_n > 0 | end_bam_normal_split_n > 0) |> count(group_id)
#
# Events matching a gnomAD-SV site by both breakpoints, with the site IDs:
#   pass |> filter(gnomad_sv_bp100_n > 0) |> select(event_key, info_SVTYPE, gnomad_sv_bp100_ids, gnomad_sv_bp100_max_af)
