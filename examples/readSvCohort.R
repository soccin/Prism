# Bootstrap for analysis of the Prism SV cohort table (postOA.sh):
# reads res/prism/db/sv_cohort_filtered.parquet and its column dictionary,
# and defines helpers for the list columns. Source it at the top of an
# analysis script, from the project folder:
#
#   source("Prism/examples/readSvCohort.R")
#   # or, for another location:
#   SV_DB_DIR <- "path/to/res/prism/db"; source("Prism/examples/readSvCohort.R")
#
# Defines:
#   sv        the table: one row per somatic SV event, typed (integer,
#             double, logical, character; list columns hold vectors)
#   sv_dict   the column dictionary: column, description, type, list_sep,
#             list_count_column, non_empty
#   sv_find(pattern)        dictionary rows whose column or description
#                           matches PATTERN (regex, case-insensitive)
#   sv_transcripts(events)  LINX breakend rows: one row per event, breakend
#                           side and transcript
#   sv_fusions(events)      LINX fusion rows: one row per event and fusion
#   sv_links(events)        LINX link rows: one row per event and link
#   sv_cluster_drivers(events)
#                           LINX driver rows of the event's cluster
#
# Layout (details in Prism/docs/prism_methods.md and the dictionary):
#   no prefix         per event (identity, LINX, and VCF fields equal on
#                     both breakends: qual, filter, info_*, tumor_*, ...)
#   start_* / end_*   per breakend; end_* is empty for single and inferred
#                     breakends
#   linx_*            LINX svs.tsv row, cluster, and list columns
#   xs_*, own_germline_*  the same breakend or pair in other files
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
# PREFIX is stripped from the column names.
unnest_family <- function(events, prefix, keep = c("group_id", "sample_id", "event_key")) {
  count_col <- str_c(prefix, "n")
  list_cols <- names(events)[str_starts(names(events), prefix) & map_lgl(events, is.list)]
  events |>
    filter(!is.na(.data[[count_col]]), .data[[count_col]] > 0) |>
    select(all_of(keep), all_of(list_cols)) |>
    unnest(all_of(list_cols)) |>
    rename_with(\(x) str_remove(x, str_c("^", prefix)), all_of(list_cols))
}

sv_transcripts <- function(events = sv) {
  bind_rows(
    unnest_family(events, "start_linx_bnd_") |> mutate(side = "start", .after = event_key),
    unnest_family(events, "end_linx_bnd_") |> mutate(side = "end", .after = event_key)
  )
}
sv_fusions <- function(events = sv) unnest_family(events, "linx_fusion_")
sv_links <- function(events = sv) unnest_family(events, "linx_link_")
sv_cluster_drivers <- function(events = sv) unnest_family(events, "linx_cluster_driver_")

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
# Breakends in GIAB low-mappability regions:
#   sv |> count(start_region_giab_lowmappability, filter)
