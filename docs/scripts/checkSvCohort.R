#!/usr/bin/env Rscript
# Check sv_cohort_filtered (scripts/svCohort.R, scripts/tsvToParquet.py)
# against the files it was built from. Evidence for docs/prism_methods.md;
# rerun after a change to svCohort.R or a pipeline upgrade.
#
#   Rscript --vanilla Prism/docs/scripts/checkSvCohort.R [OUT_DIR] [DB_DIR]
#
# OUT_DIR defaults to out. DB_DIR defaults to res/prism/db in the folder
# that holds OUT_DIR (where postOA.sh writes). It reads the PURPLE SV
# VCFs, the ESVEE breakend.tsv files, the table (TSV and Parquet) and its
# column dictionary.
#
# The table is unpacked (svCohort.R, "unpack"): fixed pairs split, PURPLE's
# per-breakend lists reduced to each breakend's own value, fields equal on
# both breakends stored once per event, columns that repeat another column
# dropped, and the VCF missing value "." stored as empty. This script undoes
# each step to rebuild every VCF record and every ESVEE row, then compares
# them with the source files. A VCF "." compares equal to empty.
#
# Per group, a check passes where its column is TRUE:
#   records_once    every record of purple.sv.vcf.gz is the start or the
#                   end of exactly one event, and no other ID is in the
#                   table
#   fields_equal    CHROM, POS, REF, ALT, QUAL, FILTER of each breakend
#                   equal the VCF record
#   info_equal      the rebuilt INFO key=value pairs equal the VCF INFO of
#                   each record (flags as "true")
#   format_equal    the rebuilt FORMAT values of the tumor and the normal
#                   equal the VCF
#   esvee_equal     every column of every ESVEE breakend.tsv row of a
#                   record equals the rebuilt value (StrandBias: FORMAT SB
#                   of each sample to two decimals); INFERRED records have
#                   no ESVEE row
#   linx_ok         in_linx events = LINX svs.tsv rows, and every event
#                   with linx_start_is_vcf_start FALSE has both breakends
#                   at one position
#   linx_bnd_coords_ok
#                   every LINX breakend row (placed by isStart) has the
#                   coords of the end it is placed on; end rows of
#                   one-breakend events are counted, not compared
#                   (sgl_events_end_rows, sgl_end_rows): LINX writes them
#                   at a location that is not a VCF record
# Then: parquet_equal, every cell of the Parquet copy equals the TSV cell
# converted to the type in the column dictionary (an independent
# reimplementation in R of tsvToParquet.py); the column dictionary, the
# events by type and FILTER, and counts of the region and cross-sample
# columns.

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
})
options(width = 200)

# SV VCF lines (assembled insert sequences) are longer than readr's
# default buffer for compressed files.
Sys.setenv(VROOM_CONNECTION_SIZE = 2^26)

args <- commandArgs(trailingOnly = TRUE)
out_dir <- path_abs(if (length(args) > 0) args[1] else "out")
db_dir <- if (length(args) > 1) path_abs(args[2]) else path(path_dir(out_dir), "res", "prism", "db")

read_chr <- function(file, ...) {
  read_tsv(
    file, col_types = cols(.default = col_character()), quote = "",
    na = character(), show_col_types = FALSE, lazy = FALSE, trim_ws = FALSE, ...
  )
}

read_vcf_body <- function(file) {
  n_skip <- sum(str_starts(read_lines(file, n_max = 20000), "##"))
  read_chr(file, skip = n_skip) |> rename(CHROM = `#CHROM`)
}

ev <- read_chr(path(db_dir, "sv_cohort_filtered.tsv.gz"))
dict <- read_chr(path(db_dir, "sv_cohort_filtered.columns.tsv"))
samples <- read_chr(path(db_dir, "samples.tsv"))
lx_svs <- read_chr(path(db_dir, "data", "linx_somatic_svs.tsv.gz"))

# ---- undo the unpacking ---------------------------------------------------

OTHER <- c(start = "end", end = "start")
PAIRS <- tribble(
  ~field, ~part1, ~part2,
  "info_CIPOS", "lower", "upper",
  "info_IHOMPOS", "lower", "upper",
  "info_TPR", "five_prime_span", "three_prime_span",
  "tumor_AD", "ref", "alt",
  "normal_AD", "ref", "alt"
)
PURPLE_KEYS <- c("info_PURPLE_AF", "info_PURPLE_CN", "info_PURPLE_CN_CHANGE")
EVENT_RE <- "^(qual|filter|info_.+|tumor_.+|normal_.+|esvee_bnd_.+)$"
# ESVEE breakend.tsv columns that svCohort.R drops as repeats: the rebuilt
# record field holding the same value (mate_*: the other breakend), and the
# separator ESVEE uses where the VCF uses ",".
ESVEE_FROM <- tribble(
  ~esvee, ~from, ~sep,
  "AssemblyId", "info_ASMID", NA,
  "AssemblyInfo", "info_ASMINFO", NA,
  "Type", "info_SVTYPE", NA,
  "Chromosome", "chrom", NA,
  "Position", "pos", NA,
  "Orientation", "orient", NA,
  "MateId", "mate_id", NA,
  "MateChr", "mate_chrom", NA,
  "MatePos", "mate_pos", NA,
  "MateOrient", "mate_orient", NA,
  "ConfidenceInterval", "info_CIPOS", NA,
  "InexactOffset", "info_IHOMPOS", NA,
  "Homology", "info_HOMSEQ", NA,
  "Qual", "qual", NA,
  "SplitFragments", "info_SF", NA,
  "DiscFragments", "info_DF", NA,
  "SequenceLength", "info_ASMLEN", NA,
  "SegmentIndex", "info_ASMSEG", NA,
  "SequenceIndex", "info_BEAPOS", NA,
  "AlignedBases", "info_SEGALEN", NA,
  "MapQual", "info_SEGMAPQ", NA,
  "Score", "info_SEGSCO", NA,
  "AdjAlignedBases", "info_SEGRL", NA,
  "FacingBreakendIds", "info_ASMLNKS", ";",
  "AltAlignments", "info_INSALN", NA,
  "UniqueFragPos", "info_UFP", NA,
  "ThreePrimeRange", "info_TPR", ";"
)

col_or_empty <- function(df, col) if (col %in% names(df)) df[[col]] else rep("", nrow(df))

# The VCF records of one side of the events, with every VCF field and ESVEE
# column under its original name.
rebuild_records <- function(e, side) {
  o <- OTHER[[side]]
  rec <- e |> filter(.data[[str_c(side, "_id")]] != "")
  own <- rec |>
    select(starts_with(str_c(side, "_"))) |>
    rename_with(\(x) str_remove(x, str_glue("^{side}_")))
  out <- bind_cols(own, rec |> select(matches(EVENT_RE), -any_of("normal_id"))) |>
    mutate(
      mate_id = col_or_empty(rec, str_c(o, "_id")), mate_chrom = col_or_empty(rec, str_c(o, "_chrom")),
      mate_pos = col_or_empty(rec, str_c(o, "_pos")), mate_orient = col_or_empty(rec, str_c(o, "_orient"))
    )
  for (i in seq_len(nrow(PAIRS))) {
    p1 <- str_c(PAIRS$field[i], "_", PAIRS$part1[i])
    p2 <- str_c(PAIRS$field[i], "_", PAIRS$part2[i])
    if (!p1 %in% names(out)) next
    out[[PAIRS$field[i]]] <- if_else(out[[p1]] == "" & out[[p2]] == "", "", str_c(out[[p1]], ",", out[[p2]]))
    out[[p1]] <- NULL
    out[[p2]] <- NULL
  }
  for (k in PURPLE_KEYS) {
    if (!k %in% names(out)) next
    mate <- col_or_empty(rec, str_c(o, "_", k))
    out[[k]] <- case_when(out[[k]] == "" ~ "", mate == "" ~ out[[k]], .default = str_c(out[[k]], ",", mate))
  }
  out$info_MATEID <- out$mate_id
  out$info_BEOR <- if_else(out$filter == "INFERRED", "", out$orient)
  for (i in seq_len(nrow(ESVEE_FROM))) {
    v <- col_or_empty(out, ESVEE_FROM$from[i])
    if (!is.na(ESVEE_FROM$sep[i])) v <- str_replace_all(v, fixed(","), ESVEE_FROM$sep[i])
    out[[str_c("esvee_bnd_", ESVEE_FROM$esvee[i])]] <- v
  }
  out
}

# ---- compare with the source files ------------------------------------------

dot_empty <- function(x) if_else(x == ".", "", x)

# key/value pairs of the VCF INFO, as svCohort.R splits them; "." values
# left out (stored as empty).
info_pairs <- function(id, info) {
  tibble(id = id, kv = str_split(info, ";")) |>
    unnest(kv) |>
    filter(kv != "" & kv != ".") |>
    transmute(
      id, key = str_remove(kv, "=.*$"),
      value = if_else(str_detect(kv, "="), str_remove(kv, "^[^=]*="), "true")
    ) |>
    filter(value != ".")
}

format_pairs <- function(id, format, values) {
  tibble(id = id, key = str_split(format, ":"), value = str_split(values, ":")) |>
    unnest(c(key, value)) |>
    filter(value != ".")
}

# key/value pairs of the record columns <prefix><key>, empty cells left out.
table_pairs <- function(rec, prefix) {
  rec |>
    select(id, starts_with(prefix)) |>
    pivot_longer(-id, names_to = "key", values_to = "value") |>
    filter(value != "") |>
    mutate(key = str_remove(key, str_glue("^{prefix}")))
}

same_pairs <- function(a, b, what) {
  only_a <- anti_join(a, b, by = c("id", "key", "value"))
  only_b <- anti_join(b, a, by = c("id", "key", "value"))
  if (nrow(only_a) + nrow(only_b) > 0) {
    message(what, ": only in the source ", nrow(only_a), ", only in the table ", nrow(only_b))
    print(head(bind_rows(source = only_a, table = only_b, .id = "in"), 6))
  }
  nrow(only_a) == 0 && nrow(only_b) == 0
}

check_group <- function(group_id, sample_id, normal_id, ...) {
  g <- group_id
  vcf <- read_vcf_body(path(out_dir, g, "purple", str_c(sample_id, ".purple.sv.vcf.gz")))
  esvee <- read_chr(path(out_dir, g, "esvee", str_c(sample_id, ".esvee.breakend.tsv")))
  e <- ev |> filter(group_id == g)
  ids <- c(e$start_id, e$end_id[e$end_id != ""])

  rec <- bind_rows(rebuild_records(e, "start"), rebuild_records(e, "end")) |>
    mutate(across(everything(), \(x) replace_na(x, "")))
  joined <- rec |> inner_join(vcf, by = c(id = "ID"))
  fields_equal <- nrow(joined) == nrow(rec) &&
    all(joined$chrom == joined$CHROM & joined$pos == joined$POS & joined$ref == joined$REF &
          joined$alt == joined$ALT & joined$qual == dot_empty(joined$QUAL) &
          joined$filter == dot_empty(joined$FILTER))

  info_equal <- same_pairs(info_pairs(vcf$ID, vcf$INFO), table_pairs(rec, "info_"), str_c(g, " INFO"))
  format_equal <- same_pairs(
    format_pairs(vcf$ID, vcf$FORMAT, vcf[[sample_id]]), table_pairs(rec, "tumor_"), str_c(g, " tumor FORMAT")
  )
  if (normal_id != "") {
    format_equal <- format_equal && same_pairs(
      format_pairs(vcf$ID, vcf$FORMAT, vcf[[normal_id]]), table_pairs(rec, "normal_"), str_c(g, " normal FORMAT")
    )
  }

  # Every ESVEE column of the records' rows: rebuilt value = file value.
  es <- esvee |> filter(Id %in% rec$id)
  es_rec <- rec |> filter(id %in% es$Id)
  es_rec <- es_rec[match(es$Id, es_rec$id), ]
  inferred_ok <- all((rec$id %in% es$Id) == (rec$filter != "INFERRED"))
  es_cols <- setdiff(names(esvee), c("Id", "StrandBias"))
  missing_cols <- setdiff(str_c("esvee_bnd_", es_cols), names(es_rec))
  col_ok <- map_lgl(es_cols, \(c) {
    tc <- str_c("esvee_bnd_", c)
    tc %in% names(es_rec) && all(es_rec[[tc]] == es[[c]])
  })
  if (length(missing_cols) > 0 || !all(col_ok)) {
    message(g, " ESVEE columns not rebuilt or different: ", str_c(c(missing_cols, es_cols[!col_ok]), collapse = ", "))
  }
  sb_expect <- if_else(
    col_or_empty(es_rec, "normal_SB") == "", es_rec$tumor_SB, str_c(es_rec$tumor_SB, ";", es_rec$normal_SB)
  )
  sb_ok <- all(map2_lgl(es$StrandBias, sb_expect, \(a, b) {
    x <- as.numeric(str_split_1(a, ";"))
    y <- as.numeric(str_split_1(b, ";"))
    length(x) == length(y) && all(abs(x - y) <= 0.005 + 1e-9)
  }))
  esvee_equal <- inferred_ok && length(missing_cols) == 0 && all(col_ok) && sb_ok

  flipped <- e |> filter(linx_start_is_vcf_start == "FALSE")
  linx_ok <- sum(e$in_linx == "TRUE") == sum(lx_svs$group_id == g) &&
    all(flipped$start_pos == flipped$end_pos)

  # LINX breakend rows are placed by isStart; their coords must be those of
  # the end they are placed on (one-breakend events: end rows have none).
  bnd_coords_ok <- function(side) {
    e |>
      filter(.data[[str_c(side, "_linx_bnd_coords")]] != "", .data[[str_c(side, "_id")]] != "") |>
      mutate(c = str_split(.data[[str_c(side, "_linx_bnd_coords")]], fixed("|"))) |>
      select(all_of(str_c(side, "_coords")), c) |>
      unnest(c) |>
      summarize(ok = all(c == .data[[str_c(side, "_coords")]])) |>
      pull(ok)
  }
  sgl_end_rows <- e |> filter(end_id == "", end_linx_bnd_n != "")

  tibble(
    group_id = g,
    vcf_records = nrow(vcf),
    events = nrow(e),
    records_once = length(ids) == nrow(vcf) && !anyDuplicated(ids) && setequal(ids, vcf$ID),
    fields_equal = fields_equal,
    info_equal = info_equal,
    format_equal = format_equal,
    esvee_equal = esvee_equal,
    linx_events = sum(e$in_linx == "TRUE"),
    linx_start_is_vcf_end = nrow(flipped),
    linx_ok = linx_ok,
    linx_bnd_coords_ok = bnd_coords_ok("start") && bnd_coords_ok("end"),
    sgl_events_end_rows = nrow(sgl_end_rows),
    sgl_end_rows = sum(as.integer(sgl_end_rows$end_linx_bnd_n))
  )
}

# ---- LINX tables = the table's LINX columns ----------------------------------
#
# Every LINX value in the table, compared as text with the stacked LINX
# tables of postOA.R (data/linx_somatic_*.tsv.gz), row for row: each source
# row must come back exactly once (as a multiset), with nothing extra.

lx <- c("svs", "clusters", "breakend", "fusion", "links", "drivers") |>
  set_names() |>
  map(\(t) read_chr(path(db_dir, "data", str_c("linx_somatic_", t, ".tsv.gz"))))

# Same rows, as multisets: same columns, same rows with the same counts.
same_rows <- function(a, b, what) {
  b <- b |> select(all_of(names(a)))
  key <- \(d) d |> arrange(across(everything())) |> as.data.frame()
  ok <- nrow(a) == nrow(b) && identical(key(a), key(b))
  if (!ok) {
    message(what, ": table rows ", nrow(a), ", LINX rows ", nrow(b), "; examples that differ:")
    print(head(bind_rows(table = anti_join(a, b, by = names(a)), linx = anti_join(b, a, by = names(a)), .id = "in"), 6))
  }
  ok
}

# The |-lists of the columns starting with PREFIX, one row per element,
# with the columns of KEEP; PREFIX removed from the names.
unnest_lists <- function(e, prefix, keep) {
  count_col <- str_c(prefix, "n")
  cols <- setdiff(str_subset(names(e), str_c("^", prefix)), count_col)
  e <- e |> filter(.data[[count_col]] != "")
  n <- as.integer(e[[count_col]])
  long <- map(set_names(cols), \(c) {
    v <- str_split(e[[c]], fixed("|"))
    v <- map2(v, n, \(x, k) if (length(x) == k) x else stop(str_glue("{c}: list length differs from {count_col}")))
    unlist(v)
  })
  e |>
    select(all_of(keep)) |>
    slice(rep(seq_len(nrow(e)), n)) |>
    bind_cols(as_tibble(long)) |>
    rename_with(\(x) str_remove(x, str_c("^", prefix)), all_of(cols))
}

check_linx <- function() {
  in_lx <- ev |> filter(in_linx == "TRUE")
  res <- list()

  # svs.tsv: the linx_* columns of every LINX event.
  svs_cols <- setdiff(names(lx$svs), c("group_id", "sample_id"))
  res$svs <- same_rows(
    in_lx |> select(group_id, sample_id, all_of(str_c("linx_", svs_cols))) |>
      rename_with(\(x) str_remove(x, "^linx_"), -c(group_id, sample_id)),
    lx$svs, "LINX svs"
  )

  # clusters.tsv: the linx_cluster_* columns, once per cluster.
  cl_cols <- setdiff(names(lx$clusters), c("group_id", "sample_id", "clusterId"))
  res$clusters <- same_rows(
    in_lx |>
      distinct(group_id, sample_id, clusterId = linx_clusterId, across(all_of(str_c("linx_cluster_", cl_cols)))) |>
      rename_with(\(x) str_remove(x, "^linx_cluster_"), starts_with("linx_cluster_")),
    lx$clusters, "LINX clusters"
  )

  # breakend.tsv: start and end lists together, svId from the event.
  bnd <- bind_rows(
    unnest_lists(in_lx, "start_linx_bnd_", c("group_id", "sample_id", svId = "linx_svId")),
    unnest_lists(in_lx, "end_linx_bnd_", c("group_id", "sample_id", svId = "linx_svId"))
  )
  res$breakend <- same_rows(bnd, lx$breakend, "LINX breakend")

  # fusion.tsv: listed on the SV of the 5' breakend (side 5p) and of the 3'
  # breakend (3p); each side must give every fusion once, on the right SV.
  fus <- unnest_lists(in_lx, "linx_fusion_", c("group_id", "sample_id", svId = "linx_svId"))
  bnd_sv <- lx$breakend |> select(group_id, sample_id, id, bnd_sv = svId)
  for (s in c("5p", "3p")) {
    id_col <- if (s == "5p") "fivePrimeBreakendId" else "threePrimeBreakendId"
    f <- fus |> filter(side == s)
    on_sv <- f |>
      left_join(bnd_sv, by = c("group_id", "sample_id", setNames("id", id_col))) |>
      summarize(ok = all(!is.na(bnd_sv) & bnd_sv == svId)) |>
      pull(ok)
    res[[str_c("fusion_", s)]] <- on_sv &&
      same_rows(f |> select(-svId, -side), lx$fusion, str_c("LINX fusion (", s, ")"))
  }

  # links.tsv: listed on lowerSvId, and on upperSvId when it differs.
  lnk <- unnest_lists(in_lx, "linx_link_", c("group_id", "sample_id", svId = "linx_svId"))
  res$links_on_sv <- all(lnk$svId == lnk$lowerSvId | lnk$svId == lnk$upperSvId)
  res$links_lower <- same_rows(lnk |> filter(svId == lowerSvId) |> select(-svId), lx$links, "LINX links (lower)")
  res$links_upper <- same_rows(
    lnk |> filter(svId == upperSvId, upperSvId != lowerSvId) |> select(-svId),
    lx$links |> filter(upperSvId != lowerSvId), "LINX links (upper)"
  )

  # drivers.tsv: the drivers of the event's cluster, once per cluster;
  # drivers outside any cluster (clusterId -1) belong to no event.
  one_per_cluster <- in_lx |> distinct(group_id, sample_id, linx_clusterId, .keep_all = TRUE)
  drv <- unnest_lists(one_per_cluster, "linx_cluster_driver_", c("group_id", "sample_id", clusterId = "linx_clusterId"))
  same_lists <- in_lx |>
    summarize(n = n_distinct(linx_cluster_driver_gene, linx_cluster_driver_eventType), .by = c(group_id, sample_id, linx_clusterId))
  if (any(same_lists$n > 1)) message("LINX drivers: events of one cluster carry different driver lists")
  res$drivers_same_in_cluster <- all(same_lists$n == 1)
  drv_src <- lx$drivers |> semi_join(in_lx, by = c("group_id", "sample_id", clusterId = "linx_clusterId"))
  res$drivers <- same_rows(drv, drv_src, "LINX drivers")
  message("LINX drivers outside any event's cluster (not in the table): ", nrow(lx$drivers) - nrow(drv_src))

  res <- unlist(res)
  print(res)
  all(res)
}

# ---- Parquet copy = TSV, typed by the dictionary ---------------------------

typed <- function(x, type) {
  x <- na_if(x, "")
  switch(type,
    string = x,
    integer = as.integer(x),
    double = as.numeric(x),
    boolean = toupper(x) == "TRUE",
    stop("unknown type ", type)
  )
}

# Doubles: equal to 1 part in a million. The check is for conversion errors
# (a wrong column, row or parse), which are far larger; it also absorbs R's
# as.numeric, which does not always give the nearest double ("5.754e-03"
# lands one unit in the last place from pyarrow's correct value).
same_values <- function(a, b) {
  if (!is.double(a) || !is.double(b)) return(identical(a, b))
  ok <- !is.na(a)
  length(a) == length(b) && identical(is.na(a), is.na(b)) && identical(is.nan(a), is.nan(b)) &&
    all(a[ok] == b[ok] | abs(a[ok] - b[ok]) <= 1e-6 * abs(a[ok]))
}

parquet_equal <- function() {
  pq <- nanoparquet::read_parquet(path(db_dir, "sv_cohort_filtered.parquet"))
  if (!identical(names(pq), dict$column) || nrow(pq) != nrow(ev)) return(FALSE)
  ok <- map_lgl(seq_len(nrow(dict)), \(i) {
    col <- dict$column[i]
    type <- dict$type[i]
    if (!str_starts(type, "list<")) {
      return(same_values(unname(typed(ev[[col]], type)), unname(pq[[col]])))
    }
    elem <- str_match(type, "^list<(.+)>$")[, 2]
    count_col <- dict$list_count_column[i]
    has_list <- if (count_col != "") coalesce(as.integer(na_if(ev[[count_col]], "")), 0L) > 0 else ev[[col]] != ""
    expect <- map2(ev[[col]], has_list, \(cell, has) {
      if (!has) return(NULL)
      typed(str_split_1(cell, fixed(dict$list_sep[i])), elem)
    })
    got <- map(pq[[col]], \(v) if (is.null(v)) NULL else unname(v))
    length(expect) == length(got) && all(map2_lgl(expect, got, \(x, y) {
      (is.null(x) && is.null(y)) || (!is.null(x) && !is.null(y) && same_values(x, y))
    }))
  })
  if (!all(ok)) message("Parquet columns that differ from the TSV: ", str_c(dict$column[!ok], collapse = ", "))
  all(ok)
}

# ---- run ----------------------------------------------------------------------

message("events ", nrow(ev), " x ", ncol(ev), "; dictionary rows ", nrow(dict),
        "; columns without description ", sum(dict$description == ""),
        "; columns empty in every row ", sum(dict$non_empty == "0"))
print(dict |> filter(non_empty == "0") |> pull(column))
print(count(dict, type))

res <- samples |> pmap(check_group) |> list_rbind()
print(res, n = Inf)
checks <- c(
  "records_once", "fields_equal", "info_equal", "format_equal", "esvee_equal", "linx_ok",
  "linx_bnd_coords_ok"
)
linx_ok <- check_linx()
message("linx_tables_equal: ", linx_ok)
pq_ok <- parquet_equal()
message("parquet_equal: ", pq_ok)
all_pass <- all(as.matrix(res[checks])) && linx_ok && pq_ok
message("all checks pass: ", all_pass)

message("events by breakends, LINX and FILTER:")
print(count(ev, n_breakends, in_linx, linx_start_is_vcf_start, filter), n = Inf)

message("breakends (start) in each region track (TRUE):")
print(
  ev |>
    summarize(across(matches("^start_region_(giab|encode_blacklist$|hmf)"), \(x) sum(x == "TRUE"))) |>
    pivot_longer(everything()),
  n = Inf
)
message("events with a hit in another group (start, pair) or in the own germline file:")
print(
  ev |>
    summarize(across(
      matches("^start_xs_.*_groups_w50$|^xs_pair_.*_groups_w50$|^start_own_germline"),
      \(x) sum(as.integer(x) > 0, na.rm = TRUE)
    )) |>
    pivot_longer(everything()),
  n = Inf
)

# Exit status 1 if any check failed (postOA.sh stops on it).
if (!all_pass) quit(status = 1)
