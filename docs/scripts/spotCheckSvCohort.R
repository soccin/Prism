#!/usr/bin/env Rscript
# Spot-check the computed columns of sv_cohort_filtered (scripts/svCohort.R)
# on a random sample of events: each column is recomputed here by a second,
# independent method (plain filtering per breakend instead of interval
# overlaps; bcftools and samtools instead of the R and pysam readers) and
# compared with the table. Evidence for docs/prism_methods.md.
#
#   Rscript --vanilla Prism/docs/scripts/spotCheckSvCohort.R OUT_DIR DB_DIR TRACK_DIR HMF_DIR [N] [N_BAM]
#
# N events are sampled (default 300; fixed seed), N_BAM of their
# breakends (default 40) for the BAM recount. Needs bcftools (BCFTOOLS,
# default below) and samtools on PATH. Checks, each over the sampled
# breakends or events:
#   region_bed        every GIAB and ENCODE track column: the position is
#                     in a BED interval (start < POS <= end)
#   region_labels     ENCODE blacklist reason and RepeatMasker name, class
#                     and family: the set of labels of the intervals
#                     containing POS
#   region_hmf        Hartwig sv_prep_blacklist (.bed, 0-based) and
#                     unmap_regions (1-based)
#   centromere        bases strictly between POS and the nearest centromere
#                     model interval (0 inside or adjacent)
#   xs_breakend       xs_<set>_groups_exact, _w50, _group_list_w50 and
#                     own_germline_*: brute force over every breakend of
#                     the sets (bcftools query of the VCFs)
#   xs_pair           xs_pair_<set>_*: brute force over MATEID-linked pairs
#   bam_own           depth, clip_n and split_n of the own tumor, from
#                     samtools view (SAM text parsed here)
#   bam_others        bam_other_* summaries recomputed from the per-BAM
#                     files
#   pop_breakend      <source>_w100_n and _max_af: brute force over sites
#   pop_event         <source>_bp100_* and _ro50_*: brute force over sites
# The run fails if any value differs; the differences are printed.

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
})
options(width = 200)
Sys.setenv(VROOM_CONNECTION_SIZE = 2^26)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: spotCheckSvCohort.R OUT_DIR DB_DIR TRACK_DIR HMF_DIR [N] [N_BAM]", call. = FALSE)
out_dir <- path_abs(args[1])
db_dir <- path_abs(args[2])
track_dir <- path_abs(args[3])
hmf_dir <- path_abs(args[4])
n_events <- if (length(args) >= 5) as.integer(args[5]) else 300L
n_bam <- if (length(args) >= 6) as.integer(args[6]) else 40L
BCFTOOLS <- Sys.getenv("BCFTOOLS", "/admin/software/bcftools/bcftools-1.20/bin/bcftools")
SEP <- "|"
XS_WINDOW <- 50
POP_WINDOW <- 100
set.seed(20261122)

read_chr <- function(file, ...) {
  read_tsv(file, col_types = cols(.default = col_character()), quote = "", na = character(),
           show_col_types = FALSE, lazy = FALSE, trim_ws = FALSE, ...)
}
log_msg <- function(...) message(format(Sys.time(), "%H:%M:%S"), " ", ...)
blank <- function(x) {
  x <- as.character(x)
  if_else(is.na(x) | x == "NA", "", x)
}
set_str <- function(x) map_chr(str_split(blank(x), fixed(SEP)), \(v) str_c(sort(unique(v[v != ""])), collapse = SEP))
num_eq <- function(a, b) {
  a <- suppressWarnings(as.numeric(blank(a)))
  b <- suppressWarnings(as.numeric(blank(b)))
  (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= 1e-6 * pmax(1, abs(a), abs(b)))
}
results <- list()
record <- function(check, column, ok, detail) {
  results[[length(results) + 1]] <<- tibble(check, column, n = length(ok), n_equal = sum(ok))
  if (any(!ok)) {
    message("DIFF ", check, " ", column, ": ", sum(!ok), " of ", length(ok))
    print(head(detail[!ok, , drop = FALSE], 5))
  }
}

message("reading the table")
tab <- read_chr(path(db_dir, "sv_cohort_filtered.tsv.gz"))
samples <- read_chr(path(db_dir, "samples.tsv"))
ev <- tab |> slice_sample(n = min(n_events, nrow(tab)))
pts <- bind_rows(
  ev |> transmute(event_key, group_id, side = "start", id = start_id, chrom = start_chrom,
                  pos = as.integer(start_pos), orient = start_orient),
  ev |> filter(n_breakends == "2") |>
    transmute(event_key, group_id, side = "end", id = end_id, chrom = end_chrom,
              pos = as.integer(end_pos), orient = end_orient)
)
tab_val <- function(name) {
  vals <- character(nrow(pts))
  for (s in c("start", "end")) {
    i <- pts$side == s
    vals[i] <- ev[[str_c(s, "_", name)]][match(pts$event_key[i], ev$event_key)]
  }
  vals
}
message(nrow(ev), " events, ", nrow(pts), " breakends sampled")

# ---- region tracks -------------------------------------------------------

tracks <- read_chr(path(track_dir, "tracks.tsv"))
# B: a BED table with X2 and X3 already integer (converted once per table).
in_bed <- function(b, chrom, pos) {
  map2_lgl(chrom, pos, \(c, p) any(b$X1 == c & b$X2 < p & b$X3 >= p))
}
labels_at <- function(b, chrom, pos, label_col) {
  map2_chr(chrom, pos, \(c, p) {
    hit <- b[[label_col]][b$X1 == c & b$X2 < p & b$X3 >= p]
    str_c(sort(unique(hit)), collapse = SEP)
  })
}
read_bed <- function(file) read_tsv(file, col_names = FALSE, col_types = cols(.default = "c"), show_col_types = FALSE)
log_msg("region tracks")
for (name in tracks$name[str_starts(tracks$name, "giab_") | tracks$name == "encode_blacklist"]) {
  b <- read_bed(path(track_dir, tracks$file[tracks$name == name]))
  b <- b |> filter(X1 %in% pts$chrom) |> mutate(X2 = as.integer(X2), X3 = as.integer(X3))
  got <- in_bed(b, pts$chrom, pts$pos)
  want <- as.logical(tab_val(str_c("region_", name)))
  record("region_bed", name, got == want, tibble(pts, got, want))
  if (name == "encode_blacklist") {
    got_l <- labels_at(b, pts$chrom, pts$pos, "X4")
    want_l <- set_str(tab_val("region_encode_blacklist_reason"))
    record("region_labels", "encode_blacklist_reason", got_l == want_l, tibble(pts, got_l, want_l))
  }
}
rmsk <- read_tsv(path(track_dir, tracks$file[tracks$name == "ucsc_rmsk"]), col_names = FALSE,
                 col_select = c(6, 7, 8, 11, 12, 13), col_types = cols(.default = "c"), show_col_types = FALSE) |>
  filter(X6 %in% pts$chrom) |>
  mutate(X7 = as.integer(X7), X8 = as.integer(X8))
for (f in c("name", "class", "family")) {
  src <- c(name = "X11", class = "X12", family = "X13")[[f]]
  got <- map2_chr(pts$chrom, pts$pos, \(c, p) {
    hit <- rmsk[[src]][rmsk$X6 == c & rmsk$X7 < p & rmsk$X8 >= p]
    str_c(sort(unique(hit)), collapse = SEP)
  })
  want <- set_str(tab_val(str_c("region_rmsk_", f)))
  record("region_labels", str_c("rmsk_", f), got == want, tibble(pts, got, want))
}
rm(rmsk)

hmf_regions <- function(file) {
  b <- read_chr(file)
  shift <- if (str_detect(file, "\\.bed(\\.gz)?$")) 1L else 0L
  tibble(chrom = b$Chromosome, s = as.integer(b$PosStart) + shift, e = as.integer(b$PosEnd))
}
for (h in list(
  list(col = "region_hmf_sv_prep_blacklist", file = path(hmf_dir, "dna", "sv", "sv_prep_blacklist.38.bed")),
  list(col = "region_hmf_unmap", file = path(hmf_dir, "common", "unmap_regions.38.tsv"))
)) {
  r <- hmf_regions(h$file)
  got <- map2_lgl(pts$chrom, pts$pos, \(c, p) any(r$chrom == c & r$s <= p & r$e >= p))
  want <- as.logical(tab_val(h$col))
  record("region_hmf", h$col, got == want, tibble(pts, got, want))
}

log_msg("centromeres")
cen <- read_bed(path(track_dir, tracks$file[tracks$name == "ucsc_centromeres"])) |>
  transmute(chrom = X2, s = as.integer(X3) + 1L, e = as.integer(X4))
got <- map2_int(pts$chrom, pts$pos, \(c, p) {
  m <- cen |> filter(chrom == c)
  if (nrow(m) == 0) return(NA_integer_)
  as.integer(min(pmax(0L, m$s - p - 1L, p - m$e - 1L)))
})
want <- tab_val("region_centromere_distance")
record("centromere", "region_centromere_distance", num_eq(got, want), tibble(pts, got, want))

# ---- cross-sample sets ---------------------------------------------------
log_msg("cross-sample sets")

query <- function(vcf, group_id) {
  out <- system2(BCFTOOLS, c("query", "-f", shQuote("%ID\\t%CHROM\\t%POS\\t%ALT\\t%FILTER\\t%INFO/MATEID\\n"), vcf),
                 stdout = TRUE)
  if (!is.null(attr(out, "status"))) stop("bcftools failed on ", vcf, call. = FALSE)
  if (length(out) == 0) return(tibble())
  m <- str_split_fixed(out, "\t", 6)
  tibble(group = group_id, id = m[, 1], chrom = m[, 2], pos = as.integer(m[, 3]),
         orient = if_else(str_detect(m[, 4], "^[ACGTNacgtn]"), "1", "-1"), filter = m[, 5],
         mate = na_if(m[, 6], "."))
}
has_normal <- !is.na(samples$normal_id) & samples$normal_id != ""
sets <- list(
  somatic = map2(samples$group_id, samples$sample_id, \(g, s) {
    query(path(out_dir, g, "purple", str_c(s, ".purple.sv.vcf.gz")), g)
  }) |> list_rbind() |> filter(filter != "INFERRED"),
  germline = map2(samples$group_id[has_normal], samples$sample_id[has_normal], \(g, s) {
    query(path(out_dir, g, "purple", str_c(s, ".purple.sv.germline.vcf.gz")), g)
  }) |> list_rbind(),
  unfiltered = map2(samples$group_id, samples$sample_id, \(g, s) {
    query(path(out_dir, g, "esvee", str_c(s, ".esvee.unfiltered.vcf.gz")), g)
  }) |> list_rbind()
)
near_hits <- function(set, i) {
  set |> filter(chrom == pts$chrom[i], orient == pts$orient[i], abs(pos - pts$pos[i]) <= XS_WINDOW)
}
for (n in names(sets)) {
  s <- sets[[n]]
  h <- map(seq_len(nrow(pts)), \(i) near_hits(s, i))
  oth <- map2(h, pts$group_id, \(x, g) x |> filter(group != g))
  got_exact <- map2_int(oth, seq_len(nrow(pts)), \(x, i) n_distinct(x$group[x$pos == pts$pos[i]]))
  got_near <- map_int(oth, \(x) n_distinct(x$group))
  got_list <- map_chr(oth, \(x) str_c(sort(unique(x$group)), collapse = SEP))
  record("xs_breakend", str_c("xs_", n, "_groups_exact"), num_eq(got_exact, tab_val(str_c("xs_", n, "_groups_exact"))),
         tibble(pts, got_exact, want = tab_val(str_c("xs_", n, "_groups_exact"))))
  w_col <- str_c("xs_", n, "_groups_w", XS_WINDOW)
  record("xs_breakend", w_col, num_eq(got_near, tab_val(w_col)), tibble(pts, got_near, want = tab_val(w_col)))
  l_col <- str_c("xs_", n, "_group_list_w", XS_WINDOW)
  record("xs_breakend", l_col, got_list == set_str(tab_val(l_col)), tibble(pts, got_list, want = tab_val(l_col)))
  if (n == "germline") {
    ownh <- map2(h, pts$group_id, \(x, g) x |> filter(group == g))
    has_g <- pts$group_id %in% samples$group_id[has_normal]
    got_oe <- if_else(has_g, map2_int(ownh, seq_len(nrow(pts)), \(x, i) sum(x$pos == pts$pos[i])), NA_integer_)
    got_on <- if_else(has_g, map_int(ownh, nrow), NA_integer_)
    record("xs_breakend", "own_germline_exact", num_eq(got_oe, tab_val("own_germline_exact")),
           tibble(pts, got_oe, want = tab_val("own_germline_exact")))
    on_col <- str_c("own_germline_w", XS_WINDOW)
    record("xs_breakend", on_col, num_eq(got_on, tab_val(on_col)), tibble(pts, got_on, want = tab_val(on_col)))
  }
  # Pairs: another group's record near the start whose MATEID record is near the end.
  pe <- ev |> filter(n_breakends == "2")
  got_pair <- map_chr(seq_len(nrow(pe)), \(j) {
    a <- s |> filter(group != pe$group_id[j], chrom == pe$start_chrom[j], orient == pe$start_orient[j],
                     abs(pos - as.integer(pe$start_pos[j])) <= XS_WINDOW, !is.na(mate))
    b <- s |> filter(group != pe$group_id[j], chrom == pe$end_chrom[j], orient == pe$end_orient[j],
                     abs(pos - as.integer(pe$end_pos[j])) <= XS_WINDOW)
    g <- inner_join(a, b, by = c("group", mate = "id")) |> pull(group)
    str_c(sort(unique(g)), collapse = SEP)
  })
  pn <- str_c("xs_pair_", n, "_groups_w", XS_WINDOW)
  pl <- str_c("xs_pair_", n, "_group_list_w", XS_WINDOW)
  got_pn <- map_int(str_split(got_pair, fixed(SEP)), \(v) sum(v != ""))
  record("xs_pair", pn, num_eq(got_pn, pe[[pn]]), tibble(pe |> select(event_key), got_pn, want = pe[[pn]]))
  record("xs_pair", pl, got_pair == set_str(pe[[pl]]), tibble(pe |> select(event_key), got_pair, want = pe[[pl]]))
}

# ---- BAM evidence --------------------------------------------------------
log_msg("BAM evidence")

bam_dir <- path(db_dir, "bam_evidence")
bam_files <- tibble(
  sample = c(samples$sample_id, samples$normal_id[has_normal]),
  group = c(samples$group_id, samples$group_id[has_normal]),
  role = c(rep("tumor", nrow(samples)), rep("normal", sum(has_normal)))
)
bam <- pmap(bam_files, \(sample, group, role) {
  read_chr(path(bam_dir, str_c(sample, ".bam_evidence.tsv.gz"))) |>
    semi_join(pts, by = c("group_id", "id")) |>
    mutate(bam_sample = sample, bam_group = group, bam_role = role)
}) |> list_rbind()
for (r in c("normal", "tumor")) {
  label <- str_c("other_", r, "s")
  o <- bam |> filter(bam_role == r, bam_group != group_id) |>
    mutate(across(c(clip_n, split_n, disc_n), \(x) suppressWarnings(as.integer(x)))) |>
    summarize(n = n(), clip_samples = sum(clip_n > 0), split_samples = sum(split_n > 0, na.rm = TRUE),
              disc_samples = sum(disc_n > 0, na.rm = TRUE), split_reads = sum(split_n, na.rm = TRUE),
              disc_reads = sum(disc_n, na.rm = TRUE), no_split = all(is.na(split_n)),
              no_disc = all(is.na(disc_n)), .by = c(group_id, id))
  o <- pts |> left_join(o, by = c("group_id", "id"))
  for (f in c("n", "clip_samples", "split_samples", "disc_samples", "split_reads", "disc_reads")) {
    got <- o[[f]]
    if (f %in% c("split_samples", "split_reads")) got[o$no_split] <- NA
    if (f %in% c("disc_samples", "disc_reads")) got[o$no_disc] <- NA
    col <- str_c("bam_", label, "_", f)
    record("bam_others", col, num_eq(got, tab_val(col)), tibble(pts, got, want = tab_val(col)))
  }
}

# Own tumor: depth, clip_n, split_n from samtools view, parsed here.
cigar_ops <- function(cigar) {
  m <- str_match_all(cigar, "(\\d+)([MIDNSHP=X])")[[1]]
  tibble(n = as.integer(m[, 2]), op = m[, 3])
}
ref_len <- function(cigar) {
  o <- cigar_ops(cigar)
  sum(o$n[o$op %in% c("M", "D", "N", "=", "X")])
}
clips <- function(cigar) {
  o <- cigar_ops(cigar)
  lead <- o$op %in% c("S", "H")
  left <- sum(o$n[cumprod(lead) == 1])
  right <- sum(rev(o$n)[cumprod(rev(lead)) == 1])
  c(left, right)
}
vcf_info <- map2(samples$group_id, samples$sample_id, \(g, s) {
  out <- system2(BCFTOOLS, c("query", "-f", shQuote("%ID\\t%INFO/CIPOS\\n"),
                             path(out_dir, g, "purple", str_c(s, ".purple.sv.vcf.gz"))), stdout = TRUE)
  m <- str_split_fixed(out, "\t", 2)
  tibble(group_id = g, id = m[, 1], cipos = m[, 2])
}) |> list_rbind()
bp <- pts |> slice_sample(n = min(n_bam, nrow(pts))) |>
  left_join(vcf_info, by = c("group_id", "id")) |>
  left_join(samples |> select(group_id, sample_id), by = "group_id")
mates <- bind_rows(
  # One-breakend events: end_id is "" (no mate), read as NA.
  tab |> transmute(group_id, id = start_id, mchrom = end_chrom, mpos = as.integer(end_pos), morient = end_orient,
                   mid = na_if(end_id, "")),
  tab |> filter(n_breakends == "2") |>
    transmute(group_id, id = end_id, mchrom = start_chrom, mpos = as.integer(start_pos), morient = start_orient, mid = start_id)
)
vcf_info_m <- vcf_info |> rename(mid = id, mcipos = cipos)
bp <- bp |> left_join(mates, by = c("group_id", "id")) |> left_join(vcf_info_m, by = c("group_id", "mid"))
ci <- function(x) {
  v <- as.integer(str_split_fixed(if_else(is.na(x) | x == ".", "0,0", x), ",", 2))
  v
}
got <- pmap(bp, function(...) {
  b <- list(...)
  bamf <- path(out_dir, b$group_id, "alignments", b$sample_id, str_c(b$sample_id, ".redux.bam"))
  cip <- ci(b$cipos)
  lo <- b$pos + cip[1] - 5L
  hi <- b$pos + cip[2] + 5L
  region <- str_c(b$chrom, ":", max(1L, min(lo, b$pos) - 1L), "-", max(hi, b$pos) + 1L)
  # Skip unmapped, secondary, QC-fail and duplicate records (0x704).
  sam <- system2("samtools", c("view", "-F", "0x704", bamf, region), stdout = TRUE)
  f <- str_split_fixed(sam, "\t", 12)
  flag <- as.integer(f[, 2])
  start <- as.integer(f[, 4])
  ends <- start + map_int(f[, 6], ref_len) - 1L
  primary <- bitwAnd(flag, 0x800) == 0
  depth <- sum(primary & start <= b$pos & ends >= b$pos)
  cl <- map(f[, 6], clips)
  left <- map_int(cl, 1)
  right <- map_int(cl, 2)
  at <- if (b$orient == "1") right >= 5 & ends >= lo & ends <= hi else left >= 5 & start >= lo & start <= hi
  names_at <- unique(f[at, 1])
  split <- NA_integer_
  if (!is.na(b$mid)) {
    mci <- ci(b$mcipos)
    mlo <- b$mpos + mci[1] - 5L
    mhi <- b$mpos + mci[2] + 5L
    sa <- str_match(f[at, 12], "SA:Z:([^\t]+)")[, 2]
    hit <- map_lgl(sa, \(x) {
      if (is.na(x)) return(FALSE)
      any(map_lgl(str_split(str_remove(x, ";$"), ";")[[1]], \(e) {
        p <- str_split(e, ",")[[1]]
        if (p[1] != b$mchrom) return(FALSE)
        s <- as.integer(p[2])
        edge <- if (b$morient == "1") s + ref_len(p[4]) - 1L else s
        edge >= mlo && edge <= mhi
      }))
    })
    split <- length(unique(f[at, 1][hit]))
  }
  tibble(depth = depth, clip_n = length(names_at), split_n = split)
}) |> list_rbind()
for (f in c("depth", "clip_n", "split_n")) {
  want <- map2_chr(bp$side, seq_len(nrow(bp)), \(s, i) ev[[str_c(s, "_bam_tumor_", f)]][match(bp$event_key[i], ev$event_key)])
  record("bam_own", str_c("bam_tumor_", f), num_eq(got[[f]], want), tibble(bp |> select(event_key, side, chrom, pos), got = got[[f]], want))
}

# ---- population SV call sets ---------------------------------------------
log_msg("population SV call sets")

POP_TYPES <- list(DEL = c("DEL", "CNV", "BND"), DUP = c("DUP", "CNV", "BND"), INV = c("INV", "BND"),
                  INS = "INS", BND = c("BND", "CTX"))
POP_RO_TYPES <- list(DEL = c("DEL", "CNV"), DUP = c("DUP", "CNV"), INV = "INV")
chroms <- unique(c(pts$chrom))
for (src in c("gnomad_sv", "kgp_sv", "dbvar_nstd166", "dbvar_nstd200")) {
  s <- read_chr(path(track_dir, src, "sites.tsv.gz")) |>
    mutate(
      pos = as.integer(pos), endn = suppressWarnings(as.integer(end)),
      c2 = if_else(chr2 == ".", chrom, chr2),
      p2 = suppressWarnings(as.integer(if_else(end2 != ".", end2, end))),
      pass = filter %in% c(".", "PASS"), afn = suppressWarnings(as.numeric(af))
    ) |>
    filter(chrom %in% chroms | c2 %in% chroms)
  # Breakend: any breakpoint within the window.
  pt <- map2(pts$chrom, pts$pos, \(c, p) {
    s |> filter((chrom == c & abs(pos - p) <= POP_WINDOW) | (c2 == c & !is.na(p2) & abs(p2 - p) <= POP_WINDOW))
  })
  got_n <- map_int(pt, nrow)
  got_af <- map_dbl(pt, \(x) if (any(x$pass & !is.na(x$afn))) max(x$afn[x$pass], na.rm = TRUE) else NA_real_)
  cn <- str_c(src, "_w", POP_WINDOW, "_n")
  ca <- str_c(src, "_w", POP_WINDOW, "_max_af")
  record("pop_breakend", cn, num_eq(got_n, tab_val(cn)), tibble(pts, got_n, want = tab_val(cn)))
  record("pop_breakend", ca, num_eq(got_af, tab_val(ca)), tibble(pts, got_af, want = tab_val(ca)))
  # Events.
  summarize_m <- function(m) {
    list(n = nrow(m), n_pass = sum(m$pass), ids = str_c(sort(m$id), collapse = SEP),
         max_af = if (any(m$pass & !is.na(m$afn))) max(m$afn[m$pass], na.rm = TRUE) else NA_real_)
  }
  pe <- ev |> filter(n_breakends == "2")
  bp_got <- map(seq_len(nrow(pe)), \(j) {
    sc <- pe$start_chrom[j]; spos <- as.integer(pe$start_pos[j])
    ec <- pe$end_chrom[j]; epos <- as.integer(pe$end_pos[j])
    t <- pe$info_SVTYPE[j]
    m <- s |> filter(
      ((chrom == sc & abs(pos - spos) <= POP_WINDOW & c2 == ec & !is.na(p2) & abs(p2 - epos) <= POP_WINDOW) |
         (chrom == ec & abs(pos - epos) <= POP_WINDOW & c2 == sc & !is.na(p2) & abs(p2 - spos) <= POP_WINDOW)),
      svtype %in% (POP_TYPES[[t]] %||% character())
    )
    summarize_m(m)
  })
  ro_rows <- pe$info_SVTYPE %in% names(POP_RO_TYPES) & pe$start_chrom == pe$end_chrom
  ro_got <- map(seq_len(nrow(pe)), \(j) {
    if (!ro_rows[j]) return(list(n = NA, n_pass = NA, ids = "", max_af = NA))
    lo <- min(as.integer(pe$start_pos[j]), as.integer(pe$end_pos[j]))
    hi <- max(as.integer(pe$start_pos[j]), as.integer(pe$end_pos[j]))
    m <- s |>
      filter(chrom == pe$start_chrom[j], c2 == chrom, !is.na(endn), svtype %in% POP_RO_TYPES[[pe$info_SVTYPE[j]]]) |>
      mutate(slo = pmin(pos, endn), shi = pmax(pos, endn),
             ov = pmin(hi, shi) - pmax(lo, slo) + 1, longer = pmax(hi - lo + 1, shi - slo + 1)) |>
      filter(ov > 0, ov >= 0.5 * longer)
    summarize_m(m)
  })
  for (rule in c(str_c("bp", POP_WINDOW), "ro50")) {
    g <- if (rule == "ro50") ro_got else bp_got
    for (f in c("n", "n_pass", "ids", "max_af")) {
      col <- str_c(src, "_", rule, "_", f)
      got <- map(g, f) |> unlist()
      ok <- if (f == "ids") blank(as.character(got)) == set_str(pe[[col]]) else num_eq(got, pe[[col]])
      record("pop_event", col, ok, tibble(pe |> select(event_key, info_SVTYPE), got, want = pe[[col]]))
    }
  }
  rm(s)
}

# ---- summary -------------------------------------------------------------

res <- list_rbind(results)
summary <- res |> summarize(columns = n(), values = sum(n), equal = sum(n_equal), .by = check)
print(summary, n = Inf)
if (any(res$n_equal != res$n)) {
  message("FAILED: ", sum(res$n - res$n_equal), " values differ")
  quit(status = 1)
}
message("all spot checks passed: ", sum(res$n), " values in ", nrow(res), " columns")
