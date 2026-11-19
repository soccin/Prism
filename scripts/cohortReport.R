#!/usr/bin/env Rscript
# Prism: cohort report from a finished oncoanalyser run. Called by
# cohortReport.sh; see that script for the usage and the output layout.
#
#   Rscript cohortReport.R --out-dir OUT --res-dir RES [--manifest CSV]
#       [--genes TXT] [--panel TSV] [--version STR]
#
# Reads the per-sample result tables under OUT and writes RES/prism_report.xlsx
# plus one CSV per table in RES/data. The report counts and lists events; it
# applies no significance test and no filter on sample class.
#
# --genes: the genes of interest, one symbol per line. Without it they are
#   the genes with a pipeline driver call of likelihood 0.8 or higher in at
#   least one sample, and that list is written to
#   RES/genes_of_interest_default.txt.

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
  library(openxlsx)
})

# SV VCF lines (assembled insert sequences) are longer than readr's
# default buffer for compressed files.
Sys.setenv(VROOM_CONNECTION_SIZE = 2^24)

# Thresholds (PURPLE's driver rules, applied to every gene).
AMP_PLOIDY_RATIO <- 3
DEEP_DEL_MAX_CN <- 0.5
LOH_MAX_MINOR_CN <- 0.5
CODING_EFFECTS <- c("MISSENSE", "NONSENSE_OR_FRAMESHIFT", "SPLICE")

# Without --genes: the lowest driver likelihood that puts a gene on the
# default list of genes of interest.
DEFAULT_GOI_MIN_LIKELIHOOD <- 0.8

# ---- arguments ---------------------------------------------------------

parse_args <- function(args) {
  opts <- list(
    out_dir = NULL, res_dir = NULL, manifest = NULL, genes = NULL,
    panel = NULL, version = "unknown"
  )
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

# Every table is read as text; numeric columns are converted where used.
read_chr <- function(file, ...) {
  read_tsv(
    file, col_types = cols(.default = col_character()), quote = "",
    na = character(), show_col_types = FALSE, lazy = FALSE, ...
  )
}

# Read FILE for each sample (FILE is a column of TBL), skipping samples
# without it, and bind the rows with sample_name in front.
read_per_sample <- function(tbl, file_col, reader = read_chr) {
  tbl |>
    filter(file_exists(.data[[file_col]])) |>
    select(sample_name, file = all_of(file_col)) |>
    pmap(\(sample_name, file) {
      reader(file) |> mutate(sample_name = sample_name, .before = 1)
    }) |>
    list_rbind()
}

first_file <- function(dir, glob) {
  if (!dir_exists(dir)) {
    return(NA_character_)
  }
  found <- dir_ls(dir, glob = glob, type = "file")
  if (length(found) == 0) NA_character_ else as.character(found[1])
}

read_vcf_body <- function(file) {
  n_skip <- sum(str_starts(read_lines(file, n_max = 20000), "##"))
  read_chr(file, skip = n_skip, comment = "") |>
    rename(CHROM = `#CHROM`)
}

info_get <- function(info, key) {
  str_match(info, str_glue("(?:^|;){key}=([^;]+)"))[, 2]
}

info_flag <- function(info, key) {
  str_detect(info, str_glue("(?:^|;){key}(?:;|$)"))
}

fmt_get <- function(format, value, key) {
  map2_chr(str_split(format, ":"), str_split(value, ":"), \(k, v) {
    i <- match(key, k)
    if (is.na(i) || i > length(v)) NA_character_ else v[i]
  })
}

parse_coords <- function(x) {
  m <- str_match(x, "^([^:]+):(\\d+):(-?1)$")
  tibble(chrom = m[, 2], pos = as.numeric(m[, 3]), orient = m[, 4])
}

is_xy_chr <- function(x) {
  str_remove(x, "^chr") %in% c("X", "Y")
}

pct <- function(n, total) {
  if_else(total > 0, round(100 * n / total, 1), NA_real_)
}

collapse_unique <- function(x, sep = ";") {
  x <- sort(unique(x[!is.na(x) & x != ""]))
  if (length(x) == 0) NA_character_ else str_c(x, collapse = sep)
}

# ---- samples -----------------------------------------------------------

group_dirs <- dir_ls(out_dir, type = "directory")
samples <- tibble(group_dir = as.character(group_dirs)) |>
  mutate(
    group = path_file(group_dir),
    purity_file = map_chr(
      group_dir, \(d) first_file(path(d, "purple"), "*.purple.purity.tsv")
    )
  ) |>
  filter(!is.na(purity_file)) |>
  mutate(
    sample_name = purity_file |> path_file() |> str_remove("\\.purple\\.purity\\.tsv$"),
    purple_qc_file = path(group_dir, "purple", str_c(sample_name, ".purple.qc")),
    somatic_vcf = path(group_dir, "purple", str_c(sample_name, ".purple.somatic.vcf.gz")),
    sv_vcf = path(group_dir, "purple", str_c(sample_name, ".purple.sv.vcf.gz")),
    gene_cn_file = path(group_dir, "purple", str_c(sample_name, ".purple.cnv.gene.tsv")),
    arm_file = path(group_dir, "purple", str_c(sample_name, ".purple.chromosome_arm.tsv")),
    linx_dir = path(group_dir, "linx", "somatic_annotations"),
    svs_file = path(linx_dir, str_c(sample_name, ".linx.svs.tsv")),
    clusters_file = path(linx_dir, str_c(sample_name, ".linx.clusters.tsv")),
    breakend_file = path(linx_dir, str_c(sample_name, ".linx.breakend.tsv")),
    fusion_file = path(linx_dir, str_c(sample_name, ".linx.fusion.tsv")),
    driver_file = path(linx_dir, str_c(sample_name, ".linx.driver.catalog.tsv")),
    virus_file = path(group_dir, "virusinterpreter", str_c(sample_name, ".virus.annotated.tsv")),
    chord_file = path(group_dir, "chord", str_c(sample_name, ".chord.prediction.tsv")),
    sigs_file = path(group_dir, "sigs", str_c(sample_name, ".sig.allocation.tsv")),
    lilac_file = path(group_dir, "lilac", str_c(sample_name, ".lilac.tsv")),
    peach_file = map_chr(
      group_dir, \(d) first_file(path(d, "peach"), "*.peach.haplotypes.best.tsv")
    )
  )
if (nrow(samples) == 0) {
  stop(str_glue("No groups with a PURPLE purity file in {out_dir}"), call. = FALSE)
}
log_msg("{nrow(samples)} samples in {out_dir}")

purity <- read_per_sample(samples, "purity_file")
samples <- samples |>
  left_join(
    purity |> transmute(sample_name, matched_normal = runMode == "TUMOR_GERMLINE"),
    by = "sample_name"
  )

# Manifest: sample_class and any other sample-level column.
manifest_cols <- character()
if (!is.null(opts$manifest)) {
  manifest <- read_csv(
    opts$manifest, col_types = cols(.default = col_character()),
    show_col_types = FALSE
  )
  if (!all(c("sample_id", "sample_class") %in% names(manifest))) {
    stop("The manifest needs the columns sample_id and sample_class", call. = FALSE)
  }
  manifest <- manifest |>
    rename(sample_name = sample_id) |>
    rename_with(str_to_lower, .cols = matches("^sex$", ignore.case = TRUE)) |>
    filter(sample_name %in% samples$sample_name)
  missing <- setdiff(samples$sample_name, manifest$sample_name)
  if (length(missing) > 0) {
    stop(
      "Samples not in the manifest: ", str_c(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (any(duplicated(manifest$sample_name))) {
    stop("The manifest lists a sample more than once", call. = FALSE)
  }
  if (any(is.na(manifest$sample_class) | manifest$sample_class == "")) {
    stop("The manifest has samples without a sample_class", call. = FALSE)
  }
  manifest_cols <- setdiff(names(manifest), c("sample_name", "sample_class"))
  samples <- samples |> left_join(manifest, by = "sample_name")
} else {
  samples <- samples |> mutate(sample_class = "Tumor")
}

# X/Y events are kept only for a sample with a sex in the manifest and a
# matched normal; for every other sample they are removed from all tables.
if ("sex" %in% names(samples)) {
  samples <- samples |>
    mutate(xy_reported = matched_normal & !is.na(sex) & sex != "")
} else {
  samples <- samples |> mutate(xy_reported = FALSE)
}
xy_ok_samples <- samples |> filter(xy_reported) |> pull(sample_name)
redacted <- list()

redact_xy <- function(df, chr_cols, label) {
  hit <- map(chr_cols, \(col) is_xy_chr(df[[col]])) |> reduce(`|`)
  drop <- hit & !(df$sample_name %in% xy_ok_samples)
  redacted[[label]] <<- sum(drop)
  df[!drop, ]
}

sample_key <- samples |> select(sample_name, sample_class, matched_normal)
classes <- sort(unique(samples$sample_class))

add_class <- function(df) {
  df |>
    left_join(sample_key, by = "sample_name") |>
    relocate(sample_class, matched_normal, .after = sample_name)
}

# ---- reference: driver gene panel, genes of interest -------------------

panel <- tibble(gene = character(), panel_role = character())
if (!is.null(opts$panel) && file_exists(opts$panel)) {
  panel <- read_chr(opts$panel) |> transmute(gene, panel_role = likelihoodType)
} else {
  log_msg("No driver gene panel; panel_role is left empty")
}

# Without a gene list the genes of interest are set from the pipeline's
# driver calls, once those are read.
genes_given <- !is.null(opts$genes)
genes_of_interest <- character()
if (genes_given) {
  genes_of_interest <- read_lines(opts$genes) |>
    str_remove("#.*$") |>
    str_trim() |>
    discard(\(g) g == "") |>
    unique()
  log_msg("{length(genes_of_interest)} genes of interest")
}

# ---- gene copy number --------------------------------------------------

log_msg("gene copy number")
ploidy_tbl <- purity |> transmute(sample_name, ploidy = as.numeric(ploidy))

gene_cn_all <- read_per_sample(samples, "gene_cn_file") |>
  filter(isCanonical == "true") |>
  mutate(across(
    c(start, end, minCopyNumber, maxCopyNumber, minMinorAlleleCopyNumber),
    as.numeric
  )) |>
  left_join(ploidy_tbl, by = "sample_name")

gene_info <- gene_cn_all |>
  distinct(gene, .keep_all = TRUE) |>
  transmute(
    gene, chromosome, band = chromosomeBand, gene_start = start,
    gene_end = end, gene_length_kb = round((end - start + 1) / 1000, 1),
    is_xy = is_xy_chr(chromosome)
  )

unknown_genes <- setdiff(genes_of_interest, gene_info$gene)
if (length(unknown_genes) > 0) {
  stop(
    "Genes of interest not in the PURPLE gene table: ",
    str_c(unknown_genes, collapse = ", "), call. = FALSE
  )
}

gene_cn <- gene_cn_all |>
  mutate(
    cn_event = case_when(
      minCopyNumber < DEEP_DEL_MAX_CN ~ "DEEP_DEL",
      minCopyNumber > AMP_PLOIDY_RATIO * ploidy ~ "AMP",
      .default = NA_character_
    ),
    loh = minCopyNumber >= DEEP_DEL_MAX_CN &
      minMinorAlleleCopyNumber < LOH_MAX_MINOR_CN
  ) |>
  redact_xy("chromosome", "gene_cn") |>
  transmute(
    sample_name, gene, chromosome, band = chromosomeBand, start, end,
    cn_event, loh,
    min_copy_number = minCopyNumber, max_copy_number = maxCopyNumber,
    min_minor_allele_copy_number = minMinorAlleleCopyNumber,
    sample_ploidy = ploidy,
    min_copy_number_over_ploidy = round(minCopyNumber / ploidy, 3),
    somatic_regions = somaticRegions,
    min_region_start_support = minRegionStartSupport,
    min_region_end_support = minRegionEndSupport,
    min_region_method = minRegionMethod,
    depth_window_count = depthWindowCount
  )
rm(gene_cn_all)
gene_cn_events <- gene_cn |> filter(!is.na(cn_event) | loh)

# ---- somatic variants --------------------------------------------------

log_msg("somatic variants")
read_somatic_vcf <- function(file, sample_name) {
  vcf <- read_vcf_body(file)
  vcf |>
    filter(FILTER == "PASS") |>
    mutate(impact = info_get(INFO, "IMPACT")) |>
    filter(!is.na(impact)) |>
    separate_wider_delim(
      impact, ",",
      names = c(
        "gene", "transcript", "canonical_effect", "coding_effect",
        "splice_region", "hgvs_coding", "hgvs_protein", "other_effects",
        "worst_coding_effect", "genes_affected"
      ),
      too_few = "align_start", too_many = "merge"
    ) |>
    filter(coding_effect %in% CODING_EFFECTS) |>
    mutate(
      tumor = .data[[sample_name]],
      tumor_ad = fmt_get(FORMAT, tumor, "AD"),
      mutation_class = case_when(
        coding_effect == "NONSENSE_OR_FRAMESHIFT" ~ "truncating",
        coding_effect == "SPLICE" ~ "splice",
        str_detect(canonical_effect, "inframe") ~ "inframe",
        .default = "missense"
      )
    ) |>
    transmute(
      gene, chromosome = CHROM, position = as.numeric(POS), ref = REF,
      alt = ALT, mutation_class, canonical_effect, hgvs_coding, hgvs_protein,
      transcript,
      hotspot = info_flag(INFO, "HOTSPOT"),
      biallelic = info_flag(INFO, "BIALLELIC"),
      reported = info_flag(INFO, "REPORTED"),
      tier = info_get(INFO, "TIER"),
      qual = as.numeric(QUAL),
      tumor_depth = as.numeric(fmt_get(FORMAT, tumor, "DP")),
      tumor_alt_reads = as.numeric(str_extract(tumor_ad, "[^,]+$")),
      tumor_vaf = as.numeric(fmt_get(FORMAT, tumor, "AF")),
      purity_adjusted_vaf = as.numeric(info_get(INFO, "PURPLE_AF")),
      copy_number = as.numeric(info_get(INFO, "PURPLE_CN")),
      variant_copy_number = as.numeric(info_get(INFO, "PURPLE_VCN")),
      minor_allele_copy_number = as.numeric(info_get(INFO, "PURPLE_MACN")),
      subclonal_likelihood = as.numeric(info_get(INFO, "SUBCL")),
      mappability = as.numeric(info_get(INFO, "MAPPABILITY")),
      gnomad_freq = as.numeric(info_get(INFO, "GND_FREQ")),
      pon_count = as.numeric(info_get(INFO, "PON_COUNT"))
    )
}

variants <- samples |>
  filter(file_exists(somatic_vcf)) |>
  select(sample_name, somatic_vcf) |>
  pmap(\(sample_name, somatic_vcf) {
    log_msg("  {sample_name}")
    read_somatic_vcf(somatic_vcf, sample_name) |>
      mutate(sample_name = sample_name, .before = 1)
  }) |>
  list_rbind() |>
  redact_xy("chromosome", "somatic_variants")

# ---- structural variants ----------------------------------------------

log_msg("structural variants")
read_sv_vcf <- function(file, sample_name) {
  vcf <- read_vcf_body(file)
  tumor <- vcf[[sample_name]]
  vcf |>
    transmute(
      vcf_id = ID,
      sv_qual = suppressWarnings(as.numeric(QUAL)),
      sv_filter = FILTER,
      tumor_fragments = as.numeric(fmt_get(FORMAT, tumor, "VF")),
      tumor_split_fragments = as.numeric(fmt_get(FORMAT, tumor, "SF")),
      tumor_discordant_fragments = as.numeric(fmt_get(FORMAT, tumor, "DF")),
      tumor_ref_reads = as.numeric(fmt_get(FORMAT, tumor, "REF")),
      tumor_ref_pairs = as.numeric(fmt_get(FORMAT, tumor, "REFPAIR")),
      tumor_breakend_af = as.numeric(fmt_get(FORMAT, tumor, "AF")),
      sv_pon_count = as.numeric(info_get(INFO, "PON_COUNT"))
    )
}

sv_support <- samples |>
  filter(file_exists(sv_vcf)) |>
  select(sample_name, sv_vcf) |>
  pmap(\(sample_name, sv_vcf) {
    read_sv_vcf(sv_vcf, sample_name) |>
      mutate(sample_name = sample_name, .before = 1)
  }) |>
  list_rbind()

clusters <- read_per_sample(samples, "clusters_file") |>
  transmute(
    sample_name, clusterId, cluster_category = category,
    cluster_type = resolvedType, cluster_sv_count = as.numeric(clusterCount)
  )

svs <- read_per_sample(samples, "svs_file") |>
  left_join(clusters, by = c("sample_name", "clusterId"))

# One row per SV end.
sv_ends <- bind_rows(
  svs |> mutate(
    isStart = "true", vcf_id = vcfIdStart, coords = coordsStart,
    other_end = coordsEnd, fragile_site = fragileSiteStart,
    line_type = lineTypeStart
  ),
  svs |> filter(coordsEnd != "") |> mutate(
    isStart = "false", vcf_id = vcfIdEnd, coords = coordsEnd,
    other_end = coordsStart, fragile_site = fragileSiteEnd,
    line_type = lineTypeEnd
  )
) |>
  mutate(parse_coords(coords), other_chrom = parse_coords(other_end)$chrom) |>
  transmute(
    sample_name, svId, isStart, vcf_id, breakpoint = coords, chrom, pos,
    orient, other_end, other_chrom, sv_type = type,
    cluster_id = clusterId, cluster_category, cluster_type, cluster_sv_count,
    junction_cn_min = as.numeric(junctionCopyNumberMin),
    junction_cn_max = as.numeric(junctionCopyNumberMax),
    fragile_site, line_type, sv_genes_start = geneStart, sv_genes_end = geneEnd
  ) |>
  left_join(sv_support, by = c("sample_name", "vcf_id")) |>
  mutate(weak_support = cluster_category == "ARTIFACT")

sv_per_sample <- svs |>
  summarize(
    sv_total = n(),
    sv_weak_support = sum(cluster_category == "ARTIFACT"),
    .by = sample_name
  )

sv_ends <- sv_ends |> redact_xy(c("chrom", "other_chrom"), "sv_breakends")

# The same breakpoint (chromosome, position, orientation) in more than one
# sample. Breakends inferred from copy number (INF) have no base-pair
# position and are left out.
same_bp <- sv_ends |>
  filter(sv_type != "INF", !is.na(pos)) |>
  summarize(
    same_bp_n_samples = n_distinct(sample_name),
    same_bp_samples = collapse_unique(sample_name),
    .by = breakpoint
  )
sv_ends <- sv_ends |>
  left_join(same_bp, by = "breakpoint") |>
  mutate(
    same_bp_n_samples = replace_na(same_bp_n_samples, 1L),
    same_bp_flag = same_bp_n_samples > 1
  )

sv_same_breakpoint <- sv_ends |>
  filter(same_bp_flag) |>
  arrange(desc(same_bp_n_samples), chrom, pos, sample_name) |>
  select(
    breakpoint, same_bp_n_samples, same_bp_samples, sample_name, other_end,
    sv_type, cluster_category, cluster_type, cluster_sv_count,
    junction_cn_min, junction_cn_max, weak_support, sv_qual, sv_filter,
    tumor_fragments, tumor_split_fragments, tumor_discordant_fragments,
    tumor_ref_reads, tumor_ref_pairs, tumor_breakend_af, sv_pon_count,
    sv_genes_start, sv_genes_end
  )

# Gene disruptions: breakends LINX calls disruptive (canonical transcript),
# plus any breakend the pipeline reported.
breakends <- read_per_sample(samples, "breakend_file")
disruptions <- breakends |>
  filter(canonical == "true", disruptive == "true" | reportedStatus == "REPORTED") |>
  inner_join(
    sv_ends |> select(-vcf_id),
    by = c("sample_name", "svId", "isStart")
  ) |>
  left_join(gene_info |> select(gene, gene_length_kb), by = "gene") |>
  mutate(row_id = row_number())

nearest <- disruptions |>
  select(row_id, gene, sample_name, pos) |>
  inner_join(
    disruptions |> select(gene, other_sample = sample_name, other_pos = pos),
    by = "gene", relationship = "many-to-many"
  ) |>
  filter(sample_name != other_sample) |>
  summarize(nearest_other_sample_bp = min(abs(pos - other_pos)), .by = row_id)

disruptions <- disruptions |>
  left_join(nearest, by = "row_id") |>
  arrange(gene, pos, sample_name) |>
  transmute(
    sample_name, gene, gene_length_kb, breakpoint, other_end, sv_type,
    disruptive = disruptive == "true",
    pipeline_reported = reportedStatus == "REPORTED",
    weak_support, same_bp_flag, same_bp_n_samples, nearest_other_sample_bp,
    region = regionType, coding = codingType, exon_up = exonUp,
    exon_down = exonDown, gene_orientation = geneOrientation,
    undisrupted_copy_number = as.numeric(undisruptedCopyNumber),
    junction_cn_min, junction_cn_max, cluster_id, cluster_category,
    cluster_type, cluster_sv_count, sv_qual, sv_filter, tumor_fragments,
    tumor_split_fragments, tumor_discordant_fragments, tumor_ref_reads,
    tumor_ref_pairs, tumor_breakend_af, sv_pon_count, fragile_site,
    line_type, transcript = transcriptId
  )

# Fusions.
be_gene <- breakends |> select(sample_name, id, gene)
end_support <- sv_ends |>
  select(sample_name, vcf_id, tumor_fragments, weak_support, same_bp_n_samples)
fusions <- read_per_sample(samples, "fusion_file")
if (nrow(fusions) > 0) {
  fusions <- fusions |>
    left_join(be_gene, by = c("sample_name", fivePrimeBreakendId = "id")) |>
    rename(gene_5p = gene) |>
    left_join(be_gene, by = c("sample_name", threePrimeBreakendId = "id")) |>
    rename(gene_3p = gene) |>
    mutate(
      chrom_5p = parse_coords(fivePrimeCoords)$chrom,
      chrom_3p = parse_coords(threePrimeCoords)$chrom
    ) |>
    redact_xy(c("chrom_5p", "chrom_3p"), "fusions") |>
    left_join(end_support, by = c("sample_name", fivePrimeVcfId = "vcf_id")) |>
    left_join(
      end_support, by = c("sample_name", threePrimeVcfId = "vcf_id"),
      suffix = c("_5p", "_3p")
    ) |>
    transmute(
      sample_name, fusion = name, gene_5p, gene_3p,
      pipeline_reported = reported == "true", reported_type = reportedType,
      phased, likelihood, breakpoint_5p = fivePrimeCoords,
      breakpoint_3p = threePrimeCoords, fused_exon_up = fusedExonUp,
      fused_exon_down = fusedExonDown, chain_length = chainLength,
      chain_links = chainLinks, chain_terminated = chainTerminated,
      skipped_exons_up = skippedExonsUp, skipped_exons_down = skippedExonsDown,
      domains_kept = domainsKept, domains_lost = domainsLost,
      tumor_fragments_5p, tumor_fragments_3p,
      weak_support = weak_support_5p | weak_support_3p,
      same_bp_n_samples = pmax(same_bp_n_samples_5p, same_bp_n_samples_3p, na.rm = TRUE),
      reportable_reasons = reportableReasons
    )
} else {
  fusions <- tibble(
    sample_name = character(), fusion = character(), gene_5p = character(),
    gene_3p = character(), pipeline_reported = logical()
  )
}

# ---- drivers, viruses, signatures, HLA, QC, DPYD -----------------------

log_msg("drivers and other tables")
drivers <- read_per_sample(samples, "driver_file") |>
  redact_xy("chromosome", "drivers") |>
  transmute(
    sample_name, gene, chromosome, band = chromosomeBand, driver,
    panel_role = category, likelihood_method = likelihoodMethod,
    reported_status = reportedStatus,
    driver_likelihood = as.numeric(driverLikelihood),
    missense, nonsense, splice, inframe, frameshift, biallelic,
    min_copy_number = as.numeric(minCopyNumber),
    max_copy_number = as.numeric(maxCopyNumber)
  )

# The default genes of interest: the genes with a likely somatic driver
# call from the pipeline, those called in the most samples first. The list
# is written next to the report, to be edited and given back with --genes.
default_genes_file <- path(res_dir, "genes_of_interest_default.txt")
if (!genes_given) {
  genes_of_interest <- drivers |>
    filter(
      driver_likelihood >= DEFAULT_GOI_MIN_LIKELIHOOD,
      !str_detect(driver, "^GERMLINE"), gene %in% gene_info$gene
    ) |>
    summarize(n_samples = n_distinct(sample_name), .by = gene) |>
    arrange(desc(n_samples), gene) |>
    pull(gene)
  write_lines(
    c(
      "# Genes of interest of the Prism cohort report: the default list,",
      "# written because the report was run without a gene list (-g).",
      as.character(str_glue(
        "# Rule: a pipeline driver call of likelihood ",
        "{DEFAULT_GOI_MIN_LIKELIHOOD} or higher in at least one sample."
      )),
      "# Order: the number of samples with such a call, highest first.",
      "# To change the list, edit a copy of this file and give it with -g.",
      "# A rerun without -g replaces this file.",
      genes_of_interest
    ),
    default_genes_file
  )
  log_msg(
    "{length(genes_of_interest)} genes of interest (default list): ",
    "{default_genes_file}"
  )
}

viruses <- read_per_sample(samples, "virus_file")
if (nrow(viruses) > 0) {
  viruses <- viruses |>
    transmute(
      sample_name, virus = name, interpretation,
      integrations = as.numeric(integrations),
      percent_covered = as.numeric(percentageCovered),
      mean_coverage = as.numeric(meanCoverage),
      expected_clonal_coverage = as.numeric(expectedClonalCoverage),
      qc_status = qcStatus, pipeline_reported = reported == "true",
      driver_likelihood = driverLikelihood
    )
} else {
  viruses <- tibble(
    sample_name = character(), virus = character(),
    interpretation = character()
  )
}

signatures <- read_per_sample(samples, "sigs_file")
if (nrow(signatures) > 0) {
  signatures <- signatures |>
    transmute(
      sample_name, signature, allocation = as.numeric(allocation),
      percent = round(100 * as.numeric(percent), 2)
    )
}

hla <- read_per_sample(samples, "lilac_file")
if (nrow(hla) > 0) {
  hla <- hla |>
    transmute(
      sample_name, hla_gene = Genes, allele = Allele,
      ref_fragments = as.numeric(RefTotal),
      tumor_fragments = as.numeric(TumorTotal),
      tumor_copy_number = as.numeric(TumorCopyNumber),
      somatic_missense = as.numeric(SomaticMissense),
      somatic_nonsense_or_frameshift = as.numeric(SomaticNonsenseOrFrameshift),
      somatic_splice = as.numeric(SomaticSplice),
      somatic_inframe_indel = as.numeric(SomaticInframeIndel)
    )
}

qc <- samples |>
  mutate(file = map_chr(
    group_dir, \(d) first_file(path(d, "qsee"), "*.qsee.status.tsv.gz")
  )) |>
  read_per_sample("file")
if (nrow(qc) > 0) {
  qc <- qc |>
    filter(str_detect(QcStatus, "WARN|FAIL")) |>
    transmute(
      sample_name, qc_sample = SampleId, qc_sample_type = SampleType,
      tool = SourceTool, metric = FeatureName,
      value = as.numeric(FeatureValue), qc_status = QcStatus,
      flagged_when = FailCondition
    )
}

purple_qc <- samples |>
  filter(file_exists(purple_qc_file)) |>
  select(sample_name, purple_qc_file) |>
  pmap(\(sample_name, purple_qc_file) {
    read_chr(purple_qc_file, col_names = c("key", "value")) |>
      mutate(sample_name = sample_name)
  }) |>
  list_rbind() |>
  filter(key %in% c("QCStatus", "Contamination", "TincLevel", "LohPercent")) |>
  pivot_wider(names_from = key, values_from = value)

dpyd <- read_per_sample(samples, "peach_file")
if (nrow(dpyd) > 0) {
  dpyd <- dpyd |>
    filter(gene == "DPYD") |>
    summarize(
      dpyd_haplotypes = str_c(haplotype, " x", count, collapse = "; "),
      dpyd_function = collapse_unique(`function`, "; "),
      .by = sample_name
    )
} else {
  dpyd <- tibble(
    sample_name = character(), dpyd_haplotypes = character(),
    dpyd_function = character()
  )
}

chord <- read_per_sample(samples, "chord_file")
if (nrow(chord) > 0) {
  chord <- chord |>
    transmute(sample_name, hr_status, p_hrd = as.numeric(p_hrd))
} else {
  chord <- tibble(sample_name = character(), hr_status = character(), p_hrd = numeric())
}

virus_per_sample <- viruses |>
  summarize(
    viruses = collapse_unique(virus, "; "),
    ebv = any(interpretation == "EBV"),
    .by = sample_name
  )

# ---- samples sheet -----------------------------------------------------

samples_sheet <- samples |>
  select(sample_name, sample_class, matched_normal, group, xy_reported, all_of(manifest_cols)) |>
  left_join(
    purity |> transmute(
      sample_name,
      purity = as.numeric(purity), ploidy = as.numeric(ploidy),
      min_purity = as.numeric(minPurity), max_purity = as.numeric(maxPurity),
      min_ploidy = as.numeric(minPloidy), max_ploidy = as.numeric(maxPloidy),
      purple_fit = status,
      whole_genome_duplication = wholeGenomeDuplication == "true",
      ms_status = msStatus, ms_indels_per_mb = as.numeric(msIndelsPerMb),
      tmb_per_mb = as.numeric(tmbPerMb), tmb_status = tmbStatus,
      tml = as.numeric(tml), tml_status = tmlStatus,
      sv_tmb = as.numeric(svTumorMutationalBurden)
    ),
    by = "sample_name"
  ) |>
  left_join(
    purple_qc |> transmute(
      sample_name, purple_qc = QCStatus,
      contamination = as.numeric(Contamination),
      tumor_in_normal = as.numeric(TincLevel),
      loh_fraction = as.numeric(LohPercent)
    ),
    by = "sample_name"
  ) |>
  left_join(chord, by = "sample_name") |>
  left_join(virus_per_sample, by = "sample_name") |>
  left_join(dpyd, by = "sample_name") |>
  left_join(sv_per_sample, by = "sample_name") |>
  left_join(
    variants |> count(sample_name, name = "coding_variants"),
    by = "sample_name"
  ) |>
  mutate(
    ebv = replace_na(ebv, FALSE),
    dpyd_haplotypes = replace_na(dpyd_haplotypes, "not assessed"),
    dpyd_function = replace_na(dpyd_function, "not assessed"),
    sv_weak_support_pct = pct(sv_weak_support, sv_total),
    coding_variants = replace_na(coding_variants, 0L)
  ) |>
  relocate(sv_weak_support_pct, .after = sv_weak_support)

# ---- cohort summary ----------------------------------------------------

summarize_group <- function(df, group_label) {
  flag <- function(metric, x) {
    tibble(
      metric = metric, n = sum(x, na.rm = TRUE), N = sum(!is.na(x)),
      pct = pct(n, N)
    )
  }
  num <- function(metric, x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(tibble(metric = metric, N = 0L))
    }
    tibble(
      metric = metric, N = length(x), median = signif(median(x), 6), min = min(x),
      max = max(x)
    )
  }
  assessed <- df$dpyd_function != "not assessed"
  bind_rows(
    tibble(metric = "samples", n = nrow(df), N = nrow(samples_sheet), pct = pct(n, N)),
    flag("matched normal", df$matched_normal),
    flag("whole genome duplication", df$whole_genome_duplication),
    flag("MSI", df$ms_status == "MSI"),
    flag("TMB high", df$tmb_status == "HIGH"),
    flag("TML high", df$tml_status == "HIGH"),
    flag("HR deficient (CHORD; N = samples with a CHORD result)",
         if_else(is.na(df$hr_status), NA, df$hr_status == "HR_deficient")),
    flag("EBV detected", df$ebv),
    flag("any virus detected", !is.na(df$viruses)),
    flag("DPYD assessed (needs a matched normal)", assessed),
    flag("DPYD normal function (N = assessed)",
         if_else(assessed, df$dpyd_function == "Normal Function", NA)),
    flag("PURPLE QC PASS", df$purple_qc == "PASS"),
    num("purity", df$purity),
    num("ploidy", df$ploidy),
    num("TMB per Mb", df$tmb_per_mb),
    num("MS indels per Mb", df$ms_indels_per_mb),
    num("SV count (PURPLE)", df$sv_tmb),
    num("SV count (LINX)", df$sv_total),
    num("SVs with weak support, %", df$sv_weak_support_pct),
    num("coding variants", df$coding_variants)
  ) |>
    mutate(group = group_label, .after = metric)
}

summary_full <- bind_rows(
  summarize_group(samples_sheet, "all samples"),
  map(classes, \(cl) {
    summarize_group(
      samples_sheet |> filter(sample_class == cl),
      str_glue("sample_class: {cl}")
    )
  }),
  summarize_group(samples_sheet |> filter(matched_normal), "matched normal: yes"),
  summarize_group(samples_sheet |> filter(!matched_normal), "matched normal: no")
) |>
  mutate(metric = fct_inorder(metric)) |>
  arrange(metric) |>
  mutate(metric = as.character(metric))

# Wide form for reading: one row per yes/no metric, one column per
# sample class, each cell "n of N (pct%)". summary_full keeps the
# numeric metrics and the matched-normal groups.
cohort_summary <- summary_full |>
  filter(!is.na(n), !str_detect(group, "^matched normal:")) |>
  mutate(
    group = str_remove(group, "^sample_class: "),
    value = if_else(N == 0, "0 of 0", as.character(str_glue("{n} of {N} ({pct}%)")))
  ) |>
  select(metric, group, value) |>
  pivot_wider(names_from = group, values_from = value)

# ---- chromosome arms ---------------------------------------------------

arms <- read_per_sample(samples, "arm_file")
if (nrow(arms) > 0) {
  arms <- arms |>
    redact_xy("chromosome", "arm_cn") |>
    left_join(ploidy_tbl, by = "sample_name") |>
    mutate(
      arm = str_c(chromosome, str_to_lower(arm)),
      value = round(as.numeric(medianCopyNumber) - ploidy, 2)
    ) |>
    mutate(arm = fct_inorder(arm)) |>
    select(arm, sample_name, value) |>
    pivot_wider(names_from = sample_name, values_from = value) |>
    arrange(arm) |>
    mutate(arm = as.character(arm))
}

# ---- genes of interest -------------------------------------------------

# Every event in a gene of interest, one row per event. Each kind of event
# is written into the same columns: what it is (event), where (location),
# the evidence (support) and the reasons for caution (caution).
GOI_EVENT_TYPES <- c(
  "mutation", "amp", "partial_amp", "deep_del", "loh", "disruption", "fusion"
)

# Join text parts with "; ", leaving out the missing ones.
join_parts <- function(...) {
  pmap_chr(list(...), \(...) {
    parts <- discard(c(...), is.na)
    if (length(parts) == 0) NA_character_ else str_c(parts, collapse = "; ")
  })
}

# Numbers as text: rounded, and whole numbers without an exponent.
num_chr <- function(x, digits = 1) as.character(round(x, digits))
int_chr <- function(x) as.character(as.integer(x))

goi_events <- NULL
if (length(genes_of_interest) > 0) {
  goi_mutations <- variants |>
    filter(gene %in% genes_of_interest) |>
    transmute(
      sample_name, gene, event_type = "mutation",
      event = str_c(
        if_else(is.na(hgvs_protein) | hgvs_protein == "", hgvs_coding, hgvs_protein),
        " (", mutation_class, if_else(hotspot, ", hotspot", ""),
        if_else(biallelic, ", biallelic", ""), ")"
      ),
      location = str_c(
        chromosome, ":", int_chr(position), " ", str_trunc(ref, 20), ">",
        str_trunc(alt, 20)
      ),
      support = join_parts(
        str_c("VAF ", num_chr(tumor_vaf, 2)),
        str_c(int_chr(tumor_alt_reads), " of ", int_chr(tumor_depth), " reads"),
        str_c("copy number ", num_chr(copy_number))
      ),
      caution = NA_character_, pipeline_reported = reported, pos = position
    )

  # Copy number events. The location is the span of the gene. A partial
  # amplification (the highest copy number of the gene passes the
  # amplification threshold, the lowest does not) is listed here only;
  # gene_cn has no row for it.
  goi_cn <- gene_cn |>
    filter(gene %in% genes_of_interest) |>
    mutate(
      location = str_c(chromosome, ":", int_chr(start), "-", int_chr(end)),
      ploidy_txt = str_c("ploidy ", num_chr(sample_ploidy, 2)),
      segment_txt = join_parts(
        str_c("segment ends ", min_region_start_support, "/", min_region_end_support),
        str_c(depth_window_count, " depth windows")
      ),
      pos = start
    )
  goi_cn_events <- bind_rows(
    goi_cn |>
      filter(cn_event == "AMP") |>
      mutate(
        event_type = "amp",
        event = str_c(
          "min copy number ", num_chr(min_copy_number), " (",
          num_chr(min_copy_number / sample_ploidy), " x ploidy)"
        ),
        support = join_parts(
          str_c("max copy number ", num_chr(max_copy_number)), ploidy_txt,
          segment_txt
        )
      ),
    goi_cn |>
      filter(
        coalesce(cn_event != "AMP", TRUE),
        max_copy_number > AMP_PLOIDY_RATIO * sample_ploidy
      ) |>
      mutate(
        event_type = "partial_amp",
        event = str_c(
          "max copy number ", num_chr(max_copy_number), " (",
          num_chr(max_copy_number / sample_ploidy), " x ploidy)"
        ),
        support = join_parts(
          str_c("min copy number ", num_chr(min_copy_number)), ploidy_txt
        )
      ),
    goi_cn |>
      filter(cn_event == "DEEP_DEL") |>
      mutate(
        event_type = "deep_del",
        event = str_c("min copy number ", num_chr(min_copy_number)),
        support = join_parts(
          str_c("max copy number ", num_chr(max_copy_number)), ploidy_txt,
          segment_txt
        )
      ),
    goi_cn |>
      filter(loh) |>
      mutate(
        event_type = "loh",
        event = str_c(
          "minor allele copy number ", num_chr(min_minor_allele_copy_number)
        ),
        support = join_parts(
          str_c(
            "copy number ", num_chr(min_copy_number), " to ",
            num_chr(max_copy_number)
          ),
          ploidy_txt
        )
      )
  ) |>
    select(sample_name, gene, event_type, event, location, support, pos)

  goi_disruptions <- disruptions |>
    filter(gene %in% genes_of_interest) |>
    transmute(
      sample_name, gene, event_type = "disruption",
      event = str_c(
        sv_type, " breakend ",
        case_when(
          region == "INTRONIC" ~ str_c("between exons ", exon_up, " and ", exon_down),
          region == "EXONIC" ~ str_c("in exon ", exon_up),
          .default = str_c("(", str_to_lower(region), ")")
        ),
        if_else(disruptive, "", ", not disruptive")
      ),
      location = join_parts(breakpoint, str_c("other end ", na_if(other_end, ""))),
      support = join_parts(
        str_c(int_chr(tumor_fragments), " fragments"),
        str_c("AF ", num_chr(tumor_breakend_af, 2)),
        str_c(
          "junction copy number ", num_chr(junction_cn_min), " to ",
          num_chr(junction_cn_max)
        ),
        str_c("undisrupted copy number ", num_chr(undisrupted_copy_number))
      ),
      caution = join_parts(
        if_else(weak_support, "weak support", NA_character_),
        if_else(
          same_bp_flag,
          str_c("shared breakpoint (", same_bp_n_samples, " samples)"),
          NA_character_
        )
      ),
      pipeline_reported, pos = parse_coords(breakpoint)$pos
    )

  # A fusion gives one row for each of its genes that is a gene of interest.
  goi_fusions <- NULL
  if (nrow(fusions) > 0) {
    goi_fusions <- bind_rows(
      fusions |> mutate(
        gene = gene_5p, side = "5' gene", own_bp = breakpoint_5p,
        other_bp = breakpoint_3p
      ),
      fusions |> mutate(
        gene = gene_3p, side = "3' gene", own_bp = breakpoint_3p,
        other_bp = breakpoint_5p
      )
    ) |>
      filter(gene %in% genes_of_interest) |>
      transmute(
        sample_name, gene, event_type = "fusion",
        event = coalesce(
          str_c(
            fusion, " (", side, "; exon ", fused_exon_up, " to exon ",
            fused_exon_down, "; ",
            str_to_lower(str_replace_all(phased, "_", " ")), ")"
          ),
          fusion
        ),
        location = join_parts(own_bp, str_c("other end ", na_if(other_bp, ""))),
        support = join_parts(
          str_c(
            int_chr(tumor_fragments_5p), " and ", int_chr(tumor_fragments_3p),
            " fragments (5' and 3' end)"
          ),
          if_else(
            chain_links != "0", str_c("chain of ", chain_links, " links"),
            NA_character_
          )
        ),
        caution = join_parts(
          if_else(weak_support, "weak support", NA_character_),
          if_else(
            same_bp_n_samples > 1,
            str_c("shared breakpoint (", same_bp_n_samples, " samples)"),
            NA_character_
          )
        ),
        pipeline_reported, pos = parse_coords(own_bp)$pos
      )
  }

  # The pipeline's driver call for the same kind of event in the gene and
  # sample. Fusions are not in the driver list.
  driver_key <- drivers |>
    mutate(
      event_type = case_when(
        driver == "MUTATION" ~ "mutation",
        driver == "AMP" ~ "amp",
        driver == "PARTIAL_AMP" ~ "partial_amp",
        driver == "DEL" ~ "deep_del",
        driver == "LOH" ~ "loh",
        str_detect(driver, "DISRUPTION$") ~ "disruption"
      )
    ) |>
    filter(!is.na(event_type)) |>
    summarize(
      driver_likelihood = max(driver_likelihood),
      driver_reported = any(reported_status == "REPORTED"),
      .by = c(sample_name, gene, event_type)
    )

  # Copy number events have no reported flag of their own; for them
  # pipeline_reported says whether the driver list has the event.
  goi_events <- bind_rows(
    goi_mutations, goi_cn_events, goi_disruptions, goi_fusions
  ) |>
    left_join(driver_key, by = c("sample_name", "gene", "event_type")) |>
    mutate(
      pipeline_reported = coalesce(pipeline_reported, driver_reported, FALSE),
      gene = factor(gene, levels = genes_of_interest),
      event_type = factor(event_type, levels = GOI_EVENT_TYPES)
    ) |>
    arrange(gene, event_type, sample_name, pos) |>
    mutate(across(c(gene, event_type), as.character)) |>
    select(
      gene, sample_name, event_type, event, location, support, caution,
      pipeline_reported, driver_likelihood
    )
  log_msg(
    "genes of interest: {length(genes_of_interest)} genes, ",
    "{nrow(goi_events)} events"
  )
}

# ---- manifest and glossary ----------------------------------------------

# Text with values from this run filled in.
g <- function(...) as.character(str_glue(..., .envir = parent.frame()))

review_file <- "prism_sv_same_breakpoint.xlsx"
redacted_txt <- if (length(redacted) == 0) {
  "none"
} else {
  str_c(names(redacted), " ", unlist(redacted), collapse = "; ")
}

# Where the genes of interest come from.
goi_txt <- if (genes_given) {
  g("the {length(genes_of_interest)} genes of the gene list given to the report")
} else {
  g(
    "the default list, used because no gene list was given: the ",
    "{length(genes_of_interest)} genes with a pipeline driver call of ",
    "likelihood {DEFAULT_GOI_MIN_LIKELIHOOD} or higher in at least one ",
    "sample (the list is in {path_file(default_genes_file)}, next to this ",
    "workbook)"
  )
}

# manifest: what each sheet holds. One entry per sheet the report can have.
manifest_text <- tribble(
  ~sheet, ~one_row_is, ~contents,
  "manifest", "a sheet",
  g(
    "This sheet: what each sheet holds. Report made by Prism ({opts$version}) ",
    "on {Sys.Date()} from nf-core/oncoanalyser results of {nrow(samples)} ",
    "tumor samples, {sum(samples$matched_normal)} run with a matched normal ",
    "and {sum(!samples$matched_normal)} without. The report counts and lists ",
    "events; it applies no significance test and leaves no sample out. Events ",
    "on chrX and chrY are left out for {nrow(samples) - length(xy_ok_samples)} ",
    "samples (see xy_reported in the glossary); rows left out: {redacted_txt}. ",
    "The same tables are in data/ as CSV files for use in code; in Excel open ",
    "this workbook, because Excel changes some gene names to dates when it ",
    "opens a CSV."
  ),
  "genes_of_interest_events", "an event in a gene of interest, in a sample",
  g(
    "Every event in the genes of interest, one row per event, sorted by gene ",
    "and kind of event: each coding mutation, amplification, partial ",
    "amplification, deep deletion, LOH, disruptive breakend and predicted ",
    "fusion, with what it is, where it is, the evidence for it and whether ",
    "the pipeline reports it. The genes of interest are {goi_txt}. A gene ",
    "and sample with no event has no row. Partial amplifications are listed ",
    "in this sheet only; gene_cn has no row for them."
  ),
  "samples", "a tumor sample",
  "The sample's class and whether it has a matched normal, then its genome-wide results: purity and ploidy, whole genome duplication, microsatellite status, mutational burden, homologous recombination status, viruses, DPYD status, quality checks, and its number of structural variants and coding variants.",
  "cohort_summary", "a yes/no measure",
  "How many samples have each sample-level feature (whole genome duplication, MSI, high mutational burden, EBV, DPYD status and others), for all samples and for each sample class. Each cell reads 'n of N (pct%)'.",
  "somatic_variants", "a somatic variant in a sample",
  "Every somatic point mutation and small indel that passes the pipeline's filters and changes the protein of the gene's main transcript (missense, inframe, nonsense, frameshift, splice), with its position, protein change, read support and copy number context. In a sample without a matched normal, some of these can be germline variants that the filters did not remove.",
  "disruptions", "a structural variant breakend inside a gene, in a sample",
  "Every breakend that breaks a gene's main transcript, with the genomic address of the breakpoint, the other end of the structural variant, its position in the gene, and the evidence for it (read support, copy number of the junction, quality). Sorted by gene and position, so all breakends of one gene are together. weak_support and same_bp_flag mark the breakends to treat with caution.",
  "gene_cn", "a gene with a copy number event, in a sample",
  g(
    "Gene copy number events: amplification (minimum copy number of the gene ",
    "above {AMP_PLOIDY_RATIO} x the sample's ploidy) and deep deletion ",
    "(minimum copy number below {DEEP_DEL_MAX_CN}) for every gene, and LOH ",
    "(minor allele copy number below {LOH_MAX_MINOR_CN}, not a deep deletion) ",
    "for the genes of the driver panel and the genes of interest. One large ",
    "event gives one row for every gene it covers. data/gene_cn.csv also has ",
    "the LOH rows of all other genes."
  ),
  "fusions", "a gene fusion in a sample",
  "Every gene fusion the pipeline predicts from the structural variants, whether or not the pipeline reports it, with both breakpoints, the exons joined, the reading frame and the read support. pipeline_reported marks the fusions the pipeline itself reports.",
  "viruses", "a virus found in a sample",
  "Viral sequence found in the tumor reads: the virus, how much of its genome is covered and how deeply, and the number of integration sites in the human genome.",
  "drivers", "a driver event in a sample",
  "The pipeline's own list of driver events (mutation, amplification, deletion, disruption) in the genes of its driver panel, each with the pipeline's estimate of the likelihood that the event is a driver.",
  "arm_cn", "a chromosome arm",
  "Copy number of each chromosome arm in each sample, relative to the sample's ploidy: the median copy number of the arm minus the ploidy. Positive values are gains, negative values are losses. One column per sample.",
  "qc", "a quality warning for a sample",
  "The quality measures that the pipeline's QC tool flags as WARN or FAIL, for the tumor or its normal, with the value and the condition that set the flag. Measures that pass are not listed.",
  "signatures", "a mutational signature in a sample",
  "The number and percent of each sample's somatic variants assigned to each mutational signature. Only samples with a matched normal have this result.",
  "hla", "an HLA allele in a sample",
  "The HLA class I and class II alleles of each sample, with read support, the copy number of the allele in the tumor and the number of somatic variants in it.",
  "summary_full", "a measure in a group of samples",
  "The counts of cohort_summary as numbers, also for the samples with and without a matched normal, and the median, minimum and maximum of the numeric measures (purity, ploidy, mutational burden, number of structural variants, number of coding variants).",
  "glossary", "a column",
  "The meaning of every column in every sheet. A column that is in more than one sheet is explained once.",
  g("{review_file} (separate file)"), "a shared breakpoint in a sample",
  "Not in this workbook. Every structural variant breakpoint that is found at the same base pair in more than one sample, inside or outside genes. Independent events are not expected to share a breakpoint, so these are for review."
)

# glossary: every column, once. Columns named after a sample or a sample
# class are listed under a placeholder name.
glossary_text <- tribble(
  ~column, ~description,
  # sample
  "sample_name", "Name of the tumor sample.",
  "sample_class", "Class of the sample, from the sample manifest given to the report (for example Tumor or CellLine). Tumor for every sample when no manifest is given.",
  "matched_normal", "TRUE if the sample was run with a matched normal sample. Without one, somatic calls can include germline variants that the filters did not remove, and DPYD status, HR status and signatures are not available.",
  "group", "Group name of the sample in the pipeline run (the name of its result folder).",
  "xy_reported", "TRUE if events on chrX and chrY are reported for this sample. That needs both a sex in the sample manifest and a matched normal; for all other samples the X/Y events are left out of every sheet.",
  "sex", "Sex of the patient, from the sample manifest.",
  "purity", "Estimated fraction of tumor cells in the sample (PURPLE).",
  "ploidy", "Estimated average copy number of the tumor genome (PURPLE).",
  "min_purity", "Lowest purity among the fits PURPLE scores as close to the best fit.",
  "max_purity", "Highest purity among the fits PURPLE scores as close to the best fit.",
  "min_ploidy", "Lowest ploidy among the fits PURPLE scores as close to the best fit.",
  "max_ploidy", "Highest ploidy among the fits PURPLE scores as close to the best fit.",
  "purple_fit", "How PURPLE arrived at the purity and ploidy: NORMAL (from copy number), HIGHLY_DIPLOID, SOMATIC (from somatic variant allele frequencies) or NO_TUMOR.",
  "whole_genome_duplication", "TRUE if PURPLE calls a whole genome duplication.",
  "ms_status", "Microsatellite status from PURPLE: MSI (unstable) or MSS (stable).",
  "ms_indels_per_mb", "Number of indels in microsatellite regions per megabase; the measure behind ms_status.",
  "tmb_per_mb", "Tumor mutational burden: somatic variants per megabase (PURPLE).",
  "tmb_status", "HIGH or LOW tumor mutational burden, as called by PURPLE.",
  "tml", "Tumor mutational load: the number of somatic missense variants in the sample (PURPLE).",
  "tml_status", "HIGH or LOW tumor mutational load, as called by PURPLE.",
  "sv_tmb", "Number of structural variants in the sample, as counted by PURPLE.",
  "purple_qc", "PURPLE's quality status for the sample: PASS, or the warnings and failures raised.",
  "contamination", "Estimated fraction of reads from another person in the tumor sample (0 if none found).",
  "tumor_in_normal", "Estimated level of tumor DNA in the matched normal sample.",
  "loh_fraction", "Fraction of the genome with loss of heterozygosity (PURPLE).",
  "hr_status", "Homologous recombination status from CHORD: HR_proficient, HR_deficient or cannot_be_determined. Empty without a matched normal.",
  "p_hrd", "CHORD's probability that the sample is homologous recombination deficient (0 to 1).",
  "viruses", "The viruses found in the sample; details in the viruses sheet.",
  "ebv", "TRUE if Epstein-Barr virus was found in the sample.",
  "dpyd_haplotypes", "DPYD haplotypes of the patient with the number of copies of each (PEACH). 'not assessed' without a matched normal, because the genotype is read from the normal sample.",
  "dpyd_function", "Function of the DPYD haplotypes (for example Normal Function). 'not assessed' without a matched normal.",
  "sv_total", "Number of structural variants in the sample, as listed by LINX.",
  "sv_weak_support", "Number of the sample's structural variants with weak support (see weak_support).",
  "sv_weak_support_pct", "sv_weak_support as a percent of sv_total.",
  "coding_variants", "Number of the sample's rows in somatic_variants.",
  "<manifest column>", "Any other column of the sample manifest, copied as given.",
  # genes of interest: events
  "event_type", g(
    "Kind of event: mutation (a row of somatic_variants), amp ",
    "(amplification), partial_amp (partial amplification), deep_del (deep ",
    "deletion), loh (loss of heterozygosity), disruption (a row of ",
    "disruptions) or fusion (a row of fusions). amp, deep_del and loh are ",
    "defined under cn_event and loh. partial_amp: the highest copy number of ",
    "the gene is above {AMP_PLOIDY_RATIO} x the sample's ploidy but the ",
    "lowest is not."
  ),
  "event", "What the event is. Mutation: the protein change (the coding change if there is none), the mutation class, and 'hotspot' or 'biallelic' where they apply. Copy number events: the copy number that makes the call. Disruption: the type of structural variant and where the breakend falls in the transcript. Fusion: its name, which end this gene is, the exons joined and the reading frame.",
  "location", "Where the event is (GRCh38). Mutation: position and alleles. Copy number events: the span of the gene. Disruption and fusion: the breakpoint in this gene as chromosome:position:orientation (see breakpoint), then the other end of the structural variant or fusion.",
  "support", "The evidence for the event, as text. Mutation: allele fraction (VAF), reads with the variant out of all reads, copy number. Copy number events: the gene's other copy number value, the sample's ploidy, what marks the ends of the segment with the lowest copy number and its number of read depth windows. Disruption: supporting fragments, allele fraction (AF), copy number of the junction, copies of the gene left intact. Fusion: supporting fragments at each end.",
  "caution", "Reasons to treat a breakend or fusion with caution: 'weak support' (see weak_support) and 'shared breakpoint' with the number of samples that have it (see same_bp_flag). Empty if neither applies, and for mutations and copy number events.",
  # summaries
  "metric", "The measure the row is about. In qc: the name of the quality measure.",
  "all samples", "The count over all samples, as 'n of N (pct%)'.",
  "<sample_class>", "One column per sample class: the count within that class, as 'n of N (pct%)'.",
  "sample_group", "The set of samples the row counts: all samples, one sample class, or the samples with or without a matched normal.",
  "n", "Number of samples with the feature.",
  "N", "Number of samples the feature could be counted in.",
  "pct", "100 x n / N.",
  "median", "Median of a numeric measure over the samples of the group.",
  "min", "Lowest value of a numeric measure over the samples of the group.",
  "max", "Highest value of a numeric measure over the samples of the group.",
  # genes
  "gene", "Gene symbol.",
  "panel_role", "Role the driver panel gives the gene: ONCO (oncogene) or TSG (tumor suppressor). The pipeline scores a gene with the model for its role, so an activating mutation in a gene classed as TSG gets a low driver likelihood.",
  "chromosome", "Chromosome of the gene or variant.",
  "band", "Chromosome band of the gene.",
  "gene_length_kb", "Length of the gene in kilobases. Given because a longer gene can collect more breakends.",
  # somatic variants
  "position", "Position of the variant on the chromosome (GRCh38).",
  "ref", "Reference allele.",
  "alt", "Variant allele.",
  "mutation_class", "missense, inframe (indel that keeps the reading frame), truncating (nonsense or frameshift) or splice.",
  "canonical_effect", "Effect of the variant on the gene's main (canonical) transcript, in sequence ontology terms.",
  "hgvs_coding", "The change in the coding sequence, in HGVS notation.",
  "hgvs_protein", "The change in the protein, in HGVS notation.",
  "transcript", "Ensembl ID of the transcript the effect or breakend is described on (the gene's canonical transcript).",
  "hotspot", "TRUE if the variant is at a known cancer hotspot.",
  "biallelic", "In somatic_variants: TRUE if PURPLE finds no remaining wild-type copy at the variant. In drivers: the same flag for the driver event.",
  "reported", "TRUE if the pipeline marks the variant as reportable (a coding variant in a driver panel gene).",
  "tier", "Confidence tier of the variant caller (SAGE): HOTSPOT, PANEL, HIGH_CONFIDENCE or LOW_CONFIDENCE.",
  "qual", "Quality score of the variant call.",
  "tumor_depth", "Number of tumor reads at the variant position.",
  "tumor_alt_reads", "Number of tumor reads with the variant allele.",
  "tumor_vaf", "Fraction of tumor reads with the variant allele.",
  "purity_adjusted_vaf", "Variant allele fraction in the tumor cells, after correcting for purity and copy number (PURPLE).",
  "copy_number", "Copy number of the tumor genome at the variant position (PURPLE).",
  "variant_copy_number", "Number of copies carrying the variant (PURPLE).",
  "minor_allele_copy_number", "Copy number of the less frequent parental allele at the variant position (PURPLE).",
  "subclonal_likelihood", "Probability that the variant is in only part of the tumor cells (0 to 1).",
  "mappability", "How uniquely reads map at the position (1 = unique).",
  "gnomad_freq", "Frequency of the variant in the gnomAD population database. Empty if not found there.",
  "pon_count", "Number of samples in the pipeline's panel of normals with the variant. Empty if none.",
  # structural variants
  "breakpoint", "Genomic address of the breakend: chromosome:position:orientation (GRCh38). Orientation 1: the retained sequence is to the left of the position; -1: to the right.",
  "other_end", "Address of the other end of the same structural variant. Empty for a single breakend, where the other end could not be placed.",
  "sv_type", "Type of the structural variant: DEL (deletion), DUP (duplication), INV (inversion), BND (translocation), INS (insertion), SGL (single breakend, other end not placed), INF (inferred from a copy number change).",
  "disruptive", "TRUE if the breakend breaks the transcript, so that it cannot be made intact from this copy (LINX).",
  "pipeline_reported", "TRUE if the pipeline itself reports this event. For the copy number events of genes_of_interest_events: TRUE if the pipeline's driver list (drivers sheet) has the event.",
  "weak_support", "TRUE if LINX classes the structural variant as a likely artifact, in most cases because the junction is present at under about half a copy with no matching copy number change. In fusions: TRUE if either end has weak support.",
  "same_bp_flag", "TRUE if the same breakpoint (chromosome, position and orientation, to the base pair) is found in more than one sample. Independent events are not expected to share a breakpoint. Breakends inferred from copy number have no exact position and are not compared.",
  "same_bp_n_samples", "Number of samples with this same breakpoint (1 = this sample only). In fusions: the higher of the two ends.",
  "same_bp_samples", "The samples with this same breakpoint.",
  "nearest_other_sample_bp", "Distance in base pairs to the nearest breakend in the same gene in another sample. 0 means the same position. Empty if no other sample has a breakend in the gene.",
  "region", "Where in the transcript the breakend falls: INTRONIC, EXONIC or UPSTREAM.",
  "coding", "Coding context at the breakend: CODING, UTR_5P, UTR_3P, NON_CODING or ENHANCER.",
  "exon_up", "Number of the nearest exon upstream of the breakend in the transcript.",
  "exon_down", "Number of the nearest exon downstream of the breakend in the transcript.",
  "gene_orientation", "Which part of the gene stays attached at the breakend: UPSTREAM (the 5' part) or DOWNSTREAM (the 3' part).",
  "undisrupted_copy_number", "Estimated number of copies of the gene that this breakend leaves intact.",
  "junction_cn_min", "Lower estimate of the number of copies of the tumor genome that carry the junction (LINX, purity adjusted).",
  "junction_cn_max", "Upper estimate of the number of copies that carry the junction.",
  "cluster_id", "ID of the LINX cluster the structural variant belongs to, within the sample. LINX groups variants that it judges to come from one event.",
  "cluster_category", "LINX's category for the cluster: SIMPLE, COMPLEX, INCOMPLETE, ARTIFACT and others. ARTIFACT sets weak_support.",
  "cluster_type", "LINX's resolved type of the cluster (for example DEL, RECIP_TRANS, COMPLEX, LOW_VAF, SGL_PAIR_INS).",
  "cluster_sv_count", "Number of structural variants in the cluster.",
  "sv_qual", "Quality score of the structural variant call.",
  "sv_filter", "Filter status of the structural variant call (PASS if it passed all filters).",
  "tumor_fragments", "In disruptions: number of tumor read fragments supporting the breakend. In hla: number of tumor read fragments assigned to the allele.",
  "tumor_split_fragments", "Supporting fragments with a read that crosses the breakpoint.",
  "tumor_discordant_fragments", "Supporting fragments whose two reads map to either side of the junction.",
  "tumor_ref_reads", "Number of tumor reads that cross the breakpoint position and match the reference.",
  "tumor_ref_pairs", "Number of tumor read pairs that span the breakpoint position and match the reference.",
  "tumor_breakend_af", "Fraction of the tumor fragments at the breakpoint that support the breakend.",
  "sv_pon_count", "Number of samples in the pipeline's panel of normals with this structural variant. Empty if none.",
  "fragile_site", "TRUE if the breakend is in a known fragile site of the genome.",
  "line_type", "Whether the breakend is in a LINE element that is a known or suspected source of insertions: NONE, KNOWN or SUSPECT.",
  "sv_genes_start", "Genes at the first end of the structural variant.",
  "sv_genes_end", "Genes at the second end of the structural variant.",
  # copy number
  "start", "Start position of the gene (GRCh38).",
  "end", "End position of the gene.",
  "cn_event", g(
    "AMP: the minimum copy number of the gene is above {AMP_PLOIDY_RATIO} x ",
    "the sample's ploidy. DEEP_DEL: it is below {DEEP_DEL_MAX_CN}. Empty: ",
    "neither. These are the rules the pipeline uses for its driver genes, ",
    "applied here to every gene."
  ),
  "loh", g(
    "TRUE if the gene has loss of heterozygosity in the sample: minor allele ",
    "copy number below {LOH_MAX_MINOR_CN}, and not a deep deletion."
  ),
  "min_copy_number", "Lowest copy number over the gene in the sample.",
  "max_copy_number", "Highest copy number over the gene in the sample.",
  "min_minor_allele_copy_number", "Lowest copy number of the less frequent parental allele over the gene.",
  "sample_ploidy", "Ploidy of the sample (as in samples).",
  "min_copy_number_over_ploidy", "min_copy_number divided by sample_ploidy.",
  "somatic_regions", "Number of copy number segments that overlap the gene.",
  "min_region_start_support", "What marks the start of the segment with the lowest copy number: a structural variant type (DEL, DUP, INV, BND, SGL, INF), CENTROMERE, TELOMERE or NONE.",
  "min_region_end_support", "What marks the end of the segment with the lowest copy number.",
  "min_region_method", "How PURPLE set the copy number of that segment: BAF_WEIGHTED (from read depth and allele ratios), STRUCTURAL_VARIANT (from the structural variants at its ends), and others.",
  "depth_window_count", "Number of read depth windows in that segment. 0 means its copy number does not rest on read depth.",
  "arm", "Chromosome arm (for example 1p).",
  "<sample name>", "One column per sample: the median copy number of the arm minus the sample's ploidy.",
  # fusions
  "fusion", "Name of the fusion: 5' gene, then 3' gene.",
  "gene_5p", "The gene that gives the start (5' part) of the fusion.",
  "gene_3p", "The gene that gives the end (3' part) of the fusion.",
  "reported_type", "Category of the fusion in the pipeline's fusion knowledge base (for example KNOWN_PAIR, PROMISCUOUS_5, PROMISCUOUS_3); NONE if it is in none.",
  "phased", "Reading frame of the fusion: INFRAME, SKIPPED_EXONS (in frame after skipping exons) or OUT_OF_FRAME.",
  "likelihood", "The pipeline's likelihood that the fusion is a driver: HIGH, LOW or NA.",
  "breakpoint_5p", "Breakpoint in the 5' gene (chromosome:position:orientation).",
  "breakpoint_3p", "Breakpoint in the 3' gene.",
  "fused_exon_up", "Last exon of the 5' gene in the fusion.",
  "fused_exon_down", "First exon of the 3' gene in the fusion.",
  "chain_length", "Number of structural variants chained together to join the two genes (0 = one variant joins them directly).",
  "chain_links", "Number of links in that chain.",
  "chain_terminated", "TRUE if another breakend in the chain interrupts the fusion transcript.",
  "skipped_exons_up", "Number of exons of the 5' gene that must be skipped for the fusion to be in frame.",
  "skipped_exons_down", "Number of exons of the 3' gene that must be skipped.",
  "domains_kept", "Protein domains of the 3' gene kept in the fusion.",
  "domains_lost", "Protein domains of the 3' gene lost in the fusion.",
  "tumor_fragments_5p", "Number of tumor read fragments supporting the breakend in the 5' gene.",
  "tumor_fragments_3p", "Number of tumor read fragments supporting the breakend in the 3' gene.",
  "reportable_reasons", "Why the pipeline does not report the fusion (OK if it does).",
  # viruses
  "virus", "Name of the virus reference the reads match.",
  "interpretation", "The pipeline's label for viruses it tracks (for example EBV); null for others.",
  "integrations", "Number of sites where viral sequence is joined to the human genome.",
  "percent_covered", "Percent of the virus genome covered by reads.",
  "mean_coverage", "Mean read depth over the virus genome.",
  "expected_clonal_coverage", "Read depth expected for one copy of the virus in every tumor cell, given the sample's purity and depth.",
  "qc_status", "In viruses: quality status of the virus call. In qc: WARN or FAIL.",
  "driver_likelihood", "In drivers: the pipeline's estimate that the event is a driver (0 to 1). In genes_of_interest_events: that estimate for the same kind of event in the gene and sample, empty if the driver list has none; mutations are scored per gene, so all mutations of a gene in a sample show the same value. In viruses: HIGH, LOW or UNKNOWN.",
  # drivers
  "driver", "Kind of driver event: MUTATION, AMP, PARTIAL_AMP, DEL, DISRUPTION, LOH and others.",
  "likelihood_method", "How the driver likelihood was set: DNDS (mutation rate model), HOTSPOT, BIALLELIC, AMP, DEL or DISRUPTION.",
  "reported_status", "Whether the pipeline reports the driver event.",
  "missense", "Number of missense variants in the gene behind the driver event.",
  "nonsense", "Number of nonsense variants in the gene behind the driver event.",
  "splice", "Number of splice variants in the gene behind the driver event.",
  "inframe", "Number of inframe indels in the gene behind the driver event.",
  "frameshift", "Number of frameshift variants in the gene behind the driver event.",
  # qc
  "qc_sample", "The sample the quality measure was taken on: the tumor or its matched normal.",
  "qc_sample_type", "TUMOR or NORMAL.",
  "tool", "The pipeline tool the measure comes from.",
  "value", "Value of the quality measure.",
  "flagged_when", "The condition on the value that sets the flag.",
  # signatures
  "signature", "Name of the mutational signature (Sig1 to Sig30 are the COSMIC version 2 signatures; MISALLOC is the part that fits none).",
  "allocation", "Number of the sample's somatic variants assigned to the signature.",
  "percent", "allocation as a percent of the sample's variants.",
  # hla
  "hla_gene", "MHC_CLASS_1 for the HLA-A, -B and -C alleles; the gene name for class II alleles.",
  "allele", "The HLA allele.",
  "ref_fragments", "Number of read fragments of the reference sample assigned to the allele.",
  "tumor_copy_number", "Copy number of the allele in the tumor.",
  "somatic_missense", "Number of somatic missense variants in the allele.",
  "somatic_nonsense_or_frameshift", "Number of somatic nonsense or frameshift variants in the allele.",
  "somatic_splice", "Number of somatic splice variants in the allele.",
  "somatic_inframe_indel", "Number of somatic inframe indels in the allele.",
  # manifest and glossary
  "sheet", "Name of the sheet.",
  "one_row_is", "What one row of the sheet stands for.",
  "rows", "Number of rows in the sheet, without the header.",
  "columns", "Number of columns in the sheet.",
  "contents", "What the sheet holds.",
  "column", "Name of the column. Names in angle brackets stand for a set of columns, one per sample or per sample class.",
  "description", "What the column means.",
  "sheets", "The sheets that have the column."
)

# The glossary entry a column of a sheet is listed under.
glossary_key <- function(sheet, column) {
  case_when(
    sheet == "arm_cn" & column != "arm" ~ "<sample name>",
    sheet == "cohort_summary" & !column %in% c("metric", "all samples") ~ "<sample_class>",
    sheet == "samples" & column %in% setdiff(manifest_cols, "sex") ~ "<manifest column>",
    .default = column
  )
}

# sheets: every sheet the manifest lists other than itself, by the name
# used in manifest_text. The manifest's own row comes first, because it
# describes the report; the other rows follow the order of the sheets.
make_manifest <- function(sheets) {
  sizes <- tibble(
    sheet = names(sheets),
    rows = map_int(sheets, nrow),
    columns = map_int(sheets, ncol)
  )
  manifest_text |>
    filter(sheet %in% c("manifest", sizes$sheet)) |>
    left_join(sizes, by = "sheet") |>
    mutate(
      rows = if_else(sheet == "manifest", n(), rows),
      columns = if_else(sheet == "manifest", 5L, columns)
    ) |>
    arrange(match(sheet, c("manifest", sizes$sheet))) |>
    select(sheet, one_row_is, rows, columns, contents)
}

make_glossary <- function(sheets) {
  used <- imap(sheets, \(df, name) tibble(sheet = name, column = names(df))) |>
    list_rbind() |>
    mutate(column = glossary_key(sheet, column)) |>
    distinct() |>
    summarize(sheets = str_c(sheet, collapse = ", "), .by = column)
  missing <- setdiff(used$column, glossary_text$column)
  if (length(missing) > 0) {
    warning(
      "Columns without a glossary entry: ", str_c(missing, collapse = ", "),
      call. = FALSE
    )
  }
  used |>
    left_join(glossary_text, by = "column") |>
    mutate(
      description = replace_na(description, "(no definition)"),
      order = match(column, glossary_text$column)
    ) |>
    arrange(order) |>
    select(column, description, sheets)
}

# ---- write -------------------------------------------------------------

log_msg("writing")
in_gene_set <- union(panel$gene, genes_of_interest)

# Sheet order: the genes of interest, the tables a reader can use
# directly, the supporting tables, the full summary, then the manifest
# and the glossary.
sheets <- list(
  genes_of_interest_events = if (is.null(goi_events)) NULL else add_class(goi_events),
  samples = samples_sheet,
  cohort_summary = cohort_summary,
  somatic_variants = variants |> add_class() |> arrange(gene, sample_name),
  disruptions = disruptions |> add_class(),
  gene_cn = gene_cn_events |>
    filter(!is.na(cn_event) | gene %in% in_gene_set) |>
    add_class(),
  fusions = fusions |> add_class(),
  viruses = viruses |> add_class(),
  drivers = drivers |> add_class() |> arrange(gene, sample_name),
  arm_cn = arms,
  qc = if (nrow(qc) > 0) add_class(qc) else qc,
  signatures = if (nrow(signatures) > 0) add_class(signatures) else signatures,
  hla = if (nrow(hla) > 0) add_class(hla) else hla,
  summary_full = summary_full |> rename(sample_group = group)
) |>
  compact()

# Breakpoints shared between samples are review material, kept out of
# the report in a workbook of their own.
review_sheets <- list(sv_same_breakpoint = sv_same_breakpoint |> add_class())

# The glossary covers the columns of both workbooks, the manifest's and
# its own.
manifest_cols_tbl <- tibble(
  sheet = character(), one_row_is = character(), rows = integer(),
  columns = integer(), contents = character()
)
glossary_cols_tbl <- tibble(
  column = character(), description = character(), sheets = character()
)
glossary <- make_glossary(c(
  sheets,
  list(manifest = manifest_cols_tbl, glossary = glossary_cols_tbl),
  set_names(review_sheets, review_file)
))
manifest <- make_manifest(c(
  sheets,
  list(glossary = glossary),
  set_names(review_sheets, str_glue("{review_file} (separate file)"))
))
sheets <- c(sheets, list(manifest = manifest, glossary = glossary))

csv_tables <- c(sheets, review_sheets)
csv_tables$gene_cn <- gene_cn_events |> add_class()
iwalk(csv_tables, \(df, name) {
  write_csv(df, path(data_dir, str_c(name, ".csv")), na = "")
})

header_style <- createStyle(textDecoration = "bold", fgFill = "#D9D9D9")
wrap_style <- createStyle(wrapText = TRUE, valign = "top")
# Column widths of the sheets that hold text.
TEXT_WIDTHS <- list(
  manifest = c(30, 34, 9, 9, 120),
  glossary = c(32, 110, 60)
)
# Widths of the text columns of genes_of_interest_events.
COL_WIDTHS <- c(event = 46, location = 50, support = 72, caution = 34)

write_workbook <- function(sheets, file) {
  wb <- createWorkbook()
  modifyBaseFont(wb, fontSize = 10, fontName = "Arial")
  iwalk(sheets, \(df, name) {
    addWorksheet(wb, name)
    writeData(wb, name, df, headerStyle = header_style)
    freezePane(wb, name, firstRow = TRUE)
    if (name %in% names(TEXT_WIDTHS)) {
      setColWidths(wb, name, cols = seq_len(ncol(df)), widths = TEXT_WIDTHS[[name]])
      addStyle(
        wb, name, wrap_style, rows = seq_len(nrow(df)) + 1,
        cols = seq_len(ncol(df)), gridExpand = TRUE
      )
    } else {
      if (nrow(df) > 0) {
        addFilter(wb, name, rows = 1, cols = seq_len(ncol(df)))
      }
      setColWidths(
        wb, name, cols = seq_len(ncol(df)),
        widths = coalesce(
          unname(COL_WIDTHS[names(df)]),
          pmin(pmax(nchar(names(df)) + 3, 10), 40)
        )
      )
    }
    log_msg("  {name}: {nrow(df)} rows")
  })
  # Explicit window size, so Excel does not open the workbook full screen.
  wb$workbook$bookViews <- str_c(
    '<bookViews><workbookView xWindow="4000" yWindow="3000" ',
    'windowWidth="32000" windowHeight="20000" firstSheet="0" ',
    'activeTab="0"/></bookViews>'
  )
  saveWorkbook(wb, file, overwrite = TRUE)
}

report <- path(res_dir, "prism_report.xlsx")
write_workbook(sheets, report)
log_msg("Report: {report}")
review <- path(res_dir, "prism_sv_same_breakpoint.xlsx")
write_workbook(review_sheets, review)
log_msg("Shared breakpoints: {review}")
log_msg("Tables: {data_dir}")
