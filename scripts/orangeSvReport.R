# ORANGE SV report: every SV record of the ORANGE JSON (section linx:
# fusions, somaticBreakends, germlineBreakends) of every sample, with the
# fields as ORANGE wrote them, plus columns that link each record to its
# event in the SV table and in the SV report (VCF FILTER, flags,
# high_confidence). Description: Prism/docs/SV_REPORT.md section 13.
#
#   Rscript --vanilla Prism/scripts/orangeSvReport.R [OUT_DIR] [DB_DIR] [SV_DIR]
#
# OUT_DIR (oncoanalyser output) defaults to out, DB_DIR to res/prism/db,
# SV_DIR to res/prism/sv. Reads OUT_DIR/<group>/orange/<sample>.orange.json,
# DB_DIR/samples.tsv, DB_DIR/sv_cohort_filtered.parquet, the stacked LINX
# tables in DB_DIR/data, and SV_DIR/prism_sv_events.csv (run svReport.R
# first). Writes SV_DIR/orange_sv_report.xlsx, sheets: fusions,
# somatic_breakends, germline_breakends, samples, about, glossary.
#
# What ORANGE 5.0.1 selects (source at tag orange-v5.0.1):
#   fusions             LINX fusion.tsv rows with reported true
#                       (LinxDataLoader.loadFusions)
#   somaticBreakends    LINX breakend.tsv rows with reportedStatus REPORTED
#                       (LinxBreakendInterpreter.buildSomaticBreakends)
#   germlineBreakends   LINX germline breakend.tsv rows with reportedStatus
#                       REPORTED (buildGermlineBreakends); null when there
#                       is no germline LINX run (tumor-only groups)
# ORANGE tests no SV FILTER value: its SVs are those LINX analyzed.
# The run stops if the JSON and the LINX tables disagree on any of this.

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
})

args <- commandArgs(trailingOnly = TRUE)
oa_dir <- if (length(args) > 0) args[1] else "out"
db_dir <- if (length(args) > 1) args[2] else "res/prism/db"
sv_dir <- if (length(args) > 2) args[3] else "res/prism/sv"

fail <- \(...) stop(..., call. = FALSE)
read_tab <- \(f) read_tsv(path(db_dir, "data", f), col_types = cols(.default = "c"), na = character(),
                          show_col_types = FALSE, progress = FALSE)
join_set <- \(x) if (length(x) == 0) NA_character_ else str_c(sort(unique(x)), collapse = ";")

# ---- ORANGE JSON -------------------------------------------------------

samples <- read_tsv(path(db_dir, "samples.tsv"), col_types = cols(.default = "c"), na = character(),
                    show_col_types = FALSE) |>
  mutate(has_normal = normal_id != "", json = path(oa_dir, group_id, "orange", str_c(sample_id, ".orange.json")))
absent <- samples$json[!file_exists(samples$json)]
if (length(absent) > 0) fail("ORANGE JSON not found: ", str_c(absent, collapse = ", "))

# One table per JSON list. An all-null field is read as logical NA and
# made character, so that the samples bind; rnaSupport (an object, null
# without RNA) becomes two columns.
as_table <- function(x) {
  if (length(x) == 0) return(NULL)
  x <- as_tibble(x)
  if ("rnaSupport" %in% names(x)) {
    if (!all(is.na(x$rnaSupport))) fail("rnaSupport is not null and not an object")
    x <- x |> select(-rnaSupport) |> mutate(rnaSupport.alleleReadCount = NA_integer_, rnaSupport.totalReadCount = NA_integer_)
  }
  x |> mutate(across(where(\(v) is.logical(v) && all(is.na(v))), as.character))
}
js <- samples |>
  pmap(\(group_id, sample_id, has_normal, json, ...) {
    j <- jsonlite::fromJSON(json, flatten = TRUE)
    if (j$sampleId != sample_id) fail(json, ": sampleId is ", j$sampleId)
    germline_null <- is.null(j$linx$germlineBreakends)
    if (germline_null == has_normal) fail(json, ": germlineBreakends is ", if (germline_null) "null" else "present",
                                          " for a group ", if (has_normal) "with" else "without", " a normal")
    key <- \(x) if (is.null(x)) NULL else x |> mutate(group_id = group_id, sample_id = sample_id, .before = 1)
    list(fusions = key(as_table(j$linx$fusions)),
         somatic = key(as_table(j$linx$somaticBreakends)),
         germline = key(as_table(j$linx$germlineBreakends)))
  })
pick_part <- \(p) map(js, p) |> list_rbind()
fusions <- pick_part("fusions") |> rename_with(\(x) str_replace(x, "^rnaSupport\\.", "rnaSupport_"))
somatic <- pick_part("somatic")
germline <- pick_part("germline")

# ---- the same records in the LINX tables -------------------------------

lx_fusion <- read_tab("linx_somatic_fusion.tsv.gz")
lx_bnd <- read_tab("linx_somatic_breakend.tsv.gz")
lx_gbnd <- read_tab("linx_germline_breakend.tsv.gz")
lx_gdis <- read_tab("linx_germline_disruption.tsv.gz")

# Each ORANGE breakend is the REPORTED breakend.tsv row with its id, and
# no REPORTED row is left out.
check_breakends <- function(orange, lx, label) {
  rep <- lx |> filter(reportedStatus == "REPORTED") |>
    transmute(sample_id, id = as.integer(id), svId = as.integer(svId), gene, transcript = transcriptId)
  o <- if (is.null(orange)) rep[0, ] else orange |> select(sample_id, id, svId, gene, transcript)
  if (nrow(o) != nrow(rep) || nrow(anti_join(o, rep, by = names(rep))) > 0 || nrow(anti_join(rep, o, by = names(rep))) > 0) {
    fail(label, ": ORANGE breakends are not the LINX breakend.tsv rows with reportedStatus REPORTED")
  }
}
check_breakends(somatic, lx_bnd, "somaticBreakends")
check_breakends(germline, lx_gbnd |> filter(sample_id %in% samples$sample_id[samples$has_normal]), "germlineBreakends")

# Each ORANGE fusion is one reported fusion.tsv row (matched on the
# breakend transcripts, fused exons, type, phase and chain links), and no
# reported row is left out. geneUp is the 5' gene, or for enhancer types
# the 5' transcript ID (the IG or TCR region name); geneDown the 3' gene.
bnd_key <- lx_bnd |> select(sample_id, id, svId, gene, transcriptId)
lx_rep <- lx_fusion |>
  filter(reported == "true") |>
  left_join(bnd_key |> rename_with(\(x) str_c(x, "_5p"), -sample_id), by = c("sample_id", fivePrimeBreakendId = "id_5p")) |>
  left_join(bnd_key |> rename_with(\(x) str_c(x, "_3p"), -sample_id), by = c("sample_id", threePrimeBreakendId = "id_3p"))
if (anyNA(lx_rep$svId_5p) || anyNA(lx_rep$svId_3p)) fail("reported fusion with a breakend not in breakend.tsv")
fusion_by <- c("sample_id", "transcriptUp", "transcriptDown", "fusedExonUp", "fusedExonDown", "reportedType", "phased", "chainLinks")
if (is.null(fusions) || is.null(somatic)) fail("no ORANGE fusions or no somatic breakends in any sample: not handled")
fusions <- fusions |>
  mutate(.row = row_number()) |>
  left_join(
    lx_rep |> transmute(sample_id, transcriptUp = transcriptId_5p, transcriptDown = transcriptId_3p,
                        fusedExonUp = as.integer(fusedExonUp), fusedExonDown = as.integer(fusedExonDown),
                        reportedType, phased, chainLinks = as.integer(chainLinks),
                        lx_gene_5p = gene_5p, lx_gene_3p = gene_3p, svId_5p = as.integer(svId_5p), svId_3p = as.integer(svId_3p),
                        fivePrimeCoords, threePrimeCoords),
    by = fusion_by, relationship = "many-to-many"
  )
if (anyDuplicated(fusions$.row) || anyNA(fusions$svId_5p) || nrow(fusions) != nrow(lx_rep)) {
  fail("ORANGE fusions do not match the reported LINX fusion.tsv rows one to one")
}
enhancer <- fusions$reportedType %in% c("ENHANCER_KNOWN_PAIR", "ENHANCER_PROMISCUOUS")
if (any(fusions$geneUp != if_else(enhancer, fusions$transcriptUp, fusions$lx_gene_5p)) ||
    any(fusions$geneDown != fusions$lx_gene_3p)) {
  fail("ORANGE fusion geneUp/geneDown differ from the LINX breakend genes")
}

# ---- link to the SV table and the SV report -----------------------------

sv <- nanoparquet::read_parquet(
  path(db_dir, "sv_cohort_filtered.parquet"),
  col_select = c("sample_id", "event_key", "filter", "linx_svId", "start_chrom", "start_pos", "start_orient",
                 "end_chrom", "end_pos", "end_orient")
) |> as_tibble() |> filter(!is.na(linx_svId))
if (anyDuplicated(sv[c("sample_id", "linx_svId")])) fail("LINX svId not unique within a sample in the SV table")
events <- read_csv(path(sv_dir, "prism_sv_events.csv"), show_col_types = FALSE, guess_max = 1e5)
if (!setequal(events$event_key, sv$event_key[sv$filter == "PASS"])) {
  fail(path(sv_dir, "prism_sv_events.csv"), " does not hold the PASS events of the SV table: rerun svReport.R")
}
flag_cols <- c("flag_own_normal", "flag_other_patients", "flag_population", "flag_linx_artifact", "flag_region")
report_cols <- events |>
  transmute(event_key, in_sv_report = TRUE, high_confidence, n_flags, across(all_of(flag_cols)),
            flags = pmap_chr(pick(all_of(flag_cols)), \(...) join_set(str_remove(flag_cols[c(...)], "^flag_"))))

to_event <- \(df, svid) {
  df |>
    left_join(sv |> select(sample_id, linx_svId, event_key, vcf_filter = filter, start_chrom:end_orient),
              by = c("sample_id", setNames("linx_svId", svid)), relationship = "many-to-one")
}

somatic_out <- somatic |>
  to_event("svId") |>
  mutate(breakend_side = case_when(
    chromosome == start_chrom & position == start_pos & orientation == start_orient ~ "start",
    chromosome == end_chrom & position == end_pos & orientation == end_orient ~ "end"
  )) |>
  left_join(report_cols, by = "event_key", relationship = "many-to-one") |>
  mutate(in_sv_report = coalesce(in_sv_report, FALSE)) |>
  select(-(start_chrom:end_orient), -flags)
if (anyNA(somatic_out$event_key)) fail("somatic breakend whose LINX SV is not in the SV table")
if (anyNA(somatic_out$breakend_side)) fail("somatic breakend position differs from both ends of its SV table event")

end_cols <- \(df, end) {
  df |>
    left_join(sv |> select(sample_id, linx_svId, event_key, vcf_filter = filter),
              by = c("sample_id", setNames("linx_svId", str_c("svId_", end))), relationship = "many-to-one") |>
    left_join(report_cols |> select(event_key, in_sv_report, high_confidence, flags), by = "event_key",
              relationship = "many-to-one") |>
    mutate(in_sv_report = coalesce(in_sv_report, FALSE)) |>
    rename_with(\(x) str_c(x, "_", end), c(event_key, vcf_filter, in_sv_report, high_confidence, flags))
}
fusions_out <- fusions |>
  end_cols("5p") |>
  end_cols("3p") |>
  arrange(.row) |>
  select(-.row, -lx_gene_5p, -lx_gene_3p, -svId_5p, -svId_3p)
if (anyNA(fusions_out$event_key_5p) || anyNA(fusions_out$event_key_3p)) fail("fusion breakend whose LINX SV is not in the SV table")

germline_out <- germline
if (!is.null(germline)) {
  germline_out <- germline |>
    left_join(lx_gdis |>
                transmute(sample_id, svId = as.integer(svId), gene, chromosomeStart, positionStart = as.integer(positionStart),
                          orientStart = as.integer(orientStart), chromosomeEnd, positionEnd = as.integer(positionEnd),
                          orientEnd = as.integer(orientEnd), vcf_filter = filter, qualScore = as.numeric(qualScore),
                          germlineFragments = as.integer(germlineFragments),
                          germlineReferenceFragmentsStart = as.integer(germlineReferenceFragmentsStart),
                          germlineReferenceFragmentsEnd = as.integer(germlineReferenceFragmentsEnd),
                          tumorFragments = as.integer(tumorFragments),
                          tumorReferenceFragmentsStart = as.integer(tumorReferenceFragmentsStart),
                          tumorReferenceFragmentsEnd = as.integer(tumorReferenceFragmentsEnd)),
              by = c("sample_id", "svId", "gene"), relationship = "many-to-one")
  if (anyNA(germline_out$chromosomeStart)) fail("germline breakend without its LINX germline disruption.tsv row")
  at_start <- with(germline_out, chromosome == chromosomeStart & position == positionStart & orientation == orientStart)
  at_end <- with(germline_out, chromosome == chromosomeEnd & position == positionEnd & orientation == orientEnd)
  if (!all(at_start | at_end)) fail("germline breakend position differs from both ends of its LINX germline SV")
}
if (is.null(germline_out)) germline_out <- tibble(group_id = character(), sample_id = character())

# ---- samples and about ------------------------------------------------

n_by <- \(df, name) if (nrow(df) == 0) tibble(sample_id = character(), !!name := integer()) else count(df, sample_id, name = name)
samples_out <- samples |>
  select(group_id, sample_id, has_normal) |>
  left_join(n_by(fusions_out, "fusions"), by = "sample_id") |>
  left_join(n_by(somatic_out, "somatic_breakends"), by = "sample_id") |>
  left_join(somatic_out |> distinct(sample_id, event_key, high_confidence) |>
              summarize(somatic_events = n(), somatic_events_hc = sum(high_confidence %in% TRUE), .by = sample_id),
            by = "sample_id") |>
  left_join(n_by(germline_out, "germline_breakends"), by = "sample_id") |>
  left_join(count(events, sample_id, name = "sv_report_pass_events"), by = "sample_id") |>
  mutate(across(c(fusions, somatic_breakends, somatic_events, somatic_events_hc), \(x) replace_na(x, 0L)),
         germline_breakends = if_else(has_normal, replace_na(germline_breakends, 0L), NA_integer_))

som_ev <- somatic_out |> distinct(event_key, vcf_filter, in_sv_report, high_confidence)
fus_ev <- bind_rows(
  fusions_out |> select(event_key = event_key_5p, vcf_filter = vcf_filter_5p, in_sv_report = in_sv_report_5p,
                        high_confidence = high_confidence_5p),
  fusions_out |> select(event_key = event_key_3p, vcf_filter = vcf_filter_3p, in_sv_report = in_sv_report_3p,
                        high_confidence = high_confidence_3p)
) |> distinct()
about <- tribble(
  ~item, ~definition, ~records, ~events, ~events_in_sv_report, ~events_high_confidence,
  "source", "ORANGE 5.0.1 JSON, OUT_DIR/<group>/orange/<sample>.orange.json, section linx. One row per record, with the fields as ORANGE wrote them; the columns after them are added by this report (see glossary). Source code read at tag orange-v5.0.1.", NA, NA, NA, NA,
  "fusions", "LINX somatic fusion.tsv rows with reported true (LinxDataLoader.loadFusions). The run checks that the JSON holds exactly these rows.", nrow(fusions_out), nrow(fus_ev), sum(fus_ev$in_sv_report), sum(fus_ev$high_confidence %in% TRUE),
  "somatic_breakends", "LINX somatic breakend.tsv rows with reportedStatus REPORTED (LinxBreakendInterpreter.buildSomaticBreakends). LINX sets REPORTED only on a disruptive breakend in a driver gene panel gene with disruption reporting on, on its canonical or listed transcript. The run checks that the JSON holds exactly these rows.", nrow(somatic_out), nrow(som_ev), sum(som_ev$in_sv_report), sum(som_ev$high_confidence %in% TRUE),
  "germline_breakends", "LINX germline breakend.tsv rows with reportedStatus REPORTED (LinxBreakendInterpreter.buildGermlineBreakends). Null in the JSON of tumor-only groups, which have no germline LINX run. Germline SVs are not in the SV table or the SV report.", nrow(germline_out), NA, NA, NA,
  "vcf_filter", as.character(str_glue("ORANGE tests no SV FILTER value: its SVs are those LINX analyzed, somatic FILTER PASS or INFERRED (PON records are not analyzed) and germline PASS or PON. LINX never marks a breakend of an inferred SV disruptive, so it is never REPORTED. This run: events of somatic breakends and fusions with a FILTER other than PASS: {n_distinct(c(som_ev$event_key[som_ev$vcf_filter != 'PASS'], fus_ev$event_key[fus_ev$vcf_filter != 'PASS']))}.")), NA, NA, NA, NA,
  "link", "Somatic breakends and fusion ends are placed on SV table events by sample and LINX svId (each SV table event analyzed by LINX has one svId per sample). Fusions are matched to their LINX fusion.tsv row on sample, transcripts, fused exons, reported type, phase and chain links; the 5' and 3' breakends of that row give the two events (the same event if one SV makes the fusion). Germline breakends are matched to the LINX germline disruption.tsv row on sample, svId and gene.", NA, NA, NA, NA
)

# ---- write ---------------------------------------------------------------

sheets <- list(
  fusions = fusions_out,
  somatic_breakends = somatic_out,
  germline_breakends = germline_out,
  samples = samples_out,
  about = about
)

glossary_text <- c(
  # keys
  group_id = "oncoanalyser group (one tumor, with its normal if any).",
  sample_id = "Tumor sample ID (ORANGE sampleId).",
  # ORANGE fusion fields
  geneUp = "ORANGE: 5' gene (gene of the fusion's 5' breakend in LINX breakend.tsv); for reportedType ENHANCER_KNOWN_PAIR or ENHANCER_PROMISCUOUS, the 5' breakend's transcript ID, which LINX sets to the IG or TCR region name (e.g. @IGH).",
  clusterId = "ORANGE: LINX cluster ID (svs.tsv clusterId) of the breakend's SV; in fusions, of the 5' breakend's SV.",
  contextUp = "ORANGE: where the 5' breakend is: 'Exon N' (N = fusedExonUp) for an exonic or intronic breakend, 'Promoter Region' upstream of the gene, 'Post-coding' downstream, 'IG' in an IG region.",
  transcriptUp = "ORANGE: transcript of the 5' breakend (breakend.tsv transcriptId).",
  geneDown = "ORANGE: 3' gene (gene of the fusion's 3' breakend in LINX breakend.tsv).",
  contextDown = "ORANGE: as contextUp, for the 3' breakend (N = fusedExonDown).",
  transcriptDown = "ORANGE: transcript of the 3' breakend (breakend.tsv transcriptId).",
  reportedType = "ORANGE: LINX fusion.tsv reportedType: the fusion's category in the known fusion file (KNOWN_PAIR, PROMISCUOUS_5, PROMISCUOUS_3, ENHANCER_PROMISCUOUS and others).",
  phased = "ORANGE: LINX fusion.tsv phased: INFRAME, SKIPPED_EXONS or OUT_OF_FRAME.",
  driverInterpretation = "ORANGE: LINX fusion.tsv likelihood HIGH or LOW; UNKNOWN for any other value.",
  fusedExonUp = "ORANGE: LINX fusion.tsv fusedExonUp: last exon of the 5' gene in the fusion.",
  fusedExonDown = "ORANGE: LINX fusion.tsv fusedExonDown: first exon of the 3' gene in the fusion.",
  chainLinks = "ORANGE: LINX fusion.tsv chainLinks: number of links in the chain of SVs joining the two genes (0: one SV).",
  chainTerminated = "ORANGE: LINX fusion.tsv chainTerminated: TRUE if another breakend in the chain interrupts the fusion transcript.",
  domainsKept = "ORANGE: LINX fusion.tsv domainsKept: protein domains of the 3' gene kept in the fusion.",
  domainsLost = "ORANGE: LINX fusion.tsv domainsLost: protein domains of the 3' gene lost in the fusion.",
  junctionCopyNumber = "ORANGE: breakends: mean of the SV's LINX svs.tsv junctionCopyNumberMin and junctionCopyNumberMax; fusions: mean of that value over the 5' and 3' breakends' SVs.",
  rnaSupport_alleleReadCount = "ORANGE rnaSupport.alleleReadCount: RNA fragments supporting the fusion (Isofox). Empty: no RNA data (null in the JSON).",
  rnaSupport_totalReadCount = "ORANGE rnaSupport.totalReadCount: RNA depth at the fusion (Isofox). Empty: no RNA data (null in the JSON).",
  plotFilename = "ORANGE: LINX plot of the record's cluster, if LINX wrote one; empty for germline breakends.",
  # added to fusions
  fivePrimeCoords = "Added: LINX fusion.tsv fivePrimeCoords of the matched row: 5' breakend as chromosome:position:orientation.",
  threePrimeCoords = "Added: LINX fusion.tsv threePrimeCoords of the matched row: 3' breakend as chromosome:position:orientation.",
  event_key_5p = "Added: event_key of the SV table event (sample_id:VCF ID of its start breakend) holding the fusion's 5' breakend.",
  vcf_filter_5p = "Added: PURPLE somatic SV VCF FILTER of that event (SV table filter).",
  in_sv_report_5p = "Added: TRUE if that event is in the SV report (prism_sv_report.xlsx, PASS events only).",
  high_confidence_5p = "Added: SV report high_confidence of that event (no flag set); empty if not in the SV report.",
  flags_5p = "Added: SV report flags set on that event, without the flag_ prefix, joined with ;. Empty: no flag, or not in the SV report.",
  event_key_3p = "Added: as event_key_5p, for the 3' breakend.",
  vcf_filter_3p = "Added: as vcf_filter_5p, for the 3' breakend's event.",
  in_sv_report_3p = "Added: as in_sv_report_5p, for the 3' breakend's event.",
  high_confidence_3p = "Added: as high_confidence_5p, for the 3' breakend's event.",
  flags_3p = "Added: as flags_5p, for the 3' breakend's event.",
  # ORANGE breakend fields
  id = "ORANGE: LINX breakend.tsv id (one row per breakend and transcript).",
  svId = "ORANGE: LINX SV ID (svs.tsv svId), unique within a sample.",
  gene = "ORANGE: LINX breakend.tsv gene.",
  position = "ORANGE: breakend position, from LINX svs.tsv coordsStart or coordsEnd of its SV.",
  orientation = "ORANGE: breakend orientation (1 or -1), from LINX svs.tsv coordsStart or coordsEnd.",
  chromosome = "ORANGE: breakend chromosome, from LINX svs.tsv coordsStart or coordsEnd.",
  chromosomeBand = "ORANGE: cytoband at the breakend (ORANGE look-up).",
  transcript = "ORANGE: LINX breakend.tsv transcriptId.",
  isCanonical = "ORANGE: LINX breakend.tsv canonical.",
  geneOrientation = "ORANGE: LINX breakend.tsv geneOrientation: UPSTREAM (the 5' part of the gene stays at the breakend) or DOWNSTREAM.",
  disruptive = "ORANGE: LINX breakend.tsv disruptive.",
  reportedStatus = "ORANGE: LINX breakend.tsv reportedStatus; always REPORTED, because ORANGE lists only those rows.",
  undisruptedCopyNumber = "ORANGE: LINX breakend.tsv undisruptedCopyNumber: copies of the gene this breakend leaves intact.",
  type = "ORANGE: LINX svs.tsv type of the SV: DEL, DUP, INV, INS, BND, SGL or INF.",
  regionType = "ORANGE: LINX breakend.tsv regionType: INTRONIC, EXONIC, UPSTREAM and others.",
  codingType = "ORANGE: LINX breakend.tsv codingType: CODING, UTR_5P, UTR_3P, NON_CODING, ENHANCER.",
  nextSpliceExonRank = "ORANGE: LINX breakend.tsv nextSpliceExonRank.",
  exonUp = "ORANGE: LINX breakend.tsv exonUp: nearest exon upstream of the breakend in the transcript.",
  exonDown = "ORANGE: LINX breakend.tsv exonDown: nearest exon downstream of the breakend in the transcript.",
  driverType = "ORANGE: driver type of the first LINX driver catalog entry with the same gene (somatic: LINX somatic driver types; germline: LINX germline catalog); DISRUPTION if there is none or it is GERMLINE_DISRUPTION.",
  driverLikelihood = "ORANGE: driverLikelihood of that driver catalog entry; 0 if there is none.",
  # added to somatic breakends
  event_key = "Added: event_key of the SV table event (sample_id:VCF ID of its start breakend) of the breakend's SV.",
  vcf_filter = "Added. somatic_breakends: PURPLE somatic SV VCF FILTER of the event (SV table filter). germline_breakends: LINX germline disruption.tsv filter of the SV.",
  breakend_side = "Added: which end of the SV table event the breakend is (start or end), by chromosome, position and orientation.",
  in_sv_report = "Added: TRUE if the event is in the SV report (prism_sv_report.xlsx, PASS events only).",
  high_confidence = "Added: SV report high_confidence (no flag set); empty if not in the SV report.",
  n_flags = "Added: SV report n_flags (0 to 5); empty if not in the SV report.",
  flag_own_normal = "Added: SV report flag_own_normal (definition: prism_sv_report.xlsx sheet rule).",
  flag_other_patients = "Added: SV report flag_other_patients (definition: prism_sv_report.xlsx sheet rule).",
  flag_population = "Added: SV report flag_population (definition: prism_sv_report.xlsx sheet rule).",
  flag_linx_artifact = "Added: SV report flag_linx_artifact (definition: prism_sv_report.xlsx sheet rule).",
  flag_region = "Added: SV report flag_region (definition: prism_sv_report.xlsx sheet rule).",
  # added to germline breakends (LINX germline disruption.tsv)
  chromosomeStart = "Added: LINX germline disruption.tsv chromosomeStart of the SV.",
  positionStart = "Added: LINX germline disruption.tsv positionStart.",
  orientStart = "Added: LINX germline disruption.tsv orientStart.",
  chromosomeEnd = "Added: LINX germline disruption.tsv chromosomeEnd.",
  positionEnd = "Added: LINX germline disruption.tsv positionEnd.",
  orientEnd = "Added: LINX germline disruption.tsv orientEnd.",
  qualScore = "Added: LINX germline disruption.tsv qualScore.",
  germlineFragments = "Added: LINX germline disruption.tsv germlineFragments: fragments supporting the SV in the normal.",
  germlineReferenceFragmentsStart = "Added: LINX germline disruption.tsv germlineReferenceFragmentsStart: reference fragments in the normal at the start breakend.",
  germlineReferenceFragmentsEnd = "Added: LINX germline disruption.tsv germlineReferenceFragmentsEnd.",
  tumorFragments = "Added: LINX germline disruption.tsv tumorFragments: fragments supporting the SV in the tumor.",
  tumorReferenceFragmentsStart = "Added: LINX germline disruption.tsv tumorReferenceFragmentsStart.",
  tumorReferenceFragmentsEnd = "Added: LINX germline disruption.tsv tumorReferenceFragmentsEnd.",
  # samples
  has_normal = "TRUE if the group has a matched normal (germline records exist only then).",
  fusions = "Number of ORANGE fusions of the sample.",
  somatic_breakends = "Number of ORANGE somatic breakend records of the sample.",
  somatic_events = "Number of SV table events behind those records.",
  somatic_events_hc = "Number of those events that are high confidence in the SV report.",
  germline_breakends = "Number of ORANGE germline breakend records; empty for tumor-only groups (null in the JSON).",
  sv_report_pass_events = "Number of PASS events of the sample in the SV report, for scale.",
  # about
  item = "Record type or term described on that row.",
  definition = "Definition of the item.",
  records = "Number of ORANGE records of that type, all samples.",
  events = "Number of distinct SV table events behind them.",
  events_in_sv_report = "Number of those events in the SV report.",
  events_high_confidence = "Number of those events that are high confidence in the SV report.",
  # glossary
  columnName = "Column name.",
  sheets = "Sheets that have the column.",
  description = "What the column holds. 'ORANGE:' a field of the ORANGE JSON as written (source at orange-v5.0.1); 'Added:' a column this report adds."
)
glossary_cols <- c("columnName", "sheets", "description")
glossary <- imap(c(sheets, list(glossary = set_names(vector("list", 3), glossary_cols))),
                 \(d, nm) tibble(columnName = names(d), sheet = nm)) |>
  bind_rows() |>
  group_by(columnName) |>
  summarize(sheets = str_c(sheet, collapse = ", "), first = min(cur_group_rows()), .groups = "drop") |>
  arrange(first) |>
  select(-first)
missing <- setdiff(glossary$columnName, names(glossary_text))
unused <- setdiff(names(glossary_text), glossary$columnName)
if (length(missing) > 0) fail("glossary: no description for ", str_c(missing, collapse = ", "))
if (length(unused) > 0) fail("glossary: description for no column: ", str_c(unused, collapse = ", "))
sheets$glossary <- glossary |> mutate(description = unname(glossary_text[columnName]))

wb <- openxlsx::createWorkbook()
openxlsx::modifyBaseFont(wb, fontName = "Arial", fontSize = 10)
for (nm in names(sheets)) {
  openxlsx::addWorksheet(wb, nm)
  openxlsx::writeData(wb, nm, sheets[[nm]], withFilter = !nm %in% c("about", "glossary"),
                      headerStyle = openxlsx::createStyle(textDecoration = "bold"))
  openxlsx::freezePane(wb, nm, firstRow = TRUE)
  openxlsx::setColWidths(wb, nm, cols = seq_along(sheets[[nm]]), widths = "auto")
}
for (nm in c("about", "glossary")) {
  openxlsx::setColWidths(wb, nm, cols = 2 + (nm == "glossary"), widths = 100)
  openxlsx::addStyle(wb, nm, openxlsx::createStyle(wrapText = TRUE, valign = "top"),
                     rows = seq_len(nrow(sheets[[nm]]) + 1), cols = seq_along(sheets[[nm]]), gridExpand = TRUE)
}
# Explicit window size: without it Excel for Mac opens the file full screen.
wb$workbook$bookViews <- str_c('<bookViews><workbookView xWindow="4000" yWindow="3000" ',
                               'windowWidth="32000" windowHeight="20000"/></bookViews>')
xlsx_file <- path(sv_dir, "orange_sv_report.xlsx")
openxlsx::saveWorkbook(wb, xlsx_file, overwrite = TRUE)

print(about |> select(item, records:events_high_confidence))
message("-> ", xlsx_file)
