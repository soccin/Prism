# SV report for people: the PASS somatic SV events of the cohort with
# LINX's genes and fusions, five evidence flags, and a high-confidence
# set (no flag). Reads only the SV table that postOA.sh writes
# (sv_cohort_filtered.parquet). Full description, with the LINX rules
# behind each gene column: Prism/docs/SV_REPORT.md.
#
#   Rscript --vanilla Prism/scripts/svReport.R [DB_DIR] [OUT_DIR]
#
# DB_DIR defaults to res/prism/db, OUT_DIR to res/prism/sv. Writes
#   prism_sv_report.xlsx   sheets: fusions_real, genes, events_hc,
#                          events_all, rule, glossary (every column of
#                          every sheet, once per column name)
#
# User decisions of 2026-11-23: a LINX fusion with the same gene at both
# ends joins two parts of one gene (treated by the user as a cryptic
# splice product), not two genes, and is left out of
# fusions_real and of the fusion counts of genes (kept on the events
# sheets); the genes sheet lists only genes with 2 or more tumors in
# at least one count; fusions_real and genes have a samples column.
#   prism_sv_events.csv    the events_all sheet
#
# Genes come from LINX only: start_gene/end_gene are LINX svs.tsv
# geneStart/geneEnd (genes whose span contains the breakend, or with a
# transcript LINX kept starting 1 to 10,000 bases after it in the
# direction of transcription); *_disrupted_gene: genes whose canonical
# transcript LINX marks disruptive at that breakend; *_driver_disruption:
# genes LINX reports as driver disruptions at that breakend. Fusions are
# LINX fusion.tsv names, each listed once per event.
#
# Flags (user decision, 2026-11-23; a flag is evidence, not proof):
#   flag_own_normal      split reads for the break in the group's own
#                        normal BAM (bam_normal_split_n > 0, either end);
#                        never set for tumor-only groups (has_normal FALSE)
#   flag_other_patients  the same event in another group: two-breakend
#                        events, both breakends within 50 bases (same
#                        orientations, one SV) in another group's somatic
#                        calls (PURPLE SV VCF records other than INFERRED,
#                        so PASS and PON) or germline calls
#                        (xs_pair_*_groups_w50 > 0); one-breakend events,
#                        the breakend within 50 bases
#                        (start_xs_*_groups_w50 > 0); or, any event, split
#                        reads for the break in another group's normal or
#                        tumor BAM (bam_other_*_split_samples > 0, either
#                        end)
#   flag_population      a population SV site matching both breakpoints
#                        within 100 bases (gnomAD-SV v4.1, 1000 Genomes,
#                        dbVar nstd166, nstd200; any FILTER; *_bp100_n > 0)
#   flag_linx_artifact   LINX cluster category ARTIFACT
#   flag_region          a breakend in GIAB low mappability, the ENCODE
#                        blacklist, the Hartwig SV blacklist or Hartwig
#                        unmappable regions
# high_confidence: no flag set.

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
})

args <- commandArgs(trailingOnly = TRUE)
db_dir <- if (length(args) > 0) args[1] else "res/prism/db"
out_dir <- if (length(args) > 1) args[2] else "res/prism/sv"
dir_create(out_dir)

sv <- nanoparquet::read_parquet(path(db_dir, "sv_cohort_filtered.parquet")) |> as_tibble()
pass <- sv |> filter(filter == "PASS")
if (!all(pass$in_linx)) stop("PASS events not in LINX", call. = FALSE)
if (anyDuplicated(pass$event_key)) stop("event_key not unique", call. = FALSE)

pos0 <- \(x) !is.na(x) & x > 0
join_set <- \(x) if (length(x) == 0) NA_character_ else str_c(sort(unique(x)), collapse = ";")

# ---- genes per breakend (LINX breakend.tsv lists, aligned by position) ---

bnd_long <- map(c("start", "end"), \(s) {
  pass |>
    select(event_key, gene = str_c(s, "_linx_bnd_gene"), role = str_c(s, "_gene_role"),
           canonical = str_c(s, "_linx_bnd_canonical"), disruptive = str_c(s, "_linx_bnd_disruptive"),
           reported = str_c(s, "_linx_bnd_reportedStatus")) |>
    unnest(c(gene, role, canonical, disruptive, reported)) |>
    mutate(side = s)
}) |>
  bind_rows()

side_genes <- bnd_long |>
  group_by(event_key, side) |>
  summarize(
    disrupted_gene = join_set(gene[canonical & disruptive]),
    driver_disruption = join_set(gene[reported == "REPORTED"]),
    .groups = "drop"
  ) |>
  pivot_wider(names_from = side, values_from = c(disrupted_gene, driver_disruption), names_glue = "{side}_{.value}")

gene_role <- bnd_long |>
  filter(!is.na(role)) |>
  distinct(gene, role)
if (anyDuplicated(gene_role$gene)) stop("a gene with two driver panel roles", call. = FALSE)

# ---- fusions per event, each once ----------------------------------------

fusion_long <- pass |>
  select(group_id, sample_id, event_key, name = linx_fusion_name, reported = linx_fusion_reported,
         phased = linx_fusion_phased, gene_5p = fusion_gene_5p, gene_3p = fusion_gene_3p) |>
  unnest(c(name, reported, phased, gene_5p, gene_3p)) |>
  distinct(group_id, sample_id, event_key, name, .keep_all = TRUE)
if (anyNA(fusion_long$gene_5p) || anyNA(fusion_long$gene_3p)) stop("fusion without a 5' or 3' gene", call. = FALSE)

# The events sheets list every LINX fusion of the event (user decision,
# 2026-11-23).
event_fusions <- fusion_long |>
  group_by(event_key) |>
  summarize(fusions = join_set(name), reported_fusions = join_set(name[reported]), .groups = "drop")

# A LINX fusion whose 5' and 3' breakends are in the same gene joins two
# parts of one gene (an intragenic rearrangement; the user treats these
# as cryptic splice products), not two genes (user decision,
# 2026-11-23): it is left out of fusions_real and of the fusion counts of
# genes (counted in the rule sheet), and kept on the events sheets.
same_gene <- fusion_long |> filter(gene_5p == gene_3p)
fusion_long <- fusion_long |> filter(gene_5p != gene_3p)

# ---- events with flags ---------------------------------------------------

any_of_cols <- \(df, cols) reduce(map(cols, \(c) pos0(df[[c]])), `|`)
in_region <- \(df, s) {
  reduce(map(str_c(s, c("_region_giab_lowmappability", "_region_encode_blacklist",
                         "_region_hmf_sv_prep_blacklist", "_region_hmf_unmap")),
             \(c) coalesce(df[[c]], FALSE)), `|`)
}

events <- pass |>
  mutate(
    has_normal = !is.na(normal_id) & normal_id != "",
    normal_split_reads = pmax(start_bam_normal_split_n, end_bam_normal_split_n, na.rm = TRUE),
    flag_own_normal = pos0(start_bam_normal_split_n) | pos0(end_bam_normal_split_n),
    flag_other_patients = if_else(
      n_breakends == 2,
      pos0(xs_pair_somatic_groups_w50) | pos0(xs_pair_germline_groups_w50),
      pos0(start_xs_somatic_groups_w50) | pos0(start_xs_germline_groups_w50)
    ) |
      any_of_cols(pick(everything()), c("start_bam_other_normals_split_samples", "end_bam_other_normals_split_samples",
                                        "start_bam_other_tumors_split_samples", "end_bam_other_tumors_split_samples")),
    flag_population = any_of_cols(pick(everything()), c("gnomad_sv_bp100_n", "kgp_sv_bp100_n",
                                                        "dbvar_nstd166_bp100_n", "dbvar_nstd200_bp100_n")),
    flag_linx_artifact = coalesce(linx_cluster_category == "ARTIFACT", FALSE),
    flag_region = in_region(pick(everything()), "start") | in_region(pick(everything()), "end"),
    other_patients_with_event = map2_chr(xs_pair_somatic_group_list_w50, xs_pair_germline_group_list_w50,
                                         \(a, b) join_set(c(a, b))),
    start_gene = map_chr(linx_geneStart, join_set),
    end_gene = map_chr(linx_geneEnd, join_set)
  ) |>
  mutate(n_flags = flag_own_normal + flag_other_patients + flag_population + flag_linx_artifact + flag_region,
         high_confidence = n_flags == 0) |>
  left_join(side_genes, by = "event_key", relationship = "one-to-one") |>
  left_join(event_fusions, by = "event_key", relationship = "one-to-one") |>
  transmute(
    group_id, sample_id, event_key, sv_type = linx_type,
    start_chrom, start_pos, start_orient, end_chrom, end_pos, end_orient, length,
    start_gene, start_disrupted_gene, start_driver_disruption,
    end_gene, end_disrupted_gene, end_driver_disruption,
    fusions, reported_fusions,
    cluster_type = linx_cluster_resolvedType, cluster_size = linx_cluster_clusterCount,
    tumor_fragments = tumor_VF, start_af = start_info_PURPLE_AF, end_af = end_info_PURPLE_AF,
    has_normal, normal_split_reads, other_patients_with_event,
    flag_own_normal, flag_other_patients, flag_population, flag_linx_artifact, flag_region,
    n_flags, high_confidence
  )
stopifnot(nrow(events) == nrow(pass), !anyNA(events$high_confidence))

# ---- genes and fusions, counted in tumors (sample_id) --------------------

hc_keys <- events$event_key[events$high_confidence]
gene_hits <- bind_rows(
  events |> select(sample_id, event_key, g = start_gene) |> mutate(kind = "at_breakend"),
  events |> select(sample_id, event_key, g = end_gene) |> mutate(kind = "at_breakend"),
  events |> select(sample_id, event_key, g = start_disrupted_gene) |> mutate(kind = "disrupted"),
  events |> select(sample_id, event_key, g = end_disrupted_gene) |> mutate(kind = "disrupted"),
  fusion_long |> select(sample_id, event_key, g = gene_5p) |> mutate(kind = "fusion"),
  fusion_long |> select(sample_id, event_key, g = gene_3p) |> mutate(kind = "fusion")
) |>
  filter(!is.na(g)) |>
  separate_longer_delim(g, ";") |>
  distinct() |>
  mutate(hc = event_key %in% hc_keys)

count_tumors <- \(d, k, only_hc) n_distinct(d$sample_id[d$kind == k & (d$hc | !only_hc)])
genes <- gene_hits |>
  rename(gene = g) |>
  group_by(gene) |>
  group_modify(\(d, key) tibble(
    tumors_hc_disrupted = count_tumors(d, "disrupted", TRUE),
    tumors_all_disrupted = count_tumors(d, "disrupted", FALSE),
    tumors_hc_fusion = count_tumors(d, "fusion", TRUE),
    tumors_all_fusion = count_tumors(d, "fusion", FALSE),
    tumors_hc_at_breakend = count_tumors(d, "at_breakend", TRUE),
    tumors_all_at_breakend = count_tumors(d, "at_breakend", FALSE),
    samples = join_set(d$sample_id),
    hc_events = join_set(d$event_key[d$hc])
  )) |>
  ungroup() |>
  # Only genes with 2 or more tumors in at least one count (user decision).
  filter(if_any(starts_with("tumors_"), \(x) x >= 2)) |>
  left_join(gene_role, by = "gene", relationship = "one-to-one") |>
  relocate(role, .after = gene) |>
  arrange(desc(tumors_hc_disrupted), desc(tumors_hc_fusion), desc(tumors_hc_at_breakend), gene)

fusions <- fusion_long |>
  mutate(hc = event_key %in% hc_keys) |>
  group_by(name) |>
  summarize(
    gene_5p = join_set(gene_5p), gene_3p = join_set(gene_3p),
    reported = any(reported),
    tumors_hc = n_distinct(sample_id[hc]), tumors_all = n_distinct(sample_id),
    samples = join_set(sample_id),
    hc_events = join_set(event_key[hc]), .groups = "drop"
  ) |>
  arrange(desc(tumors_hc), desc(tumors_all), name)

# ---- rule sheet: definitions and what each flag removes ------------------

flag_cols <- c("flag_own_normal", "flag_other_patients", "flag_population", "flag_linx_artifact", "flag_region")
flag_text <- c(
  flag_own_normal = "Split reads for the break in the group's own normal BAM (Prism BAM recount), either end. Never set for tumor-only groups (has_normal FALSE).",
  flag_other_patients = "Same event in another group. Two-breakend events: both breakends within 50 bases (same orientations, one SV) in another group's somatic calls (PURPLE SV VCF records other than INFERRED, so PASS and PON) or germline calls. One-breakend events: the breakend within 50 bases (same orientation) in those calls. Any event: split reads for the break in another group's normal or tumor BAM (Prism BAM recount), either end.",
  flag_population = "Population SV site matching both breakpoints within 100 bases: gnomAD-SV v4.1, 1000 Genomes, dbVar nstd166 or nstd200 (any FILTER).",
  flag_linx_artifact = "LINX puts the event in a cluster of category ARTIFACT (LINX resolved type LOW_VAF, DUP_BE, PAIR_INF or SGL_MAPPED_INF; rules in Prism/docs/SV_REPORT.md).",
  flag_region = "A breakend in GIAB low mappability, the ENCODE blacklist, the Hartwig SV blacklist or Hartwig unmappable regions."
)
only_one <- events$n_flags == 1
rule <- tibble(
  item = c(flag_cols, "high_confidence"),
  definition = c(flag_text[flag_cols], str_glue(
    "A high-confidence event is a PASS somatic SV with no sign of being an artifact: its break has no support in ",
    "the patient's normal, isn't seen in other patients or in population SV databases, isn't called an artifact by ",
    "LINX, and isn't in a hard-to-map region. In this cohort {sum(events$high_confidence)} of {nrow(events)} PASS ",
    "events qualify. Technically: no flag set. A flag is evidence, not proof; the flags are a rule chosen for this ",
    "report (2026-11-23)."
  )),
  pass_events_flagged = c(map_int(flag_cols, \(c) sum(events[[c]])), sum(events$high_confidence)),
  flagged_by_this_flag_only = c(map_int(flag_cols, \(c) sum(events[[c]] & only_one)), NA_integer_)
) |>
  add_row(item = "pass_events", definition = "Somatic SV events with FILTER PASS (PURPLE), all analyzed by LINX.",
          pass_events_flagged = nrow(events), .before = 1) |>
  add_row(item = "genes", definition = "start_gene/end_gene: LINX svs.tsv geneStart/geneEnd (genes whose span contains the breakend, or with a transcript starting 1 to 10,000 bases after it in the direction of transcription; any biotype). *_disrupted_gene: canonical transcript marked disruptive by LINX. *_driver_disruption: LINX reported driver disruption. Gene and fusion counts are tumors (sample_id). Sheet genes lists only genes with 2 or more tumors in at least one of its six counts. Details: Prism/docs/SV_REPORT.md.") |>
  add_row(item = "same_gene_fusions",
          definition = str_glue(
            "LINX fusions whose 5' and 3' breakends are in the same gene (gene_5p = gene_3p) join two parts of one gene ",
            "(intragenic rearrangements; treated as cryptic splice products by user decision), not two genes: they are ",
            "left out of fusions_real and of the fusion counts of genes (the genes still count as disrupted and at a ",
            "breakend), and kept in ",
            "the fusions columns of the events sheets. This run: {n_distinct(same_gene$sample_id, same_gene$name)} such ",
            "fusions in {n_distinct(same_gene$sample_id)} tumors."
          ),
          pass_events_flagged = n_distinct(same_gene$event_key))

# ---- write ---------------------------------------------------------------

sheets <- list(
  fusions_real = fusions,
  genes = genes,
  events_hc = events |> filter(high_confidence),
  events_all = events,
  rule = rule
)

# Glossary: one row per column name, with the sheets that have it. The
# run stops if a column has no description or a description no column.
glossary_text <- c(
  # events_hc, events_all
  group_id = "oncoanalyser group (one tumor, with its normal if any).",
  sample_id = "Tumor sample ID. Gene and fusion counts are numbers of distinct sample_id.",
  event_key = "Unique ID of the SV event: sample_id:VCF ID of the start breakend.",
  sv_type = "LINX SV type: DEL, DUP, INV, INS, BND (translocation) or SGL (single breakend: only one side mapped).",
  start_chrom = "Chromosome of the start breakend (PURPLE SV VCF CHROM).",
  start_pos = "Position of the start breakend (PURPLE SV VCF POS).",
  start_orient = "Orientation of the start breakend: 1 if the VCF ALT starts with the reference base, else -1 (as ESVEE and LINX).",
  end_chrom = "Chromosome of the end breakend; empty for single breakends (SGL).",
  end_pos = "Position of the end breakend; empty for single breakends (SGL).",
  end_orient = "Orientation of the end breakend (as start_orient); empty for single breakends.",
  length = "end_pos - start_pos for events with both breakends on one chromosome; empty otherwise.",
  start_gene = "LINX svs.tsv geneStart: genes whose span contains the start breakend, or with a transcript LINX kept starting 1 to 10,000 bases after it in the direction of transcription (any biotype); sorted, several joined with ;. Empty: LINX names no gene. Not the nearest gene: see Prism/docs/SV_REPORT.md.",
  start_disrupted_gene = "Genes whose canonical transcript LINX marks disruptive at the start breakend (breakend.tsv canonical and disruptive TRUE); sorted, several joined with ;.",
  start_driver_disruption = "Genes LINX reports as a driver disruption at the start breakend (breakend.tsv reportedStatus REPORTED: disruptive, a driver gene panel gene with disruption reporting on, its canonical or panel-listed transcript).",
  end_gene = "As start_gene, for the end breakend (LINX svs.tsv geneEnd). For a single breakend (SGL, no end position) LINX fills it from the alignments of the inserted sequence (PURPLE SV VCF INFO INSALN), not from a breakend.",
  end_disrupted_gene = "As start_disrupted_gene, for the end breakend.",
  end_driver_disruption = "As start_driver_disruption, for the end breakend.",
  fusions = "LINX fusions this event takes part in (fusion.tsv name, 5' gene_3' gene; @IGH etc. for enhancer partners), each once, including fusions of a gene with itself (left out of fusions_real; see rule sheet, same_gene_fusions); sorted, several joined with ;. A fusion made by two SVs of a chain is listed on both events.",
  reported_fusions = "Those of fusions that LINX reports (fusion.tsv reported TRUE).",
  cluster_type = "LINX resolved type of the event's cluster (clusters.tsv resolvedType), e.g. DEL, DUP, UNBAL_TRANS, LINE, COMPLEX, LOW_VAF.",
  cluster_size = "Number of SVs in the event's LINX cluster (clusters.tsv clusterCount).",
  tumor_fragments = "Tumor fragments supporting the event (PURPLE SV VCF FORMAT VF, tumor).",
  start_af = "PURPLE purity-adjusted allele fraction at the start breakend (INFO PURPLE_AF).",
  end_af = "PURPLE purity-adjusted allele fraction at the end breakend; empty for single breakends.",
  has_normal = "TRUE if the group has a matched normal; flag_own_normal can only be set when TRUE.",
  normal_split_reads = "Split reads for the break in the group's own normal BAM (Prism BAM recount), the larger of the two ends; empty for tumor-only groups and single breakends.",
  other_patients_with_event = "Other groups with the same two-breakend event (both breakends within 50 bases, same orientations, one SV) in their somatic calls (PASS and PON records) or germline calls; sorted, several joined with ;. Empty for one-breakend events, and it does not list groups found only by BAM split reads.",
  flag_own_normal = "See sheet rule, item flag_own_normal.",
  flag_other_patients = "See sheet rule, item flag_other_patients.",
  flag_population = "See sheet rule, item flag_population.",
  flag_linx_artifact = "See sheet rule, item flag_linx_artifact.",
  flag_region = "See sheet rule, item flag_region.",
  n_flags = "Number of flags set (0 to 5).",
  high_confidence = "TRUE if no flag is set (n_flags = 0). These events are the sheet events_hc.",
  # genes
  gene = "Gene name (Ensembl, as LINX writes it).",
  role = "Driver gene panel role (ONCO or TSG); empty if the gene is not on the panel.",
  tumors_hc_disrupted = "Tumors with a high-confidence event whose breakend disrupts the gene's canonical transcript (start/end_disrupted_gene).",
  tumors_all_disrupted = "As tumors_hc_disrupted, over all PASS events.",
  tumors_hc_fusion = "Tumors with a high-confidence event in a LINX fusion with the gene as 5' or 3' partner (fusion_gene_5p/3p: the gene at the fusion's breakend).",
  tumors_all_fusion = "As tumors_hc_fusion, over all PASS events.",
  tumors_hc_at_breakend = "Tumors with a high-confidence event that has the gene in start_gene or end_gene.",
  tumors_all_at_breakend = "As tumors_hc_at_breakend, over all PASS events.",
  samples = "Tumors (sample_id) behind the row's all-PASS counts (genes: gene disrupted, in a fusion or at a breakend; fusions_real: tumors with the fusion); sorted, joined with ;.",
  hc_events = "event_key of the high-confidence events behind the row's counts (genes: gene disrupted, in a fusion or at a breakend; fusions_real: events in the fusion); joined with ;.",
  # fusions
  name = "LINX fusion name: 5' gene_3' gene, or @IGH, @IGK etc. for an enhancer 5' partner.",
  gene_5p = "Gene at the fusion's 5' breakend (LINX breakend.tsv gene of fivePrimeBreakendId).",
  gene_3p = "Gene at the fusion's 3' breakend (LINX breakend.tsv gene of threePrimeBreakendId).",
  reported = "TRUE if LINX reports the fusion in any tumor (fusion.tsv reported).",
  tumors_hc = "Tumors with the fusion on a high-confidence event (a fusion of two chained SVs counts when either of its two events is high confidence).",
  tumors_all = "Tumors with the fusion on any PASS event.",
  # rule
  item = "Flag, count or term described on that row.",
  definition = "Definition of the item.",
  pass_events_flagged = "Number of PASS events with the flag set (high_confidence row: events with no flag; pass_events row: all PASS events; same_gene_fusions row: PASS events with such a fusion).",
  flagged_by_this_flag_only = "Number of PASS events with this flag and no other.",
  # glossary
  columnName = "Column name.",
  sheets = "Sheets that have the column.",
  description = "What the column holds."
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
if (length(missing) > 0) stop("glossary: no description for ", str_c(missing, collapse = ", "), call. = FALSE)
if (length(unused) > 0) stop("glossary: description for no column: ", str_c(unused, collapse = ", "), call. = FALSE)
sheets$glossary <- glossary |> mutate(description = unname(glossary_text[columnName]))

wb <- openxlsx::createWorkbook()
openxlsx::modifyBaseFont(wb, fontName = "Arial", fontSize = 10)
for (nm in names(sheets)) {
  openxlsx::addWorksheet(wb, nm)
  openxlsx::writeData(wb, nm, sheets[[nm]], withFilter = !nm %in% c("rule", "glossary"),
                      headerStyle = openxlsx::createStyle(textDecoration = "bold"))
  openxlsx::freezePane(wb, nm, firstRow = TRUE)
  openxlsx::setColWidths(wb, nm, cols = seq_along(sheets[[nm]]), widths = "auto")
}
gene_cols <- c("start_gene", "start_disrupted_gene", "start_driver_disruption",
               "end_gene", "end_disrupted_gene", "end_driver_disruption")
for (nm in c("events_hc", "events_all")) {
  openxlsx::setColWidths(wb, nm, cols = match(gene_cols, names(sheets[[nm]])), widths = 12)
}
openxlsx::setColWidths(wb, "rule", cols = 2, widths = 100)
openxlsx::addStyle(wb, "rule", openxlsx::createStyle(wrapText = TRUE, valign = "top"),
                   rows = seq_len(nrow(rule) + 1), cols = 1:4, gridExpand = TRUE)
openxlsx::setColWidths(wb, "glossary", cols = 3, widths = 100)
openxlsx::addStyle(wb, "glossary", openxlsx::createStyle(wrapText = TRUE, valign = "top"),
                   rows = seq_len(nrow(sheets$glossary) + 1), cols = 1:3, gridExpand = TRUE)
# Explicit window size: without it Excel for Mac opens the file full screen.
wb$workbook$bookViews <- str_c('<bookViews><workbookView xWindow="4000" yWindow="3000" ',
                               'windowWidth="32000" windowHeight="20000"/></bookViews>')
xlsx_file <- path(out_dir, "prism_sv_report.xlsx")
openxlsx::saveWorkbook(wb, xlsx_file, overwrite = TRUE)
write_csv(events, path(out_dir, "prism_sv_events.csv"), na = "")

message("PASS events: ", nrow(events), "; high confidence: ", sum(events$high_confidence))
print(rule |> select(item, pass_events_flagged, flagged_by_this_flag_only))
message("-> ", xlsx_file, " and ", path(out_dir, "prism_sv_events.csv"))
