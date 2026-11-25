# ORANGE JSON: a secondary file

ORANGE is the last step of oncoanalyser. For each sample it reads the
output files of the other tools and writes a PDF report and a JSON
file, `orange/<T>.orange.json`. It analyzes no sequencing data. The
JSON holds no primary data: every value is a copy of a tool file
value, a value derived from those by a fixed rule, a look-up, a run
setting or a plot path.

This page is based on the ORANGE 5.0.1 source code (`orange/` and
`orange-datamodel/` at tag `orange-v5.0.1`) and on the JSON files of all
14 samples. Detail is in `tables/orange.csv` (sheet `orange` of
`dbScheme.xlsx`). It has one row per JSON section, giving the source
file and which records ORANGE keeps, and one row per value ORANGE
computes.

## Where each section comes from

| JSON section | Source file | Records kept |
|---|---|---|
| `purple.fit` | `purple.purity.tsv`, `purple.qc` | |
| `purple.characteristics` | `purple.purity.tsv`, `purple.qc` | |
| `purple.somaticDrivers` | `purple.driver.catalog.somatic.tsv` | entries not NOT_REPORTED |
| `purple.germlineDrivers` | `purple.driver.catalog.germline.tsv` | entries not NOT_REPORTED |
| `purple.somaticVariants` | `purple.somatic.vcf.gz` | PASS records with the REPORTED flag |
| `purple.germlineVariants` | `purple.germline.vcf.gz` | records with the REPORTED flag |
| `purple.somaticGainsDels` | somatic driver catalog + `purple.cnv.gene.tsv` + `purple.chromosome_arm.tsv` | AMP, PARTIAL_AMP, DEL, HET_DEL, LOH driver entries |
| `purple.germlineGainsDels` | germline driver catalog + `purple.germline_amp_del.tsv` | REPORTED GERMLINE_DELETION and GERMLINE_AMP entries |
| `purple.armCopyNumberAbberations` | `purple.chromosome_arm.tsv` + ploidy and sex | arms called GAIN or LOSS by ORANGE |
| `linx.fusions` | `linx.fusion.tsv` + `linx.breakend.tsv` + `linx.svs.tsv` | reported = true |
| `linx.somaticBreakends` | `linx.breakend.tsv` + `linx.svs.tsv` + `linx.driver.catalog.tsv` | reportedStatus = REPORTED |
| `linx.germlineBreakends` | the germline LINX files | reportedStatus = REPORTED |
| `lilac` | `lilac.tsv`, status from `lilac.qc.tsv` | all alleles |
| `virusInterpreter` | `virus.annotated.tsv` | not blacklisted; `reportableViruses`: also reported = true |
| `chord` | `chord.prediction.tsv` | |
| `cuppa` | `cuppa.vis_data.tsv` (not `pred_summ.tsv`) | probability rows, 4 SV feature values |
| `peach` | `peach.haplotypes.best.tsv` | genes other than UGT1A1 |
| `sigAllocations` | `sig.allocation.tsv` | all |

The selection of records is the tools' own: ORANGE keeps the records
the tool files already flag as reported. In all 14 samples, the number
of somatic variants, somatic drivers, fusions, somatic breakends and
viruses in the JSON equals the number of flagged records in the tool
files. The somatic variant list holds 1 to 30 variants per sample,
under 0.1% of the PASS records in the VCF.

## What ORANGE adds

Nothing measured. These values are in no tool file, but each follows
from tool file values by a fixed rule or a look-up table:

- **Immune escape** (`immuneEscape`, 6 flags). True if a gene of a fixed
  list is hit:
  - deleted (a DEL entry in `somaticGainsDels`);
  - homozygously disrupted (a somatic breakend with undisrupted copy
    number below 0.5);
  - amplified, for the PD-L1 and SETDB1 flags (an AMP entry);
  - for HLA, also any LOH, DEL or HET_DEL driver.

  Small variants are not tested.
- **Arm gain and loss** (`purple.armCopyNumberAbberations`): arms whose
  mean copy number is above or below a threshold of the expected
  ploidy, with type, relative copy number and a driver interpretation.
- **Driver interpretation** on drivers and gains/deletions: HIGH
  (likelihood 0.8 or more), MEDIUM (0.2 or more), LOW.
- **Gains/deletions:** `geneRange` (FULL or PARTIAL). Germline gains and
  deletions: the minimum and maximum copy number over the gene's rows.
- **Fusions:** `contextUp` and `contextDown` (text such as "Exon 3"),
  `junctionCopyNumber` (mean of the SV junction copy number estimates).
- **Breakends:** `junctionCopyNumber` (mean of min and max) and
  `chromosomeBand` (cytoband look-up).
- **Variants:** `localPhaseSets` (phase set IDs shared with another
  listed variant within 1000 bases), `otherImpacts`, `reported`.
- **Signatures:** `etiology` (look-up table).
- **PEACH:** `haplotype` (haplotype with zygosity, such as `*1_HOM`).
- **CUPPA:** `mode` (WGS or WGTS).
- **Run settings:** `samplingDate` (the run date), `experimentType`,
  `refGenomeVersion`, `pipelineVersion`.
- **Plot file paths.**

## What ORANGE leaves out

- **Records not flagged as reported:** over 99.9% of the somatic
  variants, NOT_REPORTED driver candidates, non-reported fusions and
  breakends, blacklisted viruses.
- **PEACH UGT1A1:** ORANGE's code removes it.
- **Whole files:** bam-tools metrics, TEAL, CIDER, the QSEE table, PEACH
  QC, AMBER QC except contamination, the VIRUSBreakend detection
  columns.
- **Single fields:**
  - PURPLE `gender`, `AmberGender` and `CobaltGender`: ORANGE has no
    sex call.
  - CHORD remarks.
  - LILAC unique, shared and wildcard fragment counts.
  - Virus `taxid`.
  - CUPPA feature contributions and signature quantiles.

## Errors in ORANGE 5.0.1

Found in the source code and confirmed in the JSON files:

1. `cuppa.simpleDups32To200B` holds the CUPPA feature
   `sv.SIMPLE_DEL_20KB_1MB`: the number of simple deletions of 20 kb to
   1 Mb, not duplications of 32 to 200 bases.
2. The antigen presentation escape list has `TABBP`, which is not a
   gene symbol; `TAPBP` is never tested.
3. `purple.somaticGainsDels[].geneRange` is FULL only for AMP; every
   other type, including a whole-gene deletion, is labeled PARTIAL.
4. In `otherImpacts`, the codon, exon, splice region and coding effect
   are copied from the canonical transcript, not taken from the other
   transcript. The list is empty in every sample here.
5. `purple.armCopyNumberAbberations` selects arms by mean copy number
   but reports the median in `copyNumber`. The key is misspelled.
6. `pipelineVersion` is ORANGE's default value (`OA_V3_0`); ORANGE reads
   the pipeline version file but does not use it.

Misleading names: `samplingDate` holds the date ORANGE ran;
`purple.characteristics.lohPercentage` holds
a fraction (0 to 1), rounded to 3 significant digits.

## For post-processing

Read the tool files listed above, not the JSON. ORANGE is needed only
for the immune escape flags, and those can be recomputed from the
PURPLE driver catalog, gene copy numbers and LINX breakends using the
rules above. The tool files give the selection back through their own
flags (REPORTED in the VCF, `reportedStatus`, `reported`).
`scripts/orangeSvReport.R` collects the JSON's `linx` records of all
samples and checks on every run that they are exactly those flagged
LINX rows (`SV_REPORT.md` section 13). ORANGE tests no SV FILTER value.

Tumor-only samples have `referenceId`, `chord`, `cuppa`, `peach`,
`sigAllocations` and all germline sections null, because the tools
that fill them need a matched normal.
