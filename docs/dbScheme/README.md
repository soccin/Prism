# dbScheme: which file holds each per-sample value

oncoanalyser runs each sample (tumor, with or without a matched normal)
on its own and has no multi-sample step. A planned post-processing step
will collect the per-sample results into cohort tables. It should read
each value once, from the file that first holds it (the primary file),
and skip files that only repackage other files. This folder records
which file that is.

Tool names (PURPLE, AMBER, ...) are explained in `../OA_AlgoKey.md`.
Paths are relative to `out/<group>/`: `<T>` tumor ID, `<N>` normal ID,
`<S>` either one (one file per BAM).

## Summary

The per-sample outputs fall into seven groups:

| Group | What it holds | Main primary files | Status |
|---|---|---|---|
| A | One value set per sample: purity, ploidy, sex, MSI, TMB, QC, HRD, signatures, tissue of origin, HLA type, pharmacogenomics, viruses, telomere length, IG/TCR loci | `purple.purity.tsv`, `purple.qc`, `amber.qc`, bam-tools summary, `chord.prediction.tsv`, `sig.allocation.tsv`, `cuppa.vis_data.tsv`, `lilac*.tsv`, `peach*.tsv`, `virusbreakend.vcf.summary.tsv`, `teal.tellength.tsv`, `cider.locus_stats.tsv` | Described field by field (`group_a.md`) |
| B | Somatic events, many rows per sample: SNVs and indels, SVs, copy number, fusions, disruptions, drivers | `purple.somatic.vcf.gz`, `purple.sv.vcf.gz`, `purple.cnv.somatic.tsv`, `purple.cnv.gene.tsv`, `purple.chromosome_arm.tsv`, `purple.driver.catalog.somatic.tsv`, `linx.*.tsv` | Files listed; fields not yet described |
| C | Germline calls (matched samples only) | `purple.germline.vcf.gz`, `purple.germline_amp_del.tsv`, `purple.driver.catalog.germline.tsv`, `linx.germline.*.tsv` | Files listed; fields not yet described |
| D | QC distributions, for plots | `bam_metric.coverage.tsv` and other bam-tools tables, `redux.*.tsv`, `cobalt.gc.median.tsv`, `purple.purity.range.tsv` | Not described |
| E | Intermediate files read by later pipeline steps | SAGE, PAVE and ESVEE VCFs, AMBER and COBALT per-site files | Not for aggregation |
| F | Display files | ORANGE PDF, LINX HTML, QSEE PDF, PNGs | Not for aggregation |
| G | Multi-sample outputs that exist now | MultiQC report, QSEE cohort run (`qcCohort.sh`), cohort report (`cohortReport.sh`) | Outside this scheme |

The full file list with formats and row counts is in the per-sample
readout catalog that started this work; this folder takes it group by
group.

## What is known so far

**Group A.** 18 primary files hold the data. Five files hold no primary
data:

| Secondary file | What it repackages |
|---|---|
| ORANGE JSON | Every Group A value, plus the reported records of Groups B and C |
| `qsee.status.tsv.gz` | bam-tools, AMBER and PURPLE values; adds only PASS/WARN/FAIL |
| `qsee.vis.data.tsv.gz` | bam-tools, COBALT, REDUX and ESVEE distributions |
| `cuppa.pred_summ.tsv` | The top 3 per classifier of `cuppa.vis_data.tsv` |
| `virus.annotated.tsv` | VIRUSBreakend values; adds names and the reporting decision |

Detail per tool: `group_a.md`.

**ORANGE.** ORANGE is a per-sample summary, not a primary source:
- its event lists are exactly the records the tool files already flag
  as reported (checked in all 14 samples);
- it adds only a few values derived by fixed rules (immune escape
  flags, arm gain/loss calls);
- it drops UGT1A1 and has several mislabeled fields.

Detail: `orange.md`.

**Groups B and C.** Not yet worked through. Known so far:
- The ORANGE map (`tables/orange.csv`) names the source files of the
  reported events.
- The two somatic driver catalogs overlap. The LINX catalog holds
  PURPLE's reported entries plus LINX's disruption entries; PURPLE's
  catalog also holds non-reported candidates.
- Small variants pass through SAGE, then PAVE, then PURPLE. The PURPLE
  VCF is expected to hold the earlier calls plus annotations, but that
  is not yet checked.

## Roles

Each field gets one role:

| Role | Meaning |
|---|---|
| primary | First written here, by the tool that computes it |
| copy | The same value as a field of another file (renamed, rounded or recoded) |
| derived | Computed from values another file holds, by a stated rule |
| annotation | Looked up from a fixed reference table (drug links, names, thresholds); not measured |
| key | Identifies the row (sample, gene, classifier); holds no result |

A file is **primary** if it has at least one primary field, and
**secondary** if it has none. Where two files of the same tool hold the
same value, the tool's main result file counts as primary.

## Files here

| File | Contents |
|---|---|
| `group_a.md` | Group A per tool: primary file, unique fields, copies |
| `orange.md` | ORANGE: where each part comes from, what it adds, errors found |
| `dbScheme.xlsx` | The three tables below, with a `readme` sheet that defines every column |
| `tables/files.csv` | One row per file: tool, role, what it holds, field counts by role, presence |
| `tables/fields.csv` | Every field of the described files: role, description, source of copies and derived values, copy check, presence |
| `tables/orange.csv` | The ORANGE JSON by section: source file, records kept, values ORANGE computes |
| `curated/files.csv`, `curated/fields.csv`, `curated/orange_map.csv` | Hand-written inputs to the build |
| `scripts/fetchOrangeDatamodel.R` | Downloads the ORANGE data model for a release tag |
| `scripts/buildScheme.R` | Builds `tables/` and `dbScheme.xlsx` |

`tables/` and `dbScheme.xlsx` are generated; edit `curated/` and
rebuild.

## Rebuild

From the project folder (the parent of `Prism/`):

```
Rscript Prism/docs/dbScheme/scripts/fetchOrangeDatamodel.R orange-v5.0.1 DATAMODEL_DIR
Rscript Prism/docs/dbScheme/scripts/buildScheme.R --datamodel DATAMODEL_DIR --out-dir out
```

`DATAMODEL_DIR` is any scratch folder. The build reads the tool files
and the ORANGE JSON of every sample in `out/`. It:

- marks which files and fields are present in matched and tumor-only
  samples;
- compares every copy in a one-value-per-sample file with its source,
  sample by sample, and writes the result to `check` (for example
  `equal in 23 of 23`);
- assigns every ORANGE JSON path to a row of `orange_map.csv`.

It writes only all/some/none and counts: no sample names and no values.
It warns on each mismatch:
- a field missing from `fields.csv`;
- a copy that differs from its source;
- a file role that does not match its field roles;
- an ORANGE path not covered by the map.

Copies in multi-row files (Virus Interpreter) and copies of COBALT
values were checked once with separate scripts; the `notes` column says
so. On a pipeline upgrade, fetch the new ORANGE tag, rebuild, work
through the warnings, and read ORANGE's source code again: what it
computes is not in its data model.

## Sources

Built for oncoanalyser 3.0.0. Tool versions (from the pipeline's
container tags) and the documents used:

| Tool | Version | Docs |
|---|---|---|
| ORANGE | 5.0.1 | `orange/README.md`, `orange-datamodel/` and `orange/src/main/java/` at tag `orange-v5.0.1` |
| PURPLE | 4.4 | `purple/README.md` at `purple-v4.4` |
| AMBER | 4.3 | `amber/README.md` at `amber-v4.3` |
| CHORD | 2.1.2 | `chord/README.md` at `chord-v2.1.2` |
| CUPPA | 2.5.1 | `cuppa/README.md` at `cuppa-v2.5.1` |
| LILAC | 2.0 | `lilac/README.md` at `lilac-v2.0` |
| PEACH | 2.0.0 | `peach/README.md` at `peach-v2.0.0` |
| Virus Interpreter | 1.7.2 | `virus-interpreter/README.md` at `virus-interpreter-v1.7.2` |
| SIGS | 1.2.1 | `sigs/README.md` at `sigs-v1.2.1` (no output description) |
| TEAL | 1.4 | `teal/README.md` at `teal-v1.4` |
| CIDER | 1.2 | `cider/README.md` at `cider-v1.2` |
| bam-tools | 1.6.1 | `bam-tools/README.md` at `bam-tools-v1.6.1` (no column descriptions) |
| QSEE | 1.0 | `qsee/README.md` at `qsee-v1.0` |
| VIRUSBreakend | GRIDSS 2.13.2 | `VIRUSBreakend_Readme.md` in https://github.com/PapenfussLab/gridss at `v2.13.2` |

The hmftools tools are in https://github.com/hartwigmedical/hmftools.
Copies and derived values were found by comparing values across files
in all 14 samples and, for ORANGE, by reading its source code.
