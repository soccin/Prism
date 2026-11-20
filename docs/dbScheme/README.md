# dbScheme: what the pipeline writes per sample

oncoanalyser runs each sample (tumor, with or without a matched normal)
on its own and has no multi-sample step. A planned post-processing step
will collect the per-sample results into cohort tables. This folder
documents the per-sample results field by field, as the basis for those
tables.

Work so far covers Group A (one value set per sample) and ORANGE.

## Groups

The per-sample outputs fall into these groups:

| Group | Contents |
|---|---|
| A | One value set per sample: purity and ploidy, QC, MSI, TMB, HRD, signatures, tissue of origin, telomere length, HLA type, pharmacogenomics, viruses, immune receptor loci |
| B | Somatic events, many rows per sample: SNVs and indels, SVs, fusions, copy number, drivers |
| C | Germline calls (matched samples only) |
| D | QC distributions (coverage, fragment length, duplicates), for plots |
| E | Intermediate files read by later pipeline steps |
| F | Display files (PDF, HTML, PNG) |

## Files

| File | Contents |
|---|---|
| `orange.md` | The ORANGE JSON: structure, what it holds, what it lacks |
| `group_a.md` | Group A readouts, and which of their fields ORANGE holds |
| `dbScheme.xlsx` | The tables below in one workbook, with a `readme` sheet that defines every column |
| `tables/orange_fields.csv` | Every ORANGE JSON path: type, group, level, filled in matched and tumor-only samples, description, source file and column |
| `tables/orange_enums.csv` | Allowed values of the ORANGE enum fields |
| `tables/group_a_fields.csv` | Every field of every Group A tool file, its description, its ORANGE path, and whether it is present in matched and tumor-only samples |
| `curated/orange_descriptions.csv` | Hand-written descriptions of ORANGE fields, keyed by data model type and JSON key (input to the build) |
| `curated/group_a_fields.csv` | Hand-written Group A field list (input to the build) |
| `scripts/fetchOrangeDatamodel.R` | Downloads the ORANGE data model for a release tag |
| `scripts/buildScheme.R` | Builds `tables/` and `dbScheme.xlsx` |

`tables/` and `dbScheme.xlsx` are generated; edit `curated/` and
rebuild.

File names use placeholders: `<T>` tumor ID, `<N>` normal ID, `<S>`
either one (tools that write one file per BAM). Paths are relative to
`out/<group>/`.

## Rebuild

From the project folder (the parent of `Prism/`):

```
Rscript Prism/docs/dbScheme/scripts/fetchOrangeDatamodel.R orange-v5.0.1 DATAMODEL_DIR
Rscript Prism/docs/dbScheme/scripts/buildScheme.R --datamodel DATAMODEL_DIR --out-dir out
```

`DATAMODEL_DIR` is any scratch folder. The build reads the ORANGE JSON
and the Group A files of every sample in `out/` to check the curated
tables against real output and to mark which fields are filled in
matched and tumor-only samples. It writes only all, some or none per
field: no sample names and no values. It prints a warning for each
mismatch: a file field missing from `curated/group_a_fields.csv`, an
ORANGE path that does not exist, a Group A ORANGE field without a
description, a JSON key missing from the data model.

After a pipeline upgrade, look up the new ORANGE version in the
pipeline's `modules/local/orange/main.nf`, fetch that tag, rebuild, and
work through the warnings.

## Versions

Built for oncoanalyser 3.0.0. Tool versions (from the pipeline's
container tags) and the README each description comes from:

| Tool | Version | Docs |
|---|---|---|
| ORANGE | 5.0.1 | `orange/README.md`, `orange-datamodel/` at tag `orange-v5.0.1` |
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

All in https://github.com/hartwigmedical/hmftools. Where a README does
not describe a field, or describes it wrongly, the `notes` column says
so.
