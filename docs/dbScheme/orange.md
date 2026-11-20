# ORANGE JSON

ORANGE is the last step of oncoanalyser. For each sample it writes a PDF
report and a JSON file, `orange/<T>.orange.json`, that combine the
results of the other tools. This page describes the JSON as written by
ORANGE 5.0.1 (oncoanalyser 3.0.0).

Field-level detail is in `tables/orange_fields.csv` (sheet
`orange_fields` of `dbScheme.xlsx`); enum values are in
`tables/orange_enums.csv`.

## Where the structure comes from

The JSON has no published field documentation. Its structure is defined
by the Java classes in `orange-datamodel/` of the hmftools repo, and
the ORANGE README sends readers to each tool's own documentation. The
field table is therefore built in two parts:

- Structure (paths, types, nullable, enums) parsed from the data model
  at tag `orange-v5.0.1` and checked against real JSON files. Every key
  found in the files is in the data model.
- Descriptions and source columns written by hand from the tool
  READMEs, by matching each ORANGE field to the tool file column with
  the same value.

ORANGE 5.0.1 has 447 JSON paths: 121 in Group A, 145 somatic event
(Group B), 115 germline (Group C), 54 RNA, 12 plot paths. Descriptions
exist for Group A and for the top-level containers; Groups B and C are
not yet described.

## Top-level structure

| Key | Group | Holds | Tumor-only |
|---|---|---|---|
| `sampleId`, `referenceId` | A | Tumor and normal sample IDs | `referenceId` is null; ORANGE uses that to mark tumor-only mode |
| `samplingDate`, `experimentType`, `refGenomeVersion`, `pipelineVersion`, `configuredPrimaryTumor` | A | Run metadata. `samplingDate` is the run date; `configuredPrimaryTumor` is empty | same |
| `purple.fit` | A | Purity, ploidy, fit method, PURPLE QC (status, contamination, segments, deleted genes, AMBER depth, germline aberrations) | filled |
| `purple.characteristics` | A | WGD, MSI, TMB, TML, SV burden, LOH proportion | filled |
| `purple.somaticDrivers`, `somaticVariants`, `somaticGainsDels`, `armCopyNumberAbberations` | B | Somatic drivers, SNVs and indels, gene gains and deletions, arm changes | filled |
| `purple.germlineDrivers`, `germlineVariants`, `germlineGainsDels` | C | Germline drivers, variants, gene amplifications and deletions | null |
| `linx.fusions`, `linx.somaticBreakends` | B | Fusions and somatic breakends in genes | filled |
| `linx.germlineBreakends` | C | Germline breakends in genes | null |
| `lilac.alleles[]` | A | HLA alleles (class I and II) with fragment counts, copy number, somatic mutation counts, QC status | filled; `refFragments` null |
| `immuneEscape` | A | Six immune escape flags computed by ORANGE | filled |
| `virusInterpreter` | A | All detected viruses, and the reported ones | filled |
| `chord` | A | HRD probabilities and status | null |
| `cuppa` | A | Tissue of origin: every cancer type with each classifier's probability, the best call, four SV feature counts | null |
| `peach[]` | A | Pharmacogenomic haplotypes | null |
| `sigAllocations[]` | A | Mutational signature allocations with etiology | null |
| `isofox` | RNA | RNA results | null (no RNA data) |
| `plots` | display | Paths of the plot images in the PDF | filled |

## Facts that matter for aggregation

1. **The event lists are selections.** `purple.somaticVariants`,
   `somaticGainsDels`, `linx.somaticBreakends` and `linx.fusions` hold
   far fewer entries than the tool files: `somaticVariants` held under
   0.1% of the passing records of `purple.somatic.vcf.gz` in every
   sample checked, even though the ORANGE README says the JSON holds
   every mutation. ORANGE cannot be the source of the Group B event
   tables.
2. **ORANGE does not hold all of Group A.** See `group_a.md`: 75 of the
   202 Group A file fields have an ORANGE equivalent. BAM metrics,
   telomere length, immune receptor loci, the QSEE status table and the
   PURPLE sex call are missing.
3. **Tumor-only samples** have `referenceId`, `chord`, `cuppa`, `peach`,
   `sigAllocations` and all germline lists null. The tools that would
   fill them need a matched normal.
4. **PEACH is cut down.** ORANGE keeps only DPYD; PEACH also calls
   UGT1A1 (in all matched samples checked).
5. **Some JSON keys differ from the data model names**, set by
   `@SerializedName`: `chord.BRCA1Value`, `chord.BRCA2Value`,
   `cuppa.LINECount`. The key `purple.armCopyNumberAbberations` is
   misspelled in ORANGE.
6. **Misleading names.** `purple.characteristics.lohPercentage` and
   `sigAllocations[].percent` hold fractions (0 to 1).
7. **Computed field in the JSON.** `purple.somaticVariants[].reported`
   is computed by the data model (true if any transcript impact is
   reported) and is written to the JSON.
8. **No BAM metrics.** The data model still has classes for BAM metrics
   and flagstat (`OrangeSample`, `WGSMetrics`, `Flagstat`), but the
   ORANGE record no longer refers to them, so the JSON has none.
9. **ORANGE adds a few values not in the Group A files:** the immune
   escape flags, signature etiology, the full CUPPA prediction list
   (the tool's summary file keeps the top 3 per classifier), and four
   CUPPA SV feature counts taken from `cuppa_data.tsv.gz`.
