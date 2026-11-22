# oncoanalyser tool name key

What the Hartwig (HMF) tool names in oncoanalyser 3.0.0 mean. Output
folders are under `out/<group>/`.

| Name | Plain name | What it does | Output |
|---|---|---|---|
| REDUX | BAM cleanup | Marks duplicates and merges them into consensus reads; unmaps reads in known problem regions. Its BAM is the one all later steps use | `alignments/` |
| AMBER | Allele frequencies at SNPs | Tumor B-allele frequency at heterozygous SNPs (input to copy number); contamination; homozygous regions of the normal | `amber/` |
| COBALT | Read depth ratios | GC-corrected tumor/normal read depth ratios along the genome (input to copy number) | `cobalt/` |
| SAGE | Small variant caller | SNVs, MNVs and short indels, somatic and germline | `sage/` |
| PAVE | Variant effect annotation | Gene, transcript and protein effect of each variant | `pave/` |
| ESVEE | Structural variant caller | Deletions, duplications, inversions, translocations, single breakends | `esvee/` |
| PURPLE | Purity, ploidy and copy number | Purity and ploidy fit, allele-specific copy number, final variant calls, driver list, TMB, MSI, WGD, sex | `purple/` |
| LINX | Structural variant interpretation | Groups SVs into events, finds gene disruptions and fusions | `linx/` |
| CHORD | HR deficiency | Probability of homologous recombination deficiency | `chord/` |
| SIGS | Mutational signatures | Fits SNV counts to COSMIC signatures | `sigs/` |
| CUPPA | Tissue of origin | Probability of each of 40 cancer types | `cuppa/` |
| LILAC | HLA typing | HLA alleles, their copy number and somatic mutations | `lilac/` |
| NEO | Neoepitope prediction | Mutant peptides scored for presentation by the patient's HLA alleles | `neo/` |
| CIDER | IG and TCR rearrangements | CDR3 sequences and read counts per IG/TCR locus | `cider/` |
| TEAL | Telomere length | Telomere length per BAM; telomeric rearrangements | `teal/` |
| VIRUSBreakend | Virus detection | Viral reads, coverage and integration sites (part of GRIDSS, not HMF) | `virusbreakend/` |
| Virus Interpreter | Virus annotation | Names viruses, groups them (HPV, EBV, ...), decides which to report | `virusinterpreter/` |
| PEACH | Pharmacogenomics | DPYD and UGT1A1 haplotypes with linked drugs | `peach/` |
| ORANGE | Per-sample summary report | PDF and JSON of the main results of all tools; holds no primary data (`dbScheme/orange.md`) | `orange/` |
| QSEE | QC report | QC metrics against PASS/WARN/FAIL thresholds | `qsee/` |

Need a matched normal: CHORD, SIGS, CUPPA, NEO, PEACH, and the germline
runs of SAGE and LINX.

## What the names stand for

From each tool's README. The others (AMBER, SAGE, ESVEE, LINX, LILAC,
TEAL, CIDER, QSEE) give no expansion.

| Name | Stands for |
|---|---|
| PURPLE | PURity PLoidy Estimator |
| COBALT | COunt BAm Lines |
| REDUX | REcalibrate, DEduplicate, Unmap, eXtract |
| PAVE | Prediction and Annotation of Variant Effects |
| CHORD | Classifier of HOmologous Recombination Deficiency |
| CUPPA | Cancer of Unknown Primary Prediction Algorithm |
| PEACH | Pharmacogenomic Evaluator And Caller of Haplotypes |
| ORANGE | Oncogenic Results of Analyzing the Genome |
| VIRUSBreakend | Viral Integration Recognition Using Single Breakends |
| SIGS | signatures |
| NEO | neoepitopes |
