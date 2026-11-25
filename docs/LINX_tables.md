# The seven LINX somatic SV tables, and why there are several

LINX reads the structural variants (SVs) in PURPLE's `purple.sv.vcf.gz`
and describes them at several levels. Each level has its own unit, and
one row of one table does not match one row of another. LINX writes one
table per unit, with ID columns that link the tables, in the way a
database is split into tables. `postOA.sh` stacks each table over all
groups into `res/prism/db/data/linx_somatic_<table>.tsv.gz`, with
`group_id` and `sample_id` in front.

What LINX reads from PURPLE, and the checks behind the statements here,
are in `ACTUAL_METHODS.md` (section 4, LINX). The gene rules below were
checked in the LINX source at tag `linx-v2.3.1`; where the LINX README
says otherwise, `ACTUAL_METHODS.md` section 5 lists it.

LINX writes eight somatic data tables per tumor: these seven and
`neoepitope.tsv` (possible neoepitopes of fusions, not stacked by
Prism). It also writes six `vis_*` tables that only feed its plots
(`vis_sv_data`, `vis_segments`, `vis_copy_number`, `vis_fusion`,
`vis_gene_exon`, `vis_protein_domain`; not stacked).

## The seven tables at a glance

| Table | One row is | Key | Links to |
|---|---|---|---|
| `svs` | one SV | `svId` | `clusters` by `clusterId` |
| `clusters` | one group of SVs that LINX treats as one event | `clusterId` | |
| `links` | one link between two breakends in a chain | | `svs` by `lowerSvId`, `upperSvId`; `clusters` by `clusterId` |
| `breakend` | one SV end in one gene transcript | `id` | `svs` by `svId` |
| `fusion` | one fusion of two transcripts (usually of two genes) | | `breakend` by `fivePrimeBreakendId`, `threePrimeBreakendId` |
| `drivers` | one driver gene event and the cluster it is assigned to | | `clusters` by `clusterId` (-1: none) |
| `driver_catalog` | one driver call for one gene | | `gene` |

All IDs are numbered within a sample. In the aggregated files, join on
`sample_id` as well as the ID.

## Why one table is not enough

Start from one SV and see how many rows it produces.

1. **`svs`: the SV itself.** An SV has two ends (breakends), for example
   a deletion or a translocation. Some have one end only: single
   breakends (type `SGL`) and breakends PURPLE inferred from a copy
   number change (type `INF`). One row per SV, with its type, the
   coordinates of both ends, its copy number, and one gene list per end
   (`geneStart`, `geneEnd`): genes whose span contains the end, or with
   a transcript starting up to 10,000 bases after it in the direction
   of transcription; any biotype. For a single breakend, `geneEnd` comes
   from the alignments of its inserted sequence. It is not the nearest
   gene, and gives no distance to it.

2. **`clusters`: SVs that belong together.** Many SVs are not separate
   events. LINX groups SVs that it judges to come from one event, for
   example a set of breaks that rearranged one chromosome region. One
   cluster holds one or many SVs. LINX also uses clusters to mark likely
   artifacts (`category` = `ARTIFACT`): single SVs it sets aside before
   clustering, mostly isolated translocations and single breakends
   whose junction copy number maximum (`svs.tsv`
   `junctionCopyNumberMax`) is below 0.75 (resolved type `LOW_VAF`);
   the rules are in `ACTUAL_METHODS.md` and `SV_REPORT.md` section 6.4.
   Cluster facts (category, resolved type, number of SVs) belong to the
   cluster, not to each SV, so they are in their own table.

3. **`links`: the order of the SVs in a cluster.** In a complex cluster
   LINX tries to connect the SV ends into chains: the end of one SV is
   joined to the end of the next by a stretch of DNA. One row per such
   connection. It joins two SVs, so it cannot be a column of `svs`.
   Simple SVs have no links.

4. **`breakend`: what each end hits.** LINX attaches to each end every
   gene whose span contains it or whose 5' end is up to 100,000 bases
   after it (strand-aware), of any biotype; 250,000 bases near the 3'
   genes of known IG-enhancer pairs, and a few configured genes also on
   their 3' side (region type DOWNSTREAM).
   For each such gene it writes the canonical transcript, plus, for
   driver panel genes with disruption reporting, the panel's additional
   transcripts, plus any transcript used in a fusion (64,474 of the
   64,996 cohort rows are canonical). A transcript the end lies past the
   3' end of gets no row. Each row says where the end lies
   (`regionType`: EXONIC, INTRONIC, UPSTREAM, DOWNSTREAM for a few
   configured genes, or IG for the enhancer partner of an IG fusion),
   whether LINX calls it disruptive and whether it is reported. Most
   rows are UPSTREAM. An end gets no row when no gene is attached to
   it, and also, by default, the end of a single breakend (whose genes
   come from the alignments of its inserted sequence) unless it is in a
   fusion, and single and inferred breakends LINX linked to each other.

5. **`fusion`: transcripts joined.** If the two sides of a rearrangement
   join the start of one transcript to the end of another in a way that
   could make a fused transcript, LINX lists a fusion. Usually the two
   are different genes; in this cohort 96 of the 1,171 fusions join two
   parts of one gene (an intragenic rearrangement, named `GENE_GENE`).
   A fusion is made from two breakend rows (5' and 3'), which may come
   from one SV or from two SVs of a chain. One row per fusion;
   `reported` = true for the ones LINX reports.

6. **`driver_catalog`: gene-level driver calls.** One row per gene with
   a driver call: the reported rows of PURPLE's somatic driver catalog
   (mutations, amplifications, deletions), copied unchanged, plus the
   gene disruptions LINX adds (`driver` = `DISRUPTION`). This is a list
   of genes, not of SVs.

7. **`drivers`: which cluster caused a driver.** For a driver gene
   whose event LINX can explain by an SV cluster (for example an
   amplification made by a cluster), one row links the gene to that
   cluster. `clusterId` = -1 where there is no cluster, for example
   arm-level gains.

The units differ: an SV, a group of SVs, a pair of SV ends, an SV end in
a transcript, a pair of transcripts (a fusion), a gene, and a gene-cluster pair. Putting
them in one table would repeat each cluster and SV on many rows and
leave most cells empty.

## How the tables link

```
clusters ──< svs ──< breakend ──< fusion (5' and 3' breakend)
   │  │       │
   │  │       └──< links (lowerSvId, upperSvId)
   │  └──< links (clusterId)
   └──< drivers >── driver_catalog (by gene)
```

`──<` means one row on the left matches many on the right.

## Which table answers which question

| Question | Table |
|---|---|
| How many SVs, of which type, per sample | `svs` (decide first what to do with `ARTIFACT` clusters and `INF`) |
| Is an SV part of a complex event | `svs` joined to `clusters` |
| Which genes an SV breaks | `breakend` (`disruptive` = true, `canonical` = true) |
| Which samples have a disruption of gene X | `breakend`, or `driver_catalog` for reported driver genes |
| Which fusions | `fusion` (`reported` = true for those LINX reports) |
| The driver genes of a sample | `driver_catalog` |
| Which cluster made a driver | `drivers` joined to `clusters` |

## Examples in R

```r
library(tidyverse)
db <- "res/prism/db/data"
rd <- function(t) read_tsv(file.path(db, str_c("linx_somatic_", t, ".tsv.gz")),
                           col_types = cols(.default = "c"))
svs <- rd("svs")
clusters <- rd("clusters")
breakend <- rd("breakend")
fusion <- rd("fusion")

# SVs per sample by type, without SVs in ARTIFACT clusters
svs |>
  left_join(clusters, by = c("group_id", "sample_id", "clusterId")) |>
  filter(category != "ARTIFACT") |>
  count(sample_id, type)

# Genes disrupted on their canonical transcript, with the SV behind each
breakend |>
  filter(canonical == "true", disruptive == "true") |>
  left_join(svs, by = c("group_id", "sample_id", "svId")) |>
  distinct(sample_id, gene, svId, type)

# Reported fusions with the two genes
fusion |>
  filter(reported == "true") |>
  left_join(breakend |> select(group_id, sample_id, id, gene_5p = gene),
            by = c("group_id", "sample_id", fivePrimeBreakendId = "id")) |>
  left_join(breakend |> select(group_id, sample_id, id, gene_3p = gene),
            by = c("group_id", "sample_id", threePrimeBreakendId = "id")) |>
  select(sample_id, name, gene_5p, gene_3p, reportedType)
```

## Points to watch

- `svs` holds the SVs of `purple.sv.vcf.gz` with FILTER PASS or
  INFERRED. SVs with FILTER PON are not in any LINX table.
- SVs in `ARTIFACT` clusters stay in all tables. Decide for each count
  whether to include them.
- Type `INF` rows are breakends PURPLE inferred from copy number, not
  SVs seen in the reads.
- For some SVs whose two ends have the same position, `vcfIdEnd`
  repeats `vcfIdStart`. To join back to the VCF, use `vcfIdStart`.
- `breakend` has the canonical transcript of every attached gene plus
  a few others (driver panel, fusions); filter on `canonical` = true to
  count genes once. Most rows are UPSTREAM: the end is before the gene,
  not in it. Filter on `regionType` EXONIC or INTRONIC for genes the end
  falls in.
- No LINX table gives the nearest gene to an end or the distance to a
  gene's 5' end (`breakend` has only `nextSpliceDistance`, to the next
  splice site), and genes on the 3' side of an end are not attached,
  except a few configured genes (`SV_REPORT.md` sections 5.1 and 5.4).

## The germline tables

Groups with a matched normal also have six germline tables (`svs`,
`breakend`, `disruption`, `clusters`, `links`, `driver_catalog`), from
`purple.sv.germline.vcf.gz`. They follow the same plan. There is no
fusion or drivers table. `disruption` lists the germline SVs in genes,
with each SV's positions, FILTER and read support on the same row.
