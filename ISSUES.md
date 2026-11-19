# Issues

Open work on Prism. Add new items at the end.

## 1. Revisit the gene summary table in prism_report

**Status:** open. The tables were taken out of the report.

`cohortReport.sh` used to write three per-gene count tables to
`prism_report.xlsx`: `gene_summary_HC`, `gene_summary_UnFilt` and
`gene_summary_by_class`, plus `data/gene_summary_all.csv`. They gave,
for each gene, the number of samples with each kind of event. They
were removed because they were hard to read and need more work. The
report now has no per-gene count over the cohort.

The removed code is in commit `bfac29b`, in
`scripts/cohortReport.R`: the section "gene summary"
(`make_events`, `make_gene_summary`), the `-p FRAC` option
(`--min-pct`) and the matching entries of `manifest_text` and
`glossary_text`.

What made the tables hard to read:

- 32 columns, with several counts of the same thing: `n_disruption`
  and `n_disruption_incl_weak`, `n_fusion` and
  `n_fusion_incl_unreported`, six breakend and breakpoint counts.
- Three tables for one question. HC and UnFilt differed only in
  whether a breakend at a breakpoint shared with another sample
  counts toward a disruption.
- `n_any` mixed event types of unequal weight (one missense mutation
  counts the same as a deep deletion) and left out LOH.
- Samples without a matched normal were counted with the others,
  although their mutation calls can include germline variants.
- One large copy number event puts every gene it covers into the
  table, so neighboring genes fill the top rows.
- A deep deletion that the copy number caller derives from a shared
  breakend was still counted in the HC table.
- The cutoff max(2, FRAC x samples) is 2 for any cohort under 60
  samples, so the table was long (thousands of genes).

To decide before bringing a table back:

- The question it answers, and for whom.
- Which events count, and whether each event type gets its own
  ranking instead of one combined count.
- How samples without a matched normal are counted.
- Whether copy number events are counted per gene or per region.
- One table, not three.

The table of the genes of interest with one row per gene and sample
(`genes_of_interest`) was removed at the same time;
`genes_of_interest_events` lists the events of those genes.
