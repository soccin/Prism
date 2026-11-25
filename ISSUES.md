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

## 2. Check resource settings on the next full run

**Status:** open. Check after the next run from FASTQ.

The 2026-11-18 `conf/iris.config` has so far run only in a resume,
which ran ESVEE and the steps after it for a few groups. Results are
in `docs/NF_SCRIPTING_INFO.md` sections 14 and 17. Not yet tested:

- **SAGE somatic on `24 * task.attempt` CPUs.** Every SAGE somatic
  task was cached in the resume. On 12 CPUs it used 97% of them and
  took up to 5.3 h. Compare its realtime per group with the 12-CPU
  runs; about 1.7 times faster is expected.
- **`executor.queueSize = 200`.** The resume submitted 19 jobs. Look
  for submission errors in `.nextflow.log` (SLURM allows a burst of
  30 commands per user, then about 2 per second), and for queue waits
  (trace `duration` minus `realtime`) during alignment.
- **Unused CPUs, if cmobic_cpu fills.** VIRUSBREAKEND uses about 3 of
  its 12 CPUs and REDUX about 5 of 12. This did not matter at
  `queueSize` 50 (no queue waits over 2 min). If tasks now wait for
  CPUs, consider lowering their `cpus`.
- **ESVEE on 216 GB for a group larger than seen so far.** The two
  largest groups (918 M and 935 M prep BAM records) passed. A larger
  one may be killed at the prep BAM sort (exit 137) and retried on
  432 GB, which costs about 4 h. Check for ESVEE retries.
- **PURPLE on tumor-only samples** peaked at 69 to 73 GB of 96 GB.
  Check for a PURPLE retry on a sample with a larger somatic VCF.

Optional change, not applied: LINX visualiser on 12 CPUs
(`withName: '.*:LINX_PLOTTING:LINX_VISUALISER' { cpus = 12 }`). It
keeps 0.93 of its 6 CPUs busy and ORANGE waits for it; one group took
58 min. Decide after seeing its times in the next run.

How to check: join `out/pipeline_info/execution_trace_*.txt` with
`sacct` on `native_id` (`ReqCPUS`, `ReqMem`, `Elapsed`, `TotalCPU`,
`State`, `QOS`). Use the trace's `peak_rss` for memory, not `sacct`
MaxRSS (section 17).
