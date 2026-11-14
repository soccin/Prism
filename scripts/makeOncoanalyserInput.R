#!/usr/bin/env Rscript
#
# Build the nf-core/oncoanalyser samplesheet from mapping.csv and
# manifest.csv.
#
# Usage:
#   Rscript Prism/scripts/makeOncoanalyserInput.R [mapping.csv] [manifest.csv] [output.csv]
#
# Defaults: mapping.csv, manifest.csv and oncoanalyser_input.csv in the
# current directory. Templates are in Prism/examples.
#
# mapping.csv  : sample_id,fastq_dir,filetype,sequence_type
#                one row per sample per sequencing run; fastq_dir is an
#                IGO Sample_* folder holding *_L00N_R1_001.fastq.gz / R2 pairs
# manifest.csv : sample_id,sample_type,subject_id,group_id
#
# Output: one row per FASTQ pair (sample x flowcell x lane) with
#   info = library_id:<IGO id>;lane:<NNN>;flowcell:<flowcell id>
# The flowcell is needed because lane numbers repeat across flowcells and
# oncoanalyser builds read group IDs from sample, library, lane, flowcell.

suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(readr)
    library(stringr)
    library(tidyr)
    library(fs)
    library(glue)
})

args <- commandArgs(trailingOnly = TRUE)
mapping_file <- if (length(args) >= 1) args[1] else "mapping.csv"
manifest_file <- if (length(args) >= 2) args[2] else "manifest.csv"
output_file <- if (length(args) >= 3) args[3] else "oncoanalyser_input.csv"

stop_if <- function(cond, msg) {
    if (cond) stop(msg, call. = FALSE)
}

mapping <- read_csv(mapping_file, show_col_types = FALSE)
manifest <- read_csv(manifest_file, show_col_types = FALSE)

stop_if(
    !all(c("sample_id", "fastq_dir", "filetype", "sequence_type") %in% names(mapping)),
    glue("{mapping_file} must have columns sample_id,fastq_dir,filetype,sequence_type")
)
stop_if(
    !all(c("sample_id", "sample_type", "subject_id", "group_id") %in% names(manifest)),
    glue("{manifest_file} must have columns sample_id,sample_type,subject_id,group_id")
)
stop_if(
    any(mapping$filetype != "fastq"),
    "This script only handles filetype 'fastq'"
)

missing_in_manifest <- setdiff(mapping$sample_id, manifest$sample_id)
stop_if(
    length(missing_in_manifest) > 0,
    glue("Samples in mapping but not manifest: {str_c(missing_in_manifest, collapse = ', ')}")
)
missing_in_mapping <- setdiff(manifest$sample_id, mapping$sample_id)
stop_if(
    length(missing_in_mapping) > 0,
    glue("Samples in manifest but not mapping: {str_c(missing_in_mapping, collapse = ', ')}")
)

bad_groups <- manifest |>
    distinct(group_id, subject_id) |>
    count(group_id) |>
    filter(n > 1)
stop_if(
    nrow(bad_groups) > 0,
    glue("group_id with more than one subject_id: {str_c(bad_groups$group_id, collapse = ', ')}")
)

missing_dirs <- mapping$fastq_dir[!dir_exists(mapping$fastq_dir)]
stop_if(
    length(missing_dirs) > 0,
    glue("Missing fastq_dir: {str_c(missing_dirs, collapse = '\n')}")
)

# One row per R1 file found in each fastq_dir
fastqs <- mapping |>
    mutate(r1 = map(fastq_dir, \(d) as.character(dir_ls(d, regexp = "_R1_001\\.fastq\\.gz$")))) |>
    unnest(r1) |>
    mutate(
        r2 = str_replace(r1, "_R1_001\\.fastq\\.gz$", "_R2_001.fastq.gz"),
        lane = str_extract(path_file(r1), "(?<=_L)\\d{3}(?=_R1_001)"),
        library_id = str_extract(path_file(fastq_dir), "IGO_.*$"),
        flowcell = path_dir(path_dir(fastq_dir)) |>
            path_file() |>
            str_extract("[^_]+$")
    )

empty_dirs <- setdiff(mapping$fastq_dir, fastqs$fastq_dir)
stop_if(
    length(empty_dirs) > 0,
    glue("No *_R1_001.fastq.gz in: {str_c(empty_dirs, collapse = '\n')}")
)
stop_if(
    any(!file_exists(fastqs$r2)),
    glue("Missing R2 for: {str_c(fastqs$r1[!file_exists(fastqs$r2)], collapse = '\n')}")
)
stop_if(
    any(is.na(fastqs$lane) | is.na(fastqs$library_id) | is.na(fastqs$flowcell)),
    "Could not parse lane, library_id or flowcell for some FASTQ files"
)

dups <- fastqs |>
    count(sample_id, library_id, flowcell, lane) |>
    filter(n > 1)
stop_if(
    nrow(dups) > 0,
    "Duplicate sample/library/flowcell/lane combinations; read group IDs would clash"
)

samplesheet <- fastqs |>
    left_join(manifest, by = "sample_id") |>
    mutate(
        info = glue("library_id:{library_id};lane:{lane};flowcell:{flowcell}"),
        filepath = str_c(r1, ";", r2),
        type_order = match(sample_type, c("normal", "donor", "tumor"))
    ) |>
    arrange(group_id, type_order, sample_id, flowcell, lane) |>
    select(
        group_id, subject_id, sample_id, sample_type, sequence_type,
        filetype, info, filepath
    )

write_csv(samplesheet, output_file, na = "")

cat(glue("Wrote {nrow(samplesheet)} rows to {output_file}"), "\n")
print(count(samplesheet, group_id, sample_id, sample_type, name = "fastq_pairs"))
