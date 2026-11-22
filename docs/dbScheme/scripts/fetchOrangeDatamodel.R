#!/usr/bin/env Rscript
# Prism dbScheme: download the ORANGE data model (the Java classes that
# define the ORANGE JSON) for one hmftools release tag.
#
#   Rscript fetchOrangeDatamodel.R TAG DEST_DIR
#
# TAG is an hmftools release tag, e.g. orange-v5.0.1 (oncoanalyser 3.0.0
# runs ORANGE 5.0.1; the version is in the pipeline's
# modules/local/orange/main.nf). Writes one .java file per class into
# DEST_DIR plus SOURCE.txt with the tag and the download date. One GitHub
# API call (the file list), then plain downloads.

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
  library(jsonlite)
})

REPO <- "hartwigmedical/hmftools"
MODEL_DIR <- "orange-datamodel/src/main/java/"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: fetchOrangeDatamodel.R TAG DEST_DIR", call. = FALSE)
}
tag <- args[1]
dest <- path_abs(args[2])
dir_create(dest)

tree <- fromJSON(str_glue(
  "https://api.github.com/repos/{REPO}/git/trees/{tag}?recursive=1"
))
files <- tree$tree |>
  as_tibble() |>
  filter(type == "blob", str_starts(path, MODEL_DIR), str_ends(path, "\\.java")) |>
  mutate(name = path_file(path))
if (nrow(files) == 0) {
  stop(str_glue("No data model files at tag {tag}"), call. = FALSE)
}
if (anyDuplicated(files$name) > 0) {
  stop("Two data model classes share a file name", call. = FALSE)
}

walk2(files$path, files$name, \(p, name) {
  download.file(
    str_glue("https://raw.githubusercontent.com/{REPO}/{tag}/{p}"),
    path(dest, name), quiet = TRUE
  )
})
write_lines(
  c(str_glue("repo: {REPO}"), str_glue("tag: {tag}"),
    str_glue("downloaded: {Sys.Date()}"), str_glue("files: {nrow(files)}")),
  path(dest, "SOURCE.txt")
)
message(str_glue("{nrow(files)} files from {REPO} at {tag} in {dest}"))
