# Paths and settings shared by the R scripts. All scripts are run from the repository root.

suppressPackageStartupMessages(library(nanoparquet))
source("code/lca_functions.R")

read_config <- function(file = "config.txt") {
  x <- trimws(readLines(file, warn = FALSE))
  x <- x[nzchar(x) & !startsWith(x, "#") & grepl("=", x, fixed = TRUE)]
  setNames(trimws(sub("^[^=]*=", "", x)), trimws(sub("=.*$", "", x)))
}
cfg <- read_config()
out_root <- if ("OUTPUT_DIR" %in% names(cfg)) cfg[["OUTPUT_DIR"]] else "output"
derived_dir <- file.path(out_root, "derived")
res_dir <- file.path(out_root, "results")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
res_path <- function(f) file.path(res_dir, f)

items <- c("ST438Q01DA", "ST438Q02DA", "ST438Q03DA", "ST438Q04DA")
K <- 6; R <- 5; M <- 10
class_names <- c("non-user", "occasional", "learning-focused",
                 "moderate all-purpose", "frequent all-purpose", "daily all-purpose")

# Students with at least one valid answer to the four AI-use items
load_sample <- function(response_style = FALSE) {
  d <- as.data.frame(read_parquet(file.path(derived_dir, "analysis_sample.parquet")))
  d <- d[d$lca_in == 1, ]
  if (response_style) {
    rs <- as.data.frame(read_parquet(file.path(derived_dir, "response_style.parquet")))
    d <- merge(d, rs, by = "CNTSTUID", all.x = TRUE, sort = FALSE)
  }
  d
}

# Retained six-class solution (weighted pseudo-ML), written by 03_enumerate_classes.R
load_reference <- function() readRDS(res_path("lca_fits.rds"))$fits[[K]]
