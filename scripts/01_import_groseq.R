# 01_import_groseq.R — Import GRO-Seq BigWig files and save GRanges

library(rtracklayer)
library(GenomicRanges)
library(tidyverse)

# Define paths (BigWig)
mcf7_bw   <- list.files("data/GROseq/MCF7",   pattern="\\.bw$", full.names=TRUE, ignore.case=TRUE)
mcf10a_bw <- list.files("data/GROseq/MCF10A", pattern="\\.bw$", full.names=TRUE, ignore.case=TRUE)

if (length(mcf7_bw) == 0 || length(mcf10a_bw) == 0) {
  stop("Missing BigWig files in data/GROseq/MCF7 or data/GROseq/MCF10A")
}

# Import + combine strands
import_bw_pair <- function(bw_files) {
  gr_list <- lapply(bw_files, rtracklayer::import, format = "BigWig")
  combined <- reduce(do.call(c, gr_list))
  combined
}

mcf7   <- import_bw_pair(mcf7_bw)
mcf10a <- import_bw_pair(mcf10a_bw)

# Quick QC
qc <- tibble(
  sample       = c("MCF7","MCF10A"),
  n_regions    = c(length(mcf7), length(mcf10a)),
  median_width = c(median(width(mcf7)), median(width(mcf10a)))
)

# Save outputs
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(mcf7,   "data/processed/mcf7_granges.rds")
saveRDS(mcf10a, "data/processed/mcf10a_granges.rds")
readr::write_csv(qc, "data/processed/groseq_qc_summary.csv")

message("GRO-Seq BigWig import complete — GRanges saved to data/processed/")
