# 03_grohmm_calls.R
# Identify nascent transcripts from GRO-Seq data using GroHMM

library(groHMM)
library(GenomicRanges)
library(rtracklayer)

# Load GRO-Seq GRanges from earlier processed step or BEDs directly
mcf7_path   <- "data/processed/mcf7_granges.rds"
mcf10a_path <- "data/processed/mcf10a_granges.rds"

if (!file.exists(mcf7_path) || !file.exists(mcf10a_path)) {
  message("Missing GRO-Seq GRanges files. Run your import script first.")
  quit(save = "no")
}

mcf7   <- readRDS(mcf7_path)
mcf10a <- readRDS(mcf10a_path)

# Convert to coverage objects
cov_mcf7   <- coverage(mcf7)
cov_mcf10a <- coverage(mcf10a)

# Run GroHMM to detect nascent transcripts
calls_mcf7   <- detectTranscripts(cov_mcf7, threshold = 2, maxLen = 50000L)
calls_mcf10a <- detectTranscripts(cov_mcf10a, threshold = 2, maxLen = 50000L)

# Save outputs
dir.create("data/processed", showWarnings = FALSE)
saveRDS(calls_mcf7,   "data/processed/grohmm_calls_mcf7.rds")
saveRDS(calls_mcf10a, "data/processed/grohmm_calls_mcf10a.rds")

rtracklayer::export(as(calls_mcf7, "GRanges"),   "data/processed/grohmm_mcf7.bed")
rtracklayer::export(as(calls_mcf10a, "GRanges"), "data/processed/grohmm_mcf10a.bed")

message("GroHMM transcript calling complete.")
