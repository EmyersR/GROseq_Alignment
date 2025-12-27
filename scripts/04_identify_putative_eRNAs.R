# 04_identify_putative_eRNAs.R
# Identify intergenic GROseq transcripts (putative eRNAs, step 1)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})

# Load hg38 gene annotations (to remove genic transcripts)
tx <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Input files (relative to repo root)
mcf7_bed <- "data/BC_eRNA_files/MCF7.final.transcripts.bed"
mcf10a_bed <- "data/BC_eRNA_files/MCF10A.final.transcripts.bed"

mcf7_h3k27ac_bed <- "data/BC_eRNA_files/ENCFF491LQY_MCF7_H3K27AC.bed"
mcf10a_h3k27ac_bed <- "data/BC_eRNA_files/ENCFF559BLN_MCF10A_H3K27AC.bed"

# Import GROseq transcripts
mcf7_tx   <- import(mcf7_bed)
mcf10a_tx <- import(mcf10a_bed)

# Remove transcripts overlapping known genes (keep intergenic only)
mcf7_intergenic   <- subsetByOverlaps(mcf7_tx, tx, invert = TRUE)
mcf10a_intergenic <- subsetByOverlaps(mcf10a_tx, tx, invert = TRUE)

# Create output directory
dir.create("data/processed", showWarnings = FALSE)

# Save intergenic transcripts
export(
  mcf7_intergenic,
  "data/processed/MCF7_intergenic_GROseq.bed"
)

export(
  mcf10a_intergenic,
  "data/processed/MCF10A_intergenic_GROseq.bed"
)

# Quick sanity check
cat("MCF7 intergenic transcripts:", length(mcf7_intergenic), "\n")
cat("MCF10A intergenic transcripts:", length(mcf10a_intergenic), "\n")