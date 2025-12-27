# 04_identify_putative_eRNAs.R
# Goal: identify intergenic GROseq transcripts that overlap H3K27ac (putative eRNAs)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})

# Gene annotations (hg38) for removing genic transcripts
tx <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Input files (relative to repo root)
mcf7_bed        <- "data/BC_eRNA_files/MCF7.final.transcripts.bed"
mcf10a_bed      <- "data/BC_eRNA_files/MCF10A.final.transcripts.bed"
mcf7_h3k27ac_bed  <- "data/BC_eRNA_files/ENCFF491LQY_MCF7_H3K27AC.bed"
mcf10a_h3k27ac_bed <- "data/BC_eRNA_files/ENCFF559BLN_MCF10A_H3K27AC.bed"

# Import GROseq transcripts (BED -> GRanges)
mcf7_tx   <- import(mcf7_bed)
mcf10a_tx <- import(mcf10a_bed)

# Import H3K27ac peaks (ENCODE BEDs can have decimal "score" -> read as table, use first 3 cols)
mcf7_h3k27ac_df <- read.table(
  mcf7_h3k27ac_bed, sep = "\t", header = FALSE,
  comment.char = "", quote = "", stringsAsFactors = FALSE
)
mcf10a_h3k27ac_df <- read.table(
  mcf10a_h3k27ac_bed, sep = "\t", header = FALSE,
  comment.char = "", quote = "", stringsAsFactors = FALSE
)

mcf7_h3k27ac <- GRanges(
  seqnames = mcf7_h3k27ac_df[[1]],
  ranges   = IRanges(start = mcf7_h3k27ac_df[[2]] + 1, end = mcf7_h3k27ac_df[[3]])
)

mcf10a_h3k27ac <- GRanges(
  seqnames = mcf10a_h3k27ac_df[[1]],
  ranges   = IRanges(start = mcf10a_h3k27ac_df[[2]] + 1, end = mcf10a_h3k27ac_df[[3]])
)

# Keep intergenic only (remove anything overlapping known genes)
mcf7_intergenic   <- subsetByOverlaps(mcf7_tx, tx, invert = TRUE)
mcf10a_intergenic <- subsetByOverlaps(mcf10a_tx, tx, invert = TRUE)

# Putative eRNAs = intergenic transcripts overlapping H3K27ac
mcf7_putative_eRNAs   <- subsetByOverlaps(mcf7_intergenic, mcf7_h3k27ac)
mcf10a_putative_eRNAs <- subsetByOverlaps(mcf10a_intergenic, mcf10a_h3k27ac)

# Save outputs
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

export(mcf7_intergenic,        "data/processed/MCF7_intergenic_GROseq.bed")
export(mcf10a_intergenic,      "data/processed/MCF10A_intergenic_GROseq.bed")
export(mcf7_putative_eRNAs,    "data/processed/MCF7_putative_eRNAs_H3K27ac.bed")
export(mcf10a_putative_eRNAs,  "data/processed/MCF10A_putative_eRNAs_H3K27ac.bed")

# Progress check
cat("MCF7 intergenic transcripts:", length(mcf7_intergenic), "\n")
cat("MCF10A intergenic transcripts:", length(mcf10a_intergenic), "\n")
cat("MCF7 putative eRNAs (H3K27ac overlap):", length(mcf7_putative_eRNAs), "\n")
cat("MCF10A putative eRNAs (H3K27ac overlap):", length(mcf10a_putative_eRNAs), "\n")

