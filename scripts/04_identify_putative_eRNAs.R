# 04_identify_putative_eRNAs.R
# Identify putative eRNAs from groHMM output (MCF7 vs MCF10A)

library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# hg38 gene annotations
tx <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Input BED files
mcf7_grohmm_bed    <- "PATH_OR_URL_TO_MCF7_groHMM.bed"
mcf10a_grohmm_bed    <- "PATH_OR_URL_TO_MCF10A_groHMM.bed"
mcf7_h3k27ac_bed    <- "PATH_OR_URL_TO_MCF7_H3K27ac.bed"
mcf10a_h3k27ac_bed    <- "PATH_OR_URL_TO_MCF10A_H3K27ac.bed"

# Import data
mcf7_tx    <- import(mcf7_grohmm_bed)
mcf10a_tx    <- import(mcf10a_grohmm_bed)
mcf7_h3k27ac    <- import(mcf7_h3k27ac_bed)
mcf10a_h3k27ac    <- import(mc10a_h3k27ac_bed)

# Filter to intergenic transcripts
mcf7_intergenic    <- subsetByOverlaps(mcf7_tx, tx, invert = TRUE)
mcf10a_intergenic    <- subsetByOverlaps(mcf10a_tx, tx, invert = TRUE)

# Define putative eRNAs by H3K27ac overlap
mcf7_putative_eRNAs  <- subsetByOverlaps(mcf7_intergenic, mcf7_h3k27ac)
mcf10a_putative_eRNAs  <- subsetByOverlaps(mcf10a_intergenic, mcf10a_h3k27ac)

# Summary counts
cat("MCF7:", length(mcf7_tx),
    "-> intergenic:", length(mcf7_intergenic),
    "-> eRNAs:", length(mcf7_putative_eRNAs), "\n")

cat("MCF10A:", length(mcf10a_tx),
    "-> intergenic:", length(mcf10a_intergenic),
    "-> eRNAs:", length(mcf10a_putative_eRNAs), "\n")
