# 02_qc_plots.R — GRO-seq QC Plots

# Load required packages
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)

# Load processed GRO-seq GRanges objects
mcf7 <- readRDS("data/processed/mcf7_transcripts.rds")
mcf10a <- readRDS("data/processed/mcf10a_transcripts.rds")

# Summarize read distributions
df_summary <- data.frame(
  Sample = c("MCF7", "MCF10A"),
  Features = c(length(mcf7), length(mcf10a)),
  MedianLength = c(median(width(mcf7)), median(width(mcf10a))),
  MaxLength = c(max(width(mcf7)), max(width(mcf10a)))
)

print(df_summary)

# Plot read length distributions
ggplot() +
  geom_density(aes(x = width(mcf7)), color = "red", size = 1) +
  geom_density(aes(x = width(mcf10a)), color = "blue", size = 1) +
  theme_minimal() +
  labs(title = "GRO-seq Read Length Distribution",
       x = "Read Length (bp)", y = "Density",
       caption = "MCF7 (red) vs MCF10A (blue)") +
  theme(plot.title = element_text(hjust = 0.5))
# Plot read length distributions
ggplot() +
  geom_density(aes(x = width(mcf7)), color = "red", size = 1) +
  geom_density(aes(x = width(mcf10a)), color = "blue", size = 1) +
  theme_minimal() +
  labs(
    title = "GRO-seq Read Length Distribution",
    x = "Read Length (bp)",
    y = "Density",
    caption = "MCF7 (red) vs MCF10A (blue)"
  ) +
  theme(plot.title = element_text(hjust = 0.5))
