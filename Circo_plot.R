#!/usr/bin/env Rscript

# =============================
# HPV Integration Plot Generator
# =============================
# This script generates two types of circular plots:
# 1. Full circular plots with human + HPV genome integration
# 2. Simplified HPV-only circular plots
#
# Input files must be present in the working directory
# =============================

# Clean the environment
rm(list = ls())

# Load required libraries
suppressPackageStartupMessages({
  library(reshape)
  library(ggplot2)
  library(reshape2)
  library(grid)
  library(randomcoloR)
  library(magrittr)
  library(ggpubr)
  library(plotly)
  library(ggrepel)
  library(wesanderson)
  library(RColorBrewer)
  library(Formula)
  library(lattice)
  library(survival)
  library(Hmisc)
  library(ggsignif)
  library(dplyr)
  library(forcats)
  library(lsa)
  library(scales)
  library(Rcmdr)
  library(car)
  library(ggpmisc)
  library(tidyr)
  library(paletteer)
  library(circlize)
})

# Clean any existing graphics device
if (dev.cur() > 1) dev.off()
circos.clear()

# Load data
integration <- read.delim("./Inte.csv", header = TRUE, sep = "\t")
summary <- read.delim("./Resume_paris.csv", header = TRUE, sep = "\t")
clinical_data <- read.delim("./Datas_cliniques.txt", header = TRUE, sep = "\t")
sample_insertions <- read.delim("./Sample_insertions.csv", header = TRUE, sep = ",")
depth_samples <- na.omit(read.delim("./Full_hpv_coverage.txt", header = TRUE))
gene_annotations <- read.delim("./Full_infos.gb", header = FALSE)

# Merge integration and clinical data
combined <- merge(summary, integration, by.x = "best_capture", by.y = "sample")
combined <- merge(combined, clinical_data, by.x = "paire", by.y = "Paire")

# Extract unique pairs
unique_pairs <- unique(depth_samples$New)

# Define color palette
palette_colors <- paletteer_d("MoMAColors::Klein")


# =============================
# SECTION 1: Full Circular Plots (Human + HPV)
# =============================

for (pair_id in unique_pairs) {
  
  cat("Processing full circular plot for pair:", pair_id, "\n")
  
  # Get HPV genome size for this pair
  depth_adc <- depth_samples %>%
    filter(New == pair_id, Cancer == "ADC")
  
  if (nrow(depth_adc) == 0) next
  hpv_size <- max(depth_adc$Position)
  
  # Build BED-like ADC and SCC dataframes
  bed_adc <- data.frame(
    chr = "chrHPV",
    start = depth_adc$Position * (3.4e9 / hpv_size),
    stop = depth_adc$Position * (3.4e9 / hpv_size) + (3.4e9 / hpv_size),
    value1 = depth_adc$Depht
  )
  
  bed_scc <- depth_samples %>%
    filter(New == pair_id, Cancer == "SCC") %>%
    transmute(
      chr = "chrHPV",
      start = Position * (3.4e9 / hpv_size),
      stop = Position * (3.4e9 / hpv_size) + (3.4e9 / hpv_size),
      value1 = Depht
    )
  
  # Prepare cytoband data
  human_cytoband <- read.cytoband(species = "hg38")$df
  hpv_cytoband <- data.frame("chrHPV", 0, 3.4e9, "q12", "Blue")
  colnames(hpv_cytoband) <- colnames(human_cytoband)
  cytoband <- rbind(human_cytoband, hpv_cytoband) %>%
    filter(V1 != "chrY")
  
  chromosome.index <- c(paste0("chr", c(1:22, "X")), "chrHPV")
  
  # Sample insertions for this pair
  sample_adc <- sample_insertions %>%
    filter(sample == unique(depth_adc$SeqId))
  
  sample_scc <- sample_insertions %>%
    filter(sample == unique(depth_samples$SeqId[depth_samples$New == pair_id & depth_samples$Cancer == "SCC"]))
  
  # Prepare BEDs for links
  bed1 <- sample_adc[, c("chr", "chr_position")] %>%
    rename(start = chr_position) %>%
    mutate(end = start + 2e6)
  
  bed2 <- sample_adc[, c("chr", "position")] %>%
    rename(start = position) %>%
    mutate(
      chr = "chrHPV",
      start = start * (3.4e9 / hpv_size),
      end = (start + 50) * (3.4e9 / hpv_size)
    ) %>%
    select(chr, start, end)
  
  bed5 <- data.frame(
    chr = rep("chrHPV", nrow(bed2)),
    start = rep(1, nrow(bed2)),
    end = rep(3.4e9, nrow(bed2))
  )
  
  # HPV gene annotations
  strain <- unique(sample_adc$Strain)
  gene_data <- gene_annotations %>%
    filter(V1 == strain)
  
  hpv_cytoband2 <- data.frame(
    V1 = gene_data$V4,
    V2 = gene_data$V2 * (3.4e9 / hpv_size),
    V3 = gene_data$V3 * (3.4e9 / hpv_size),
    V4 = "q12",
    V5 = palette_colors[1:nrow(gene_data)]
  )
  colnames(hpv_cytoband2) <- colnames(human_cytoband)
  
  # Add y-coordinates for gene boxes
  hpv_cytoband2 <- hpv_cytoband2 %>%
    arrange(V2) %>%
    mutate(
      row_id = row_number(),
      y1 = if_else(row_id %% 2 == 0, 0, 0.5),
      y2 = if_else(row_id %% 2 == 0, 0.5, 1)
    ) %>%
    select(-row_id)
  
  # Start PDF output
  pdf(paste0("full_circular_pair_", pair_id, "_ADC.pdf"))
  
  # Setup circular plot
  circos.clear()
  circos.par(gap.after = c(rep(1, 22), 4, 4), start.degree = 90)
  circos.initializeWithIdeogram(cytoband, plotType = NULL, chromosome.index = chromosome.index)
  
  # ADC track
  circos.genomicTrackPlotRegion(bed_adc, ylim = c(0, max(bed_adc$value1, na.rm = TRUE)),
    panel.fun = function(region, value, ...) {
      circos.genomicLines(region, value, area = TRUE,
                          col = "#FFCCCC", border = "#FF6666", lwd = 0.5, ...)
    }, bg.border = NA)
  
  # SCC track
  if (nrow(bed_scc) > 0) {
    circos.genomicTrackPlotRegion(bed_scc, ylim = c(0, max(bed_scc$value1, na.rm = TRUE)),
      panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value, area = TRUE,
                            col = "#CCE5FF", border = "#3399ff", lwd = 0.5, ...)
      }, bg.border = NA)
  }
  
  # Add cytoband labels
  circos.track(ylim = c(0, 0.001), panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2),
                gsub(".*chr", "", CELL_META$sector.index), cex = 0.6,
                niceFacing = TRUE, facing = "outside")
  }, track.height = mm_h(1), cell.padding = c(0, 0, 0, 0), bg.border = NA)
  
  # Plot gene boxes
  for (y in seq_len(nrow(hpv_cytoband2))) {
    circos.rect(sector.index = "chrHPV",
                xleft = hpv_cytoband2$V2[y], xright = hpv_cytoband2$V3[y],
                ybottom = hpv_cytoband2$y1[y], ytop = hpv_cytoband2$y2[y],
                col = palette_colors[y])
  }
  
  # Gene labels
  for (z in seq_len(nrow(hpv_cytoband2))) {
    circos.text((hpv_cytoband2$V2[z] + hpv_cytoband2$V3[z]) / 2,
                hpv_cytoband2$y1[z] + 0.25,
                hpv_cytoband2$V1[z], cex = 0.7, col = "white",
                facing = "inside", niceFacing = TRUE)
  }

  # Annotations
  text(0, 1.05, paste("Pair", pair_id), cex = 1.5)
  text(0.05, 0.82, "ADC", cex = 1, col = "#FF6666")
  text(0.05, 0.62, "SCC", cex = 1, col = "#3399ff")
  text(0, 0, gsub("_", "-", strain), cex = 1.25, col = "gray")

  # Finalize
  circos.clear()
  dev.off()
}


# =============================
# SECTION 2: HPV-Only Circular Plots
# =============================

for (pair_id in unique_pairs) {
  
  cat("Processing HPV-only circular plot for pair:", pair_id, "\n")
  
  # Get HPV size for this pair (ADC only)
  depth_adc <- depth_samples %>%
    filter(New == pair_id, Cancer == "ADC")
  if (nrow(depth_adc) == 0) next
  
  hpv_size <- max(depth_adc$Position)
  
  # Build ADC and SCC tracks without rescaling to human genome
  bed_adc <- depth_adc %>%
    transmute(chr = "chrHPV",
              start = Position,
              stop = Position + 1,
              value1 = Depht)
  
  bed_scc <- depth_samples %>%
    filter(New == pair_id, Cancer == "SCC") %>%
    transmute(chr = "chrHPV",
              start = Position,
              stop = Position + 1,
              value1 = Depht)
  
  # Cytoband: HPV genome only
  human_cytoband <- read.cytoband(species = "hg38")$df
  hpv_cytoband <- data.frame("chrHPV", 0, hpv_size, "q12", "gpos1")
  colnames(hpv_cytoband) <- colnames(human_cytoband)
  
  # HPV subtype
  subtype <- unique(depth_adc$Subtype)
  gene_data <- gene_annotations %>%
    filter(V1 == subtype)
  
  hpv_cytoband2 <- data.frame(
    V1 = gene_data$V4,
    V2 = gene_data$V2,
    V3 = gene_data$V3,
    V4 = "q12",
    V5 = palette_colors[1:nrow(gene_data)]
  )
  colnames(hpv_cytoband2) <- colnames(human_cytoband)
  
  # Alternate y positions for genes
  hpv_cytoband2 <- hpv_cytoband2 %>%
    arrange(V2) %>%
    mutate(
      row_id = row_number(),
      y1 = if_else(row_id %% 2 == 0, 0, 0.5),
      y2 = if_else(row_id %% 2 == 0, 0.5, 1)
    ) %>%
    select(-row_id)
  
  # Sample insertions
  sample_adc <- sample_insertions %>%
    filter(sample == unique(depth_adc$SeqId))
  sample_scc <- sample_insertions %>%
    filter(sample == unique(depth_samples$SeqId[depth_samples$New == pair_id & depth_samples$Cancer == "SCC"]))
  
  # Build insertion arrows (ADC)
  if (nrow(sample_adc) > 0) {
    bed_adc_ins <- data.frame(chr = "chrHPV", start = sample_adc$position)
    bed_adc_ins$end <- bed_adc_ins$start + 50
    bed_adc_ins$feature <- sample_adc$feature
  }
  
  # Build insertion arrows (SCC)
  if (nrow(sample_scc) > 0) {
    bed_scc_ins <- data.frame(chr = "chrHPV", start = sample_scc$position)
    bed_scc_ins$end <- bed_scc_ins$start + 50
    bed_scc_ins$feature <- sample_scc$feature
  }
  
  # Output PDF
  pdf(paste0("hpv_only_circular_pair_", pair_id, ".pdf"))
  
  # Setup circular genome
  circos.clear()
  circos.par(gap.after = 5, start.degree = 90)
  circos.initializeCircularGenome(
    name = hpv_cytoband$V1,
    genome_size = hpv_cytoband$V3,
    plotType = c("axis", "labels")
  )
  
  # ADC track
  circos.genomicTrackPlotRegion(bed_adc, ylim = c(0, max(bed_adc$value1, na.rm = TRUE)),
    panel.fun = function(region, value, ...) {
      circos.genomicLines(region, value, area = TRUE,
                          col = "#FFCCCC", border = "#FF6666", lwd = 0.5, ...)
    }, bg.border = NA)
  
  # ADC arrows
  if (exists("bed_adc_ins")) {
    for (h in 1:nrow(bed_adc_ins)) {
      ins <- bed_adc_ins[h, ]
      circos.segments(ins$start, 0, ins$start, max(bed_adc$value1, na.rm = TRUE))
      
      if (ins$feature == "left") {
        circos.arrow(ins$start, ins$start + 100,
                     arrow.head.width = 2 * (max(bed_adc$value1, na.rm = TRUE) / 3),
                     max(bed_adc$value1, na.rm = TRUE),
                     col = "#FF6666")
      }
      
      if (ins$feature == "right") {
        circos.arrow(ins$start - 100, ins$start,
                     arrow.head.width = 2 * (max(bed_adc$value1, na.rm = TRUE) / 3),
                     max(bed_adc$value1, na.rm = TRUE),
                     arrow.position = "start",
                     col = "#FF6666")
      }
    }
  }
  
  # SCC track
  if (nrow(bed_scc) > 0) {
    circos.genomicTrackPlotRegion(bed_scc, ylim = c(0, max(bed_scc$value1, na.rm = TRUE)),
      panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value, area = TRUE,
                            col = "#CCE5FF", border = "#3399ff", lwd = 0.5, ...)
      }, bg.border = NA)
  }
  
  # SCC arrows
  if (exists("bed_scc_ins")) {
    for (h in 1:nrow(bed_scc_ins)) {
      ins <- bed_scc_ins[h, ]
      circos.segments(ins$start, 0, ins$start, max(bed_scc$value1, na.rm = TRUE))
      
      if (ins$feature == "left") {
        circos.arrow(ins$start, ins$start + 130,
                     arrow.head.width = 2 * (max(bed_scc$value1, na.rm = TRUE) / 3),
                     max(bed_scc$value1, na.rm = TRUE),
                     col = "#3399ff")
      }
      
      if (ins$feature == "right") {
        circos.arrow(ins$start - 130, ins$start,
                     arrow.head.width = 2 * (max(bed_scc$value1, na.rm = TRUE) / 3),
                     max(bed_scc$value1, na.rm = TRUE),
                     arrow.position = "start",
                     col = "#3399ff")
      }
    }
  }
  
  # Add gene rectangles
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    xlim <- CELL_META$xlim
    circos.rect(xlim[1], 0, xlim[2], 1, col = "white")
  }, track.height = 0.15, bg.border = NA)
  
  for (y in 1:nrow(hpv_cytoband2)) {
    circos.rect(sector.index = "chrHPV",
                xleft = hpv_cytoband2$V2[y], xright = hpv_cytoband2$V3[y],
                ybottom = hpv_cytoband2$y1[y], ytop = hpv_cytoband2$y2[y],
                col = palette_colors[y])
  }
  
  # Gene labels
  for (z in 1:nrow(hpv_cytoband2)) {
    circos.text((hpv_cytoband2$V2[z] + hpv_cytoband2$V3[z]) / 2,
                hpv_cytoband2$y1[z] + 0.25,
                hpv_cytoband2$V1[z], cex = 0.7, col = "white",
                facing = "inside", niceFacing = TRUE)
  }
  
  # Plot labels
  text(0, 1.05, paste("Pair", pair_id), cex = 1.5)
  text(0.05, 0.82, "ADC", cex = 1, col = "#FF6666")
  text(0.05, 0.62, "SCC", cex = 1, col = "#3399ff")
  text(0, 0, gsub("_", "-", subtype), cex = 1.25, col = "gray")
  
  # Finalize
  circos.clear()
  dev.off()
}
