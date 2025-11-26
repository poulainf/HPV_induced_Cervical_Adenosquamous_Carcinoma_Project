#!/usr/bin/Rscript

################################################################################
# Title      : Mutect2 Variant Filtering and Annotation Script
# Description: This script processes somatic variants from Mutect2 VCF files,
#              applies quality filters, annotates with CADD scores, merges
#              clinical metadata, and recovers cross-sample variants.
# Author     : Dr. Florian Poulain
# Date       : 2025-11-26
################################################################################

#---------------------------#
# Environment Setup
#---------------------------#

rm(list = ls())  # Clear R environment

# Set working directory
setwd("/media/florian/T7Shield/Projet_Liege/")
# setwd("~/Projet_Liege/")  # Alternate for local use

#---------------------------#
# Load Required Libraries
#---------------------------#

# Author: Dr. Florian Poulain | Date: 2025-11-26
library(data.table)   # Fast file reading and interval joins
library(dplyr)        # Data manipulation (filter, mutate, join, etc.)
library(stringr)      # Regex-based string processing
library(tidyr)        # Data cleaning, handling missing values

#---------------------------#
# Load and Preprocess Raw Mutect2 VCF Data
#---------------------------#

raw_vcf <- fread("/media/florian2/T7/To_send/VCF_Muect2/Raw_VCF.vcf", sep = "\t", header = FALSE, nThread = 4)
raw_vcf$Paire <- as.numeric(gsub("Paire_(\\d+)_.*", "\\1", raw_vcf$V1))

#---------------------------#
# Load MAF and Clean Sample Names
#---------------------------#

maf <- read.delim("/media/florian2/T7/To_send/VCF_Muect2/test_SNPs_MAF.txt", header = TRUE, sep = "\t")

maf$Tumor_Sample_Barcode <- gsub("CIN3|SCC2|SCC3", "SCC", maf$Tumor_Sample_Barcode)
maf$Tumor_Sample_Barcode <- gsub("ADC2", "ADC", maf$Tumor_Sample_Barcode)

maf$Paire <- gsub("Paire_(\\d+)_.*", "\\1", maf$Tumor_Sample_Barcode)
maf$Cancer <- gsub("Paire_\\d+_(.+)_GATK_somatic_filtered", "\\1", maf$Tumor_Sample_Barcode)
maf$Tumor_Sample_Barcode <- gsub("_GATK_somatic_filtered", "", maf$Tumor_Sample_Barcode)

maf <- maf[!maf$Paire %in% c(5, 7, 9, 10), ]

#---------------------------#
# Annotate Intronic SNPs as Non-Exonic
#---------------------------#

intronic_snps <- maf[maf$Variant_Classification == "" | maf$Variant_Classification == "Unknown", ]

exons_bed <- read.delim("Twist_Comprehensive_Exome_Covered_Targets_hg38.bed", header = FALSE, sep = "\t")
exons_bed$up <- exons_bed$V2 - 120
exons_bed$down <- exons_bed$V2 + 120
exons_bed$chr <- exons_bed$V1

exons <- as.data.table(exons_bed[, c("chr", "up", "down")])
setnames(exons, c("chr", "up", "down"), c("Chromosome", "start", "end"))

snps <- as.data.table(intronic_snps)
snps[, `:=`(start = Start_Position, end = Start_Position)]

setkey(exons, Chromosome, start, end)
setkey(snps, Chromosome, start, end)

non_exonic_snps <- data.frame(foverlaps(snps, exons, nomatch = 0))
non_exonic_snps <- non_exonic_snps[, intersect(colnames(maf), colnames(non_exonic_snps))]
non_exonic_snps$Variant_Classification <- "Non exonic"

maf_cleaned <- maf[!maf$Variant_Classification %in% c("", "Unknown"), ]
maf <- rbind(maf_cleaned, non_exonic_snps)

#---------------------------#
# Apply Quality Filtering on Variants
#---------------------------#

maf$Depht_tumor <- as.numeric(gsub(".*:.*:.*:([^:,]+):.*", "\\1", maf$Otherinfo13))
maf$Depht_normal <- as.numeric(gsub(".*:.*:.*:([^:,]+):.*", "\\1", maf$Otherinfo14))
maf$TLOD <- as.numeric(gsub(".*TLOD=([0-9.]+).*", "\\1", maf$Otherinfo11))

maf <- maf[maf$Otherinfo10 %in% c("PASS", "clustered_events", "haplotype"), ]
maf <- maf[maf$Depht_tumor >= 30 & maf$Depht_normal >= 20 & maf$TLOD >= 10, ]

maf$VAF_T <- as.numeric(gsub(".*:.*:(.*?):.*", "\\1", maf$Otherinfo13))
maf$VAF_N <- as.numeric(gsub(".*:.*:(.*?):.*", "\\1", maf$Otherinfo14))
maf <- maf[maf$VAF_T >= 0.08 & maf$VAF_N < 0.04, ]

maf$AD_T_ALT <- as.numeric(gsub("^[^:]+:[^,]+,([^:]+):.*", "\\1", maf$Otherinfo13))
maf <- maf[maf$AD_T_ALT >= 4, ]

#---------------------------#
# Recover Cross-Cancer Variants in Matched Pairs
#---------------------------#

raw_vcf$V1 <- gsub("CIN3|SCC2", "SCC", raw_vcf$V1)
raw_vcf$V1 <- gsub("ADC2", "ADC", raw_vcf$V1)
raw_vcf$V1[raw_vcf$V1 == "Paire_2_ADC"] <- "Paire_18_ADC"
raw_vcf$Cancer <- gsub("Paire_\\d+_(.+)", "\\1", raw_vcf$V1)

raw_vcf$combi <- paste(raw_vcf$V1, raw_vcf$V2, raw_vcf$V3, raw_vcf$V5, raw_vcf$V6)
maf$combi <- paste(paste0("Paire_", maf$Paire, "_", maf$Cancer), maf$Chromosome, maf$Start_Position, maf$Reference_Allele, maf$Tumor_Seq_Allele2)

maf <- maf[maf$Paire %in% c("11", "12", "13", "14", "16", "15", "18", "1", "3", "4", "6", "8"), ]

recovered_variants <- data.frame()
for (p in unique(maf$Paire)) {
  adc <- maf[maf$Paire == p & maf$Cancer == "ADC", ]
  scc <- maf[maf$Paire == p & maf$Cancer == "SCC", ]
  
  adc$research <- gsub("ADC", "SCC", adc$combi)
  scc$research <- gsub("SCC", "ADC", scc$combi)
  
  rec_adc <- adc[adc$research %in% raw_vcf$combi, ]
  rec_scc <- scc[scc$research %in% raw_vcf$combi, ]
  
  if (nrow(rec_adc) > 0) {
    idx <- match(rec_adc$research, raw_vcf$combi)
    rec_adc$Cancer <- "SCC"
    rec_adc$Tumor_Sample_Barcode <- raw_vcf$V1[idx]
    rec_adc$Otherinfo11 <- raw_vcf[[9]][idx]
    rec_adc$Otherinfo13 <- raw_vcf[[11]][idx]
    rec_adc$Otherinfo14 <- raw_vcf[[12]][idx]
    rec_adc$Paire <- raw_vcf[[13]][idx]
    rec_adc$research <- NULL
  }
  
  if (nrow(rec_scc) > 0) {
    idx <- match(rec_scc$research, raw_vcf$combi)
    rec_scc$Cancer <- "ADC"
    rec_scc$Tumor_Sample_Barcode <- raw_vcf$V1[idx]
    rec_scc$Otherinfo11 <- raw_vcf[[9]][idx]
    rec_scc$Otherinfo13 <- raw_vcf[[11]][idx]
    rec_scc$Otherinfo14 <- raw_vcf[[12]][idx]
    rec_scc$Paire <- raw_vcf[[13]][idx]
    rec_scc$research <- NULL
  }
  
  recovered_variants <- rbind(recovered_variants, rec_adc, rec_scc)
}

recovered_variants$Fishing <- "Yes"
maf$Fishing <- "No"
maf <- unique(rbind(maf, recovered_variants))

#---------------------------#
# Final Cleanup
#---------------------------#

# Remove known artifact
maf <- maf[maf$Reference_Allele != "GCGGCCGCCGCCGCCGCCGCTGCGGGCGGCGCGCACCAGAACTCGGCCGTGGCGGCGGCGGCGGCGGCG", ]

# Merge clinical metadata
clinical <- read.delim("Datas_cliniques.txt", header = TRUE, sep = "\t")
maf <- merge(clinical, maf, by = "Paire")
maf <- unique(maf)

#---------------------------#
# Export Final Annotated MAF
#---------------------------#

write.table(maf, file = "Mutect2_VCF0.8.maf", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
