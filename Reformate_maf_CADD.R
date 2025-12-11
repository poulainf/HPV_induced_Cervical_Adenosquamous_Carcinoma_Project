#!/usr/bin/env Rscript

library(dplyr)

# Load input files
maf_file <- "./Mutect2_VCF0.8.maf"
cadd_file <- "./Mutect2_VCF0.8.maf_cadd_FULL.txt"

Loas_files3 <- read.delim(maf_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
Loas_files_CADD <- read.delim(cadd_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Extract sample pair ID
Loas_files_CADD$Paire <- gsub("Paire_(\\d+)_.*_GATK_somatic_filtered", "\\1", Loas_files_CADD$Tumor_Sample_Barcode)

# Create a unique identifier for each mutation
Loas_files3$full_ID <- paste(Loas_files3$Start_Position,
                             Loas_files3$Chromosome,
                             Loas_files3$End_Position,
                             Loas_files3$Reference_Allele,
                             Loas_files3$Tumor_Seq_Allele2,
                             sep = "_")

Loas_files_CADD$full_ID <- paste(Loas_files_CADD$Start_Position,
                                 Loas_files_CADD$Chromosome,
                                 Loas_files_CADD$End_Position,
                                 Loas_files_CADD$Reference_Allele,
                                 Loas_files_CADD$Tumor_Seq_Allele2,
                                 sep = "_")

# Select relevant columns from the CADD file
Loas_files_CADD_simple <- Loas_files_CADD[, c("full_ID", "Paire")]

# Merge CADD annotations into the main MAF file
Loas_files3_completed <- Loas_files3 %>%
  left_join(Loas_files_CADD_simple, by = "full_ID") %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "-", .))) %>%
  unique()

# Write the updated file
write.table(Loas_files3_completed,
            file = "CADD_Annotated_Mutect2_VCF0.8.maf",
            quote = FALSE,
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE)
