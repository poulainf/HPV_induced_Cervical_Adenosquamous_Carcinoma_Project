# Whole-exome pipeline analysis / CNV, SNP, and InDel extraction / Double-capture data visualization / Figure generation
>The following pipeline complements the Materials and Methods section of the "..." research study.
>These scripts were designed to detect integrated HPV genome sequences from FFPE-fixed adenosquamous carcinoma samples.
   -Double-capture data visualization
>In addition, the cancer genome was explored through whole-exome sequencing (WES) to extract SNPs, InDels, and CNVs as part of the following analyses:
  -Whole-exome pipeline analysis
  -CNV, SNP, and InDel extraction
>In addition all scripts used for results vizualization are also reported


## Double-capture pipeline data visualization

## Whole-exome pipeline analysis
Raw FASTQ files were successively processed through deduplication using **Clumpify**, followed by **adapter trimming** and **quality filtering** with **BBDuk**. Read quality was then assessed using **FastQC**.

Filtered reads were mapped to the **hg38 reference genome** using **BWA-MEM**. **Duplicate reads** were marked with **Samtools**, and **base quality score recalibration** was performed using **GATK BQSR**. **Mapping depth** was assessed using **Mosdepth**.

All steps were performed by running `Fastq_Alignement.sh` script on the raw FASTQ directory that execute the following commandes:

```bash

# === CONFIG ===
THREADS=37
# THREADS2 n’est pas vraiment utilisé, on peut le retirer ou le garder pour GATK usage
THREADS2=23

BBMAP_ROOT="/home/florian2/bin/BBMap_39.01/bbmap"
BBDUK="${BBMAP_ROOT}/bbduk.sh"
CLUMPIFY="${BBMAP_ROOT}/clumpify.sh"

PICARD="java -jar /home/florian2/bin/picard.jar"
REFERENCE="./Homo_sapiens.GRCh38.chr.fa"
ADAPTORS="./Adaptor.fa"
BED_FILE="./exome_extended_130.bed.gz"
FASTQC="fastqc"
WORKSPACE="pon_db2"
SNPS_VCF="./1000G_phase1.snps.high_confidence.hg38.filtered.vcf.gz"
INDELS_VCF="./Mills_and_1000G_gold_standard.indels.hg38.filtered.vcf.gz"

BBDUK_XMX="8g"
CLUMPIFY_XMX="24g"

mkdir -p fastqc_reports logs tmp

# Préparer chrom_list.txt une seule fois (hors de la boucle)
if [ ! -f chrom_list.txt ]; then
  echo "Generating chrom_list.txt from BED file..."
  zcat "${BED_FILE}" | cut -f1 | sort -u > chrom_list.txt
fi

# === PIPELINE PAR ÉCHANTILLON ===
for r1 in *R1_001.fastq.gz; do
    [ -e "$r1" ] || { echo "No *R1_001.fastq.gz found. Exiting."; break; }

    i="${r1%%R1_001.fastq.gz}"
    r2="${i}R2_001.fastq.gz"

    echo "=============================="
    echo ">>> Starting sample: $i"
    echo "=============================="

    # ---------------------------
    # STEP 1: Clumpify + trimming streaming
    # ---------------------------
    if [ ! -f "clean1_${i}.fq.gz" ] || [ ! -f "clean2_${i}.fq.gz" ]; then
        echo "[Step 1] Clumpify + BBDuk streaming"
        "${CLUMPIFY}" -Xmx${CLUMPIFY_XMX} \
            in="${r1}" in2="${r2}" \
            dedupe=t \
            out=stdout.fq \
            t="${THREADS}" \
            2> "logs/${i}_clumpify.log" \
        | "${BBDUK}" -Xmx${BBDUK_XMX} \
            in=stdin.fq int=t \
            out1="clean1_${i}.fq.gz" out2="clean2_${i}.fq.gz" \
            ref="${ADAPTORS}" ktrim=r k=12 mink=10 hdist=1 tpe=t tbo=t minlength=50 maxlength=150 qtrim=rl trimq=20 \
            t="${THREADS}" \
            2> "logs/${i}_bbduk.log"
    else
        echo "[Step 1] clean1_${i}.fq.gz & clean2_${i}.fq.gz already exist. Skipping."
    fi

    # ---------------------------
    # STEP 2: FastQC on cleaned FASTQ
    # ---------------------------
    if [ ! -f "fastqc_reports/clean1_${i}_fastqc.html" ]; then
        echo "[Step 2] Running FastQC"
        "${FASTQC}" -t "${THREADS}" "clean1_${i}.fq.gz" "clean2_${i}.fq.gz" -o fastqc_reports/ \
            2> "logs/${i}_fastqc.log"
    else
        echo "[Step 2] FastQC already done for ${i}"
    fi

    # ---------------------------
    # STEP 3 → 4 combinés : BWA + markdup via samtools pipeline
    # ---------------------------
    if [ ! -f "Local_markdup_${i}_ALN.bam" ]; then
        echo "[Step 3–4] BWA MEM + fixmate + markdup (samtools)"

        # 3.1 Align to SAM
        echo "Running BWA MEM for ${i}"
        bwa mem -t "${THREADS}" \
            -R "@RG\tID:${i}\tSM:${i}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
            "$REFERENCE" "clean1_${i}.fq.gz" "clean2_${i}.fq.gz" \
            > "${i}_ALN.sam"

        # 3.2 Convert SAM -> BAM sorted by name
        echo "Sorting SAM by read name for ${i}"
        samtools view -@ 12 -Sb "${i}_ALN.sam" | samtools sort -n -@ 12 -o "${i}_querysort.bam"
        rm -f "${i}_ALN.sam"

        # 4.1 fixmate
        echo "Running fixmate for ${i}"
        samtools fixmate -m -@ 8 "${i}_querysort.bam" "${i}_fixmate.bam"
        rm -f "${i}_querysort.bam"

        # 4.2 Sort by position
        echo "Sorting fixmate BAM by position for ${i}"
        samtools sort -@ 12 -o "${i}_positionsort.bam" "${i}_fixmate.bam"
        rm -f "${i}_fixmate.bam"

        # 4.3 Mark duplicates + index
        echo "Running markdup for ${i}"
        samtools markdup -@ 8 "${i}_positionsort.bam" "Local_markdup_${i}_ALN.bam"
        samtools index -@ 8 "Local_markdup_${i}_ALN.bam"
        rm -f "${i}_positionsort.bam"

    else
        echo "[Step 3–4] Markdup BAM for ${i} already exists. Skipping."
    fi

    # ---------------------------
    # STEP 5: GATK BQSR (Parallel per chromosome)
    # ---------------------------
    if [ ! -f "Local_markdup_${i}_ALN_final.bam" ]; then
        echo "[Step 5] GATK BaseRecalibrator + ApplyBQSR"

        BAM="Local_markdup_${i}_ALN.bam"
        mkdir -p bqsr_tables

        cat chrom_list.txt | parallel -j ${THREADS} "
          echo 'Processing chromosome {} for ${i}'
          gatk BaseRecalibrator \
            -R ${REFERENCE} \
            -I ${BAM} \
            --known-sites ${SNPS_VCF} \
            --known-sites ${INDELS_VCF} \
            -L {} \
            -O bqsr_tables/${i}_recal_{}.table \
            2> logs/${i}_bqsr_chr_{}.log
        "

        recal_inputs=$(for c in $(cat chrom_list.txt); do echo -n "-I bqsr_tables/${i}_recal_${c}.table "; done)

        gatk GatherBQSRReports \
            ${recal_inputs} \
            -O "${i}_recal_data.table" \
            2> "logs/${i}_bqsr_merge.log"

        gatk ApplyBQSR \
            -R "${REFERENCE}" \
            -I "${BAM}" \
            --bqsr-recal-file "${i}_recal_data.table" \
            -O "Local_markdup_${i}_ALN_final.bam" \
            2> "logs/${i}_bqsr_apply.log"
    else
        echo "[Step 5] BQSR already done for ${i}"
    fi

    # ---------------------------
    # STEP 6: Coverage (mosdepth)
    # ---------------------------
    if [ ! -f "${i}.regions.bed.gz" ]; then
        echo "[Step 6] Running mosdepth"
        mosdepth -t "${THREADS}" -n --by "${BED_FILE}" "${i}" "Local_markdup_${i}_ALN_final.bam" \
            2> "logs/${i}_mosdepth.log"
    else
        echo "[Step 6] mosdepth already run for ${i}"
    fi

    echo ">>> Sample ${i} done."
done
```

## 🧬 SNP and INDEL Calling

SNP and INDEL calling was performed using the **Mutect2** pipeline from **GATK**.
A **Panel of Normals (PoN)** was created using **control FASTQ sequences**.

```bash

WORKSPACE="pon_db2"

# Step 1: Generate Reformat_PON_*.sh scripts from control BAM files
for bam in Local_markdup_*CTR*_ALN_sample_final.bam; do
    sample=$(basename "$bam" .bam)
    echo "Processing $sample"
    tpage --define threads="$sample" ./Reformat_PON.tt > "Reformat_PON_${sample}.sh"
done

chmod a+x *.sh

# Step 2: Run the generated scripts in parallel (limit to 4 concurrent jobs)
ls Reformat_PON_*.sh | xargs -n 1 -P 4 bash

# Step 3: Create a list of VCF files
vcf_list=$(ls Local_markdup_*CTR*__ALN_GATK.vcf.gz | sed 's/^/-V /' | tr '\n' ' ')

# Step 4: Run GenomicsDBImport with the list of VCFs
gatk GenomicsDBImport \
    -R "$REFERENCE" \
    --genomicsdb-workspace-path "$WORKSPACE" \
    $vcf_list

```

The Pon were next used with mutect2 command Next paire of normal and cancer where analyzed by mutect2 gatk. For each sammple a job submission files where generated to run on HCPs on slurm environment.  

## CNVs calling
Copy variation has been achieved based on [GATK somatic copy number variation calling pipeline](https://gatk.broadinstitute.org/hc/en-us/articles/360035535892-Somatic-copy-number-variant-discovery-CNVs)

A Panel of normal as been first built by the following script commande : 
```json

{
  "CNVSomaticPanelWorkflow.normal_bams": ["/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_33CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_39CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_42CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_45CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_50CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_51CTR2__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_53CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_55CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_71CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_74CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_75CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_78CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_80CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_95CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_96CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_97CTR3__ALN_sample_final.bam"],
 "CNVSomaticPanelWorkflow.normal_bais": ["/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_33CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_39CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_42CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_45CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_50CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_51CTR2__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_53CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_55CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_71CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_74CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_75CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_78CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_80CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_95CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_96CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_97CTR3__ALN_sample_final.bai"],
 "CNVSomaticPanelWorkflow.pon_entity_id": "wes-do-gc",
  "CNVSomaticPanelWorkflow.ref_fasta_dict": "/home/florian2/Desktop/gatk-workflows/inputs/Homo_sapiens.GRCh38.chr.dict",
  "CNVSomaticPanelWorkflow.ref_fasta": "/home/florian2/Desktop/gatk-workflows/inputs/Homo_sapiens.GRCh38.chr.fa",
  "CNVSomaticPanelWorkflow.ref_fasta_fai": "/home/florian2/Desktop/gatk-workflows/inputs/Homo_sapiens.GRCh38.chr.fa.fai",
  "CNVSomaticPanelWorkflow.intervals": "/home/florian2/Desktop/gatk-workflows/inputs/targets_C.preprocessed.interval_list",
  "CNVSomaticPanelWorkflow.blacklist_intervals": "/home/florian2/Desktop/gatk-workflows/inputs/CNV_and_centromere_blacklist.hg38liftover.list",
  "CNVSomaticPanelWorkflow.PreprocessIntervals.bin_length":"120",

  "CNVSomaticPanelWorkflow.gatk_docker": "broadinstitute/gatk:4.6.1.0",
  
  "CNVSomaticPanelWorkflow.preemptible_attempts": "3"
}

```

## Figure generation

Based on 

```R
#!/usr/bin/env Rscript

# Clean environment
rm(list = ls())
# Note: avoid explicit setwd() in scripts; consider passing paths as arguments or using here::here()
# setwd("/media/florian2/T7Shield1/Projet_Liege/")

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

# Clear any existing circular plot state
circos.clear()
if (dev.cur() > 1) dev.off()

# 1. Load and merge data
inte <- read.delim("./Inte.csv", header = TRUE, sep = "\t")
summary_df <- read.delim("./Resume_paris.csv", header = TRUE, sep = "\t")
combined_inte <- merge(summary_df, inte, by.x = "best_capture", by.y = "sample")

clinical <- read.delim("./Datas_cliniques.txt", header = TRUE, sep = "\t")
combined <- merge(combined_inte, clinical,
                  by.x = "paire", by.y = "Paire",
                  all = FALSE)

# Optionally inspect
# colnames(combined)

# Filter / annotate
combined_sel <- combined %>%
  group_by(paire, localisation, chr_position) %>%
  mutate(n_types = n_distinct(type_carcinome)) %>%
  ungroup()

# Note: corrected_Sample is not defined here; ensure it's loaded earlier
if (!exists("corrected_Sample")) {
  stop("corrected_Sample object not found — load it before subsetting.")
}
sample_ADC <- corrected_Sample %>% filter(type_carcinome == "ADC")
sample_SCC <- corrected_Sample %>% filter(type_carcinome != "ADC")

# 2. Read bed / depth / gene information
sample_insertions <- read.delim("./Sample_insertions.csv", header = TRUE, sep = ",")
depth_samples <- read.delim("./Full_hpv_coverage.txt", header = TRUE) %>%
  na.omit()

gene_infos <- read.delim("./Full_infos.gb", header = FALSE)

unique_pairs <- unique(depth_samples$New)
palette_colors <- paletteer_d("MoMAColors::Klein")

# Helper function: map to cytoband
get_cytoband <- function(chr, pos, cytoband_df) {
  hit <- cytoband_df %>%
    filter(V1 == chr, pos >= V2, pos < V3)
  if (nrow(hit) > 0) return(hit$V4[1])
  else return(NA_character_)
}

# 3. Loop over pairs for full-track circular plots
for (pair_id in unique_pairs) {
  
  # Determine HPV genome scale based on ADC data for this pair
  pos_ADC <- depth_samples %>%
    filter(New == pair_id, Cancer == "ADC") %>%
    pull(Position)
  if (length(pos_ADC) == 0) next
  HPV_size <- max(pos_ADC)
  
  # Build bed data for ADC and SCC tracks
  df_ADC <- depth_samples %>%
    filter(New == pair_id, Cancer == "ADC") %>%
    transmute(chr = "chrHPV",
              start = Position * (3.4e9 / HPV_size),
              stop  = start + (3.4e9 / HPV_size),
              value = Depht)
  df_SCC <- depth_samples %>%
    filter(New == pair_id, Cancer == "SCC") %>%
    transmute(chr = "chrHPV",
              start = Position * (3.4e9 / HPV_size),
              stop  = start + (3.4e9 / HPV_size),
              value = Depht)
  
  # Set column names
  colnames(df_ADC) <- c("chr", "start", "stop", "value")
  colnames(df_SCC) <- c("chr", "start", "stop", "value")
  
  # Cytoband: human + HPV
  human_cytoband <- read.cytoband(species = "hg38")$df
  HPV_cyto <- tibble::tibble(V1 = "chrHPV",
                             V2 = 0,
                             V3 = 3.4e9,
                             V4 = "q12",
                             stain = "Blue")
  colnames(HPV_cyto) <- colnames(human_cytoband)
  cytoband <- bind_rows(human_cytoband, HPV_cyto) %>%
    filter(V1 != "chrY")
  
  chromosome_index <- c(paste0("chr", c(1:22, "X")), "chrHPV")
  
  # Subset insertion samples
  seqid_ADC <- depth_samples %>%
    filter(New == pair_id, Cancer == "ADC") %>%
    pull(SeqId) %>%
    unique()
  seqid_SCC <- depth_samples %>%
    filter(New == pair_id, Cancer == "SCC") %>%
    pull(SeqId) %>%
    unique()
  
  sampleADC <- sample_insertions %>% filter(sample %in% seqid_ADC)
  sampleSCC <- sample_insertions %>% filter(sample %in% seqid_SCC)
  
  # Build linking beds
  bed1 <- sampleADC %>% transmute(chr = chr,
                                  start = chr_position,
                                  end   = chr_position + 2e6)
  bed2 <- sampleADC %>% transmute(chr = "chrHPV",
                                  start = position * (3.4e9 / HPV_size),
                                  end   = (position + 50) * (3.4e9 / HPV_size))
  bed3 <- sampleSCC %>% transmute(chr = chr,
                                  start = chr_position,
                                  end   = chr_position + 2e6)
  bed4 <- sampleSCC %>% transmute(chr = "chrHPV",
                                  start = position * (3.4e9 / HPV_size),
                                  end   = (position + 50) * (3.4e9 / HPV_size))
  
  bed5 <- tibble::tibble(chr = "chrHPV",
                         start = rep(1, nrow(bed2)),
                         end   = rep(3.4e9, nrow(bed2)))
  
  # Map gene annotations on HPV
  hpv_strain <- sampleADC %>% pull(Strain) %>% unique()
  hpv_genes <- gene_infos %>% filter(V1 == hpv_strain)
  HPV_cyto2 <- hpv_genes %>%
    transmute(V1 = V4,
              V2 = V2 * (3.4e9 / HPV_size),
              V3 = V3 * (3.4e9 / HPV_size),
              V4 = "q12",
              stain = palette_colors[1:nrow(hpv_genes)])
  colnames(HPV_cyto2) <- colnames(human_cytoband)
  
  # Prepare output file
  out_pdf <- paste0("circular_plot_pair_", pair_id, "_ADC.pdf")
  pdf(out_pdf)
  
  # Prepare alternating y-coordinates for annotation
  HPV_cyto2 <- HPV_cyto2 %>%
    arrange(V2) %>%
    mutate(row_id = row_number(),
           y1 = if_else(row_id %% 2 == 0, 0, 0.5),
           y2 = if_else(row_id %% 2 == 0, 0.5, 1)) %>%
    select(-row_id)
  
  # Initialize circular plot
  circos.clear()
  circos.par(start.degree = 90,
             gap.after = c(rep(1, 22), 4, 4))
  circos.initializeWithIdeogram(cytoband, plotType = NULL,
                               chromosome.index = chromosome_index)
  
  # ADC track
  circos.genomicTrackPlotRegion(df_ADC, ylim = c(0, max(df_ADC$value, na.rm = TRUE)),
    panel.fun = function(region, value, ...) {
      circos.genomicLines(region, value,
                          area = TRUE,
                          col = "#FFCCCC",
                          border = "#FF6666",
                          lwd = 0.5, ...)
    }, bg.border = NA)
  
  # SCC segments & track
  if (nrow(df_SCC) > 0) {
    circos.segments(bed4$start, 0, bed4$start, max(df_SCC$value, na.rm = TRUE))
    circos.genomicTrackPlotRegion(df_SCC, ylim = c(0, max(df_SCC$value, na.rm = TRUE)),
      panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value,
                            area = TRUE,
                            col = "#CCE5FF",
                            border = "#3399ff",
                            lwd = 0.5, ...)
      }, bg.border = NA)
    circos.segments(bed4$start, 0, bed4$start, max(df_SCC$value, na.rm = TRUE))
  }
  
  # Ideogram and text tracks
  circos.track(ylim = c(0, 0.001),
               panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2),
                             gsub(".*chr", "", CELL_META$sector.index),
                             cex = 0.6, niceFacing = TRUE, facing = "outside")
               }, track.height = mm_h(1), cell.padding = c(0, 0, 0, 0),
               bg.border = NA)
  circos.genomicIdeogram(cytoband = cytoband, track.height = 0.1)
  
  # Draw gene rectangles and labels
  for (i in seq_len(nrow(HPV_cyto2))) {
    circos.rect(sector.index = "chrHPV",
                xleft = HPV_cyto2$V2[i],
                ybottom = HPV_cyto2$y1[i],
                xright = HPV_cyto2$V3[i],
                ytop = HPV_cyto2$y2[i],
                col = HPV_cyto2$stain[i])
    circos.text((HPV_cyto2$V2[i] + HPV_cyto2$V3[i]) / 2,
                HPV_cyto2$y1[i] + 0.25,
                HPV_cyto2$V1[i],
                cex = 0.5, col = "white",
                facing = "inside", niceFacing = TRUE)
  }
  
  # Links between human and HPV
  circos.genomicLink(bed1, bed5, col = NA, border = "#696969")
  circos.genomicLink(bed5, bed5, col = "white", border = "white", lwd = 1.8)
  
  # Annotations
  text(0, 1.05, paste("Pair", pair_id), cex = 1.5)
  text(0.05, 0.84, "ADC", cex = 1, col = "#FF6666")
  text(0.05, 0.64, "SCC", cex = 1, col = "#3399ff")
  text(-0.2, 0.08, "Insertion:", cex = 1.2, col = "#696969")
  
  # Cytoband label counts
  bed1_unique <- bed1 %>% distinct(chr, start)
  bed1_unique <- bed1_unique %>%
    mutate(cytoband = mapply(get_cytoband, chr, start,
                             MoreArgs = list(cytoband_df = human_cytoband)))
  bed1_summary <- bed1_unique %>%
    group_by(cytoband) %>%
    summarise(n_sites = n(), .groups = "drop")
  for (j in seq_len(nrow(bed1_summary))) {
    y0 <- j - 1
    lab <- if (bed1_summary$n_sites[j] > 1) {
      sprintf("%d sites in %s locus", bed1_summary$n_sites[j], bed1_summary$cytoband[j])
    } else {
      sprintf("%d site in %s locus", bed1_summary$n_sites[j], bed1_summary$cytoband[j])
    }
    text(-0.2, -0.01 - (y0 * 0.06), lab, cex = 1, col = "#696969")
  }
  
  # HPV label
  text(-0.6, 0.85, gsub("_", "-", hpv_strain), cex = 1.25, col = "gray")
  
  # Finish
  circos.clear()
  dev.off()
}

# 4. (Optional) Small circular graphs & aggregated bed table
# ... (You can refactor the second big loop similarly) ...
```

