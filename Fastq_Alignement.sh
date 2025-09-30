#!/bin/bash
# ==============================
# Pipeline WES optimisé (streaming)
# ==============================

set -euo pipefail

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

# === PANEL OF NORMALS (CTR) ===
echo "=============================="
echo ">>> Panel of Normals Construction"
echo "=============================="

if command -v parallel >/dev/null 2>&1; then
  ls Local_markdup_*CTR*_ALN_final.bam 2>/dev/null | \
  parallel -j 4 "tpage --define threads={/.} ./Reformat_PON.tt | bash" \
    2> "logs/pon_parallel.log"
else
  for bam in Local_markdup_*CTR*_ALN_final.bam; do
    [ -e "$bam" ] || continue
    sample=${bam%.bam}
    echo "PON: $sample"
    tpage --define threads="${sample##*/}" ./Reformat_PON.tt | bash
  done
fi

vcf_list=$(ls Local_markdup_*CTR*__ALN_GATK.vcf.gz 2>/null | sed 's/^/-V /' | tr '\n' ' ')
if [ -n "$vcf_list" ]; then
  gatk GenomicsDBImport \
    -R "${REFERENCE}" \
    --genomicsdb-workspace-path "${WORKSPACE}" \
    $vcf_list \
    2> "logs/genomicsdbimport.log"
else
  echo "No CTR VCFs found for GenomicsDBImport."
fi
