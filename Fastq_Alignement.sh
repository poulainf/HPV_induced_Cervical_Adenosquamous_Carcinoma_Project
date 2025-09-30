#!/bin/bash

# Stop the script on errors and handle unset variables
set -euo pipefail

# Define paths to tools and reference files
BBDUK="/home/florian2/bin/BBMap_39.01/bbmap/bbduk.sh"
CLUMPIFY="/home/florian2/bin/bbmap/clumpify.sh"
PICARD="java -jar /home/florian2/bin/picard.jar"
REFERENCE="./Homo_sapiens.GRCh38.chr.fa"
ADAPTORS="./Adaptor.fa"
BED_FILE="./Twist_Comprehensive_Exome_Covered_Targets_hg38.bed.gz"
FASTQC="fastqc" # Path to FastQC (ensure it's installed and in your PATH)

# Create FastQC output directory if it doesn't already exist
mkdir -p ./fastqc_reports/

# Process each sample
for i in $(ls *R1_001.fastq* | grep -v "md5" | grep -v "dedup" | sed -e 's/R1_001.fastq.gz//' -e 's/R1_001.fastq//'); do
    echo "################### Start ${i} #############################"

	if [ ! -f "Local_markdup_${i}_ALN.bam" ]; then
	
		if [ ! -f "clean1_${i}1.fq.gz" ] &&  [ ! -f "clean2_${i}2.fq.gz" ]; then
		# Step 1: Decompress input FASTQ files
		if [ ! -f "${i}R1_001.fastq" ] && [ ! -f "dedup2_${i}R2_001.fastq.gz" ]; then
			echo "Decompressing ${i}R1_001.fastq.gz and ${i}R2_001.fastq.gz"
			gunzip -k "${i}R1_001.fastq.gz" "${i}R2_001.fastq.gz"
		else
			echo "Decompressed FASTQ files for ${i} already exist. Skipping step."
		fi

		# Step 2: Deduplication with Clumpify
		if [ ! -f "dedup1_${i}R1_001.fastq" ] && [ ! -f "dedup2_${i}R2_001.fastq.gz" ]; then
			echo "Running Clumpify for ${i}"
			$CLUMPIFY in="${i}R1_001.fastq" in2="${i}R2_001.fastq" \
				out="dedup1_${i}R1_001.fastq" out2="dedup2_${i}R2_001.fastq" dedupe
		else
			echo "Deduplicated files for ${i} already exist. Skipping step."
		fi

		# Step 2b: Remove decompressed FASTQ to save space
		if [ -f "${i}R1_001.fastq" ] && [ -f "${i}R2_001.fastq" ]; then
			echo "Removing raw decompressed FASTQ files for ${i} to save space"
			rm -f "${i}R1_001.fastq" "${i}R2_001.fastq"
		fi
		
		
		# Step 3: Adapter trimming and quality filtering with BBDuk
		if [ ! -f "clean1_${i}1.fq" ] &&  [ ! -f "clean2_${i}2.fq.gz" ]; then
			if [  -f "dedup1_${i}R1_001.fastq.gz" ] && [ ! -f "dedup2_${i}R2_001.fastq.gz" ]; then
				echo "Decompressing dedup1_${i}R1_001.fastq.gz and dedup2_${i}R2_001.fastq.gz"
				gunzip -k "dedup1_${i}R1_001.fastq.gz" "dedup2_${i}R2_001.fastq.gz"
			fi
			
			echo "Running BBDuk for ${i}"
			$BBDUK in1="dedup1_${i}R1_001.fastq" in2="dedup2_${i}R2_001.fastq" \
				out1="clean1_${i}1.fq" out2="clean2_${i}2.fq" \
				ref=$ADAPTORS ktrim=r k=12 mink=10 hdist=1 tpe tbo minlength=50 maxlength=150 qtrim=rl trimq=20
				
				
			if [  -f "dedup1_${i}R1_001.fastq.gz" ] && [ ! -f "dedup2_${i}R2_001.fastq.gz" ]; then
				echo "rm dedup1_${i}R1_001.fastq.gz and dedup2_${i}R2_001.fastq.gz"
				rm -f "dedup1_${i}R1_001.fastq" "dedup2_${i}R2_001.fastq"
			fi
		else
			echo "Cleaned FASTQ files for ${i} already exist. Skipping step."
		fi

		# Compress deduplicated and cleaned FASTQ files
		for file in "dedup1_${i}R1_001.fastq" "dedup2_${i}R2_001.fastq" "clean1_${i}1.fq" "clean2_${i}2.fq"; do
			if [ -f "$file" ]; then
				echo "Compressing $file"
				gzip "$file"
			else
				echo "File $file is already compressed or missing. Skipping compression."
			fi
		done
		fi
		# Step 4: Alignment with BWA (adding read group tags dynamically)

			echo "Running BWA MEM for ${i}"
			sample_name=$(echo "$i" | sed 's/.+markdup_\(.*\)__ALN.bam/\1/')
			bwa mem -t 35 \
				-R "@RG\tID:${sample_name}\tSM:${sample_name}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
				$REFERENCE "clean1_${i}1.fq.gz" "clean2_${i}2.fq.gz" > "${i}_ALN.sam"


		# Step 5: Convert SAM to sorted BAM
		if [ ! -f "${i}_ALN.bam" ]; then
			echo "Sorting SAM to BAM for ${i}"
			samtools sort -@ 8 -o "${i}_ALN.bam" "${i}_ALN.sam"
			rm -f "${i}_ALN.sam"
		else
			echo "Sorted BAM file for ${i} already exists. Skipping step."
		fi
 fi
		# Step 6: Mark duplicates with Picard
		if [ ! -f "Local_markdup_${i}_ALN.bam" ]; then
			echo "Marking duplicates for ${i}"
			$PICARD MarkDuplicates \
				INPUT="${i}_ALN.bam" \
				OUTPUT="Local_markdup_${i}_ALN.bam" \
				METRICS_FILE="metrics_${i}.txt" \
				CREATE_INDEX=true
			rm -f "${i}_ALN.bam"
		else
			echo "Duplicate-marked BAM file for ${i} already exists. Skipping step."
		fi
   
 NAME=$(echo $i | perl -pe 's/.+markdup_(.+)__ALN.bam/$1/');
    # Step 7: Depth coverage calculation with mosdepth
    if [ ! -f "${NAME}.regions.bed.gz" ]; then
        echo "Calculating depth coverage with mosdepth for ${i}"
       
        mosdepth -n --by $BED_FILE $NAME "Local_markdup_${i}_ALN.bam"
    else
        echo "Depth coverage for ${i} already calculated. Skipping step."
    fi

    # Step 8: FastQC on final BAM file
    if [ ! -f "./fastqc_reports/Local_markdup_${i}_ALN_fastqc.html" ]; then
        echo "Running FastQC for: Local_markdup_${i}_ALN.bam"
        $FASTQC "Local_markdup_${i}_ALN.bam" -o ./fastqc_reports/
    else
        echo "FastQC for ${i} already completed. Skipping step."
    fi

    echo "Sample ${i} processed successfully."
    
    
    # Step 9: FastQC on final BAM file
    
   if [ ! -f "Local_markdup_${i}_ALN_sample_final.bam" ]; then
   gatk BaseRecalibrator -I Local_markdup_${i}_ALN.bam -R $REFERENCE \
	--known-sites 1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
	  -O Local_markdup_${i}_recal_data.table
	  
	gatk ApplyBQSR -R $REFERENCE -I Local_markdup_${i}_ALN.bam --bqsr-recal-file Local_markdup_${i}_recal_data.table -O Local_markdup_${i}_ALN_sample_final.bam
    else
		echo "FastQC for ${i} already completed. Skipping step."
    fi
    echo "#######################################################"
done

