# Whole-exome pipeline analysis / CNV, SNP, and InDel extraction / Double-capture data visualization / Figure generation

## Whole-exome pipeline analysis
Raw FASTQ files were successively processed through deduplication using **Clumpify**, followed by **adapter trimming** and **quality filtering** with **BBDuk**. Read quality was then assessed using **FastQC**.

Filtered reads were mapped to the **hg38 reference genome** using **BWA-MEM**. **Duplicate reads** were marked with **Samtools**, and **base quality score recalibration** was performed using **GATK BQSR**. **Mapping depth** was assessed using **Mosdepth**.

All steps were performed by running `Fastq_Alignement.sh` script on the raw FASTQ directory that execute the following commandes:

```bash

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
	        # Step 4: Alignment with BWA MEM
        echo "Running BWA MEM for ${i}"
        bwa mem -t 35 \
            -R "@RG\tID:${i}\tSM:${i}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
            "$REFERENCE" "clean1_${i}1.fq.gz" "clean2_${i}2.fq.gz" \
            > "${i}_ALN.sam"

        # Step 5: Convert SAM to BAM sorted by query name (required for fixmate)
        echo "Sorting SAM by read name for ${i}"
        samtools view -@ 4 -Sb "${i}_ALN.sam" | samtools sort -n -@ 4 -o "${i}_querysort.bam"
        rm -f "${i}_ALN.sam"

        # Step 6: Run fixmate to add ms/mc tags
        echo "Running fixmate for ${i}"
        samtools fixmate -m -@ 4 "${i}_querysort.bam" "${i}_fixmate.bam"
        rm -f "${i}_querysort.bam"

        # Step 7: Sort BAM by position
        echo "Sorting fixmate BAM by position for ${i}"
        samtools sort -@ 4 -o "${i}_positionsort.bam" "${i}_fixmate.bam"
        rm -f "${i}_fixmate.bam"

        # Step 8: Mark duplicates and index final BAM
        echo "Running markdup for ${i}"
        samtools markdup -@ 4 "${i}_positionsort.bam" "Local_markdup_${i}_ALN.bam"
        samtools index -@ 4 "Local_markdup_${i}_ALN.bam"
        rm -f "${i}_positionsort.bam"
        
   
	NAME=$(echo $i | perl -pe 's/.+markdup_(.+)__ALN.bam/$1/');
    # Step 9: Depth coverage calculation with mosdepth
    if [ ! -f "${NAME}.regions.bed.gz" ]; then
        echo "Calculating depth coverage with mosdepth for ${i}"
       
        mosdepth -n --by $BED_FILE $NAME "Local_markdup_${i}_ALN.bam"
    else
        echo "Depth coverage for ${i} already calculated. Skipping step."
    fi

    # Step 10: FastQC on final BAM file
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

Next paire of Normal and cancer where analyzed by mutect2 gatk. For each sammple a job submission files where generated to run on HCPs on slurm environment.  

## CNVs calling
Copy variation has been achieved based on [GATK somatic copy number variation calling pipeline](https://gatk.broadinstitute.org/hc/en-us/articles/360035535892-Somatic-copy-number-variant-discovery-CNVs)

A Panel of normal as been first built by the following script commande : 
```json

{
  "CNVSomaticPanelWorkflow.normal_bams": ["/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_33CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_39CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_42CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_45CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_50CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_51CTR2__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_53CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_55CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_71CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_74CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_75CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_78CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_80CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_95CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_96CTR__ALN_sample_final.bam","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_97CTR3__ALN_sample_final.bam"

],
 "CNVSomaticPanelWorkflow.normal_bais": ["/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_33CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_39CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_42CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_45CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_50CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_51CTR2__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_53CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_55CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_71CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_74CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_75CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_78CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_80CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_95CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_96CTR__ALN_sample_final.bai","/home/florian2/Desktop/BAM_dataset/input/Filtred_BAMs/Local_markdup_97CTR3__ALN_sample_final.bai"],
 "CNVSomaticPanelWorkflow.pon_entity_id": "wes-do-gc",
  "CNVSomaticPanelWorkflow.ref_fasta_dict": "/home/florian2/Desktop/gatk-workflows/inputs/HOMO38_2.dict",
  "CNVSomaticPanelWorkflow.ref_fasta": "/home/florian2/Desktop/gatk-workflows/inputs/HOMO38_2.fa",
  "CNVSomaticPanelWorkflow.ref_fasta_fai": "/home/florian2/Desktop/gatk-workflows/inputs/HOMO38_2.fa.fai",
  "CNVSomaticPanelWorkflow.intervals": "/home/florian2/Desktop/gatk-workflows/inputs/targets_C.preprocessed.interval_list",
  "CNVSomaticPanelWorkflow.blacklist_intervals": "/home/florian2/Desktop/gatk-workflows/inputs/CNV_and_centromere_blacklist.hg38liftover.list",
  "CNVSomaticPanelWorkflow.PreprocessIntervals.bin_length":"0",

  "CNVSomaticPanelWorkflow.gatk_docker": "broadinstitute/gatk:4.6.1.0",
  
  "CNVSomaticPanelWorkflow.preemptible_attempts": "3"
}

```

## Double-capture pipeline data visualization

## Figure generation

