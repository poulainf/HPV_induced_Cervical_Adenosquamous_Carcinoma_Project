# Whole-exome pipeline analysis / CNV, SNP, and InDel extraction / Double-capture data visualization / Figure generation
>The following pipeline complements the Materials and Methods section of the "..." research study.
>These scripts were designed to detect integrated HPV genome sequences from FFPE-fixed adenosquamous carcinoma samples.
   -Double-capture data visualization
>In addition, the cancer genome was explored through whole-exome sequencing (WES) to extract SNPs, InDels, and CNVs as part of the following analyses:
  -Whole-exome pipeline analysis
  -CNV, SNP, and InDel extraction
>In addition all scripts used for results vizualization are also reported


## Double-capture pipeline data visualization
Fastq from HPV genome double capture hase been treated by [ViroCapt](https://github.com/maximewack/viroCapt). Briefly ViroCapt consist 
## Whole-exome pipeline analysis
Raw FASTQ files were successively processed through deduplication using **Clumpify**, followed by **adapter trimming** and **quality filtering** with **BBDuk**. Read quality was then assessed using **FastQC**.

Filtered reads were mapped to the **hg38 reference genome** using **BWA-MEM**. **Duplicate reads** were marked with **Samtools**, and **base quality score recalibration** was performed using **GATK BQSR**. **Mapping depth** was assessed using **Mosdepth**.

All steps were performed by running `Fastq_Alignement.sh` script on the raw FASTQ directory.

```bash
./Fastq_Alignement.sh
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


```bash
while read line;do

	CTR="$( echo $line | cut -f1 )";
	TEST="$( echo $line | cut -f2 )";
	PAIR="$( echo $line | cut -f3 )";
	TYPE="$( echo $TEST | perl -pe "s/\d\d//g" )"

	tpage --define CTR=$CTR --define TEST=$TEST --define PAIR=$PAIR --define TYPE=$TYPE RUN_mtect2_CLUSTER_3.0.tt  > Mutect2_RUN_${PAIR}_${TYPE}_CLUSTER.sh ;
done < Refs_samples2.txt


```


## CNVs calling
Copy variation has been achieved based on [GATK somatic copy number variation calling pipeline](https://gatk.broadinstitute.org/hc/en-us/articles/360035535892-Somatic-copy-number-variant-discovery-CNVs)

A Panel of normal as been first built by the following script commande : 
```bash

```

## Figure generation

Based on 

```R
```

