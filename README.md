# Whole-Exome and Double-Capture Analysis Pipeline  
**CNV, SNP, and InDel extraction • HPV double-capture visualization • Figure generation**

This repository contains the computational pipelines described in the *Materials and Methods* section of the study "**[...]**".  
All scripts were designed to detect **integrated HPV genomic sequences** from FFPE-fixed adenosquamous carcinoma samples and to characterize the **somatic landscape** using whole-exome sequencing.

## Overview
- **HPV double-capture data visualization**
- **Whole-exome sequencing (WES) analysis pipeline**
- **SNP, InDel, and CNV extraction**
- **Figure generation scripts**

All visualization and post-processing scripts used to generate figures for the manuscript are included.

---

## 📊 Double-Capture Pipeline (HPV)

FASTQ files from HPV genome double-capture sequencing were processed using  
**[ViroCapt](https://github.com/maximewack/viroCapt)**.  
ViroCapt performs quality control, filtering, host-sequence removal, and viral-sequence enrichment to enable the detection of integrated HPV fragments.

---

## 🧬 Whole-Exome Sequencing Pipeline

Raw FASTQ files undergo the following steps:

1. **Deduplication** — *Clumpify*  
2. **Adapter trimming & quality filtering** — *BBDuk*  
3. **Quality control** — *FastQC*  
4. **Alignment to hg38** — *BWA-MEM*  
5. **Duplicate marking** — *Samtools*  
6. **Base Quality Score Recalibration (BQSR)** — *GATK*  
7. **Depth assessment** — *Mosdepth*

All steps are executed through:

```bash
./Fastq_Alignement.sh
```
---

## 🧬 SNP and InDel Calling (Mutect2)

### PoN construction
Somatic SNP and InDel calling was performed using GATK Mutect2.
A Panel of Normals (PoN) was first built from control BAM samples.

```bash
WORKSPACE="pon_db2"

# Step 1: Generate per-sample PON scripts
for bam in Local_markdup_*CTR*_ALN_sample_final.bam; do
    sample=$(basename "$bam" .bam)
    echo "Processing $sample"
    tpage --define threads="$sample" ./Reformat_PON.tt > "Reformat_PON_${sample}.sh"
done

chmod a+x *.sh

# Step 2: Run scripts in parallel
ls Reformat_PON_*.sh | xargs -n 1 -P 4 bash

# Step 3: Create VCF list
vcf_list=$(ls Local_markdup_*CTR*__ALN_GATK.vcf.gz | sed 's/^/-V /' | tr '\n' ' ')

# Step 4: Import into GenomicsDB
gatk GenomicsDBImport \
    -R "$REFERENCE" \
    --genomicsdb-workspace-path "$WORKSPACE" \
    $vcf_list
```

### Mutect2 paired tumor/normal calling
Job scripts were automatically generated for parallel submission.
Following GATK Best Practices, somatic variants were called using Mutect2 for each tumor sample by comparison to its matched normal, with Panel of Normals (PoN) correction applied. The resulting raw VCFs were then analyzed sequentially using LearnReadOrientationModel, GetPileupSummaries, CalculateContamination, and finally FilterMutectCalls.

```bash
while read line; do
    CTR="$(echo $line | cut -f1)"
    TEST="$(echo $line | cut -f2)"
    PAIR="$(echo $line | cut -f3)"
    TYPE="$(echo $TEST | perl -pe "s/\d\d//g")"

    tpage --define CTR=$CTR --define TEST=$TEST --define PAIR=$PAIR --define TYPE=$TYPE \
        Run_mutect2_cluster.tt > Mutect2_RUN_${PAIR}_${TYPE}_CLUSTER.sh ;
done < Refs_samples2.txt

find . -type f | grep "Mutect2_RUN_" | parallel --tmuxpane '{}'
```
### Variant filtering, mutation annotation and MAF generation.
Variants were filtered by selecting those with Mutect2 filters 'PASS', 'clustered_events', or 'haplotype', tumor alternative allele depth >30, matched normal depth >20, VAF_T >0.08, VAF_N <0.04, AD_T_ALT >3, and TLOD >10. Germline variants were excluded using a gnomAD (v4) allele frequency threshold of <0.01. Filtered VCFs were annotated with ANNOVAR using refGene and avsnp150. Non-exonic mutations were defined as those within ±120 bp of exon boundaries. Mutations were considered shared between ADC and SCC from the same patient when matching in chromosomal position, mutation type, and alleles. To enhance detection of shared events, a rescue step retained unfiltered mutations if also present in filtered calls from the paired sample. Intra-pair comparisons (ADC vs. SCC) were performed within tumors; inter-pair comparisons estimated a cross-pair identity threshold. Putative deleterious variants were extracted based on CADD v1.6 Phred scores >20. 

#### Mutect2 filter and gnomad selection 
```bash
for i in `ls Paire_*GATK_somatic_filtered.vcf`; do
    Trie_vcf_by_gnomad4.pl Selected_chr_0.01_gnomad.exomes.v4.1.vcf $i 0.08 mutect2 no snps
done

for i in `ls Filtred_0.08*.vcf | sed -e s/".vcf"//`; do
    tpage --define CTR=$i ./annovar.tt > Annovar_RUN_${i}.sh
done

chmod a+x *.sh
find . -type f | grep "Annovar" | parallel --tmuxpane '{}'
```
#### Convert annotated variants to MAF:
```bash
for i in `ls Filtred_0.08_*_annotated.hg38_multianno.txt | sed -e s/'.txt'//`; do
    MyID="$(echo $i | perl -pe 's/.+(Paire_.+)_annotated.+/$1/g')"
    annovar2maf.py -t ${MyID} ${i}.txt > ${i}.maf
done

for i in `ls *.maf`; do
    ./maftools.R $i
done
```
#### Build a combined MAF file:
```bash
cat Annovar_corrected_Filtred_0.08_._GATK_somatic_filtered_annotated.hg38_multianno_maftools.maf | \
    head -n1 > MUTECT2_MAF.txt  

for i in `ls *.hg38_multianno_maftools.maf`; do
    SAMPLE="$(echo $i | perl -pe 's/.+Paire\_(\d+)\_(\w+)\_GATK.+/$1$2/g')"
    tail -n +2 $i >> MUTECT2_MAF.txt
done
```
#### Mutation sub selection and CADD Annotation (SNVs & Indels)
```bash
awk -v OFS='\t' '{start=$2-120; end=$3+120; if(start<0) start=0; print $1, start, end}' ./Twist_Comprehensive_Exome_Covered_Targets_hg38.bed > exome_extended_120.bed

SNPs_INDELs_mutation_subselection.R

Annotate_mutations_CADD.sh
```
---

## 🧬 CNV Calling (GATK CNV)
Somatic copy number variants (CNVs) were called using the GATK Somatic CNV pipeline (version 4.6.1.0). A PoN generated from 126 matched normal tissue samples was used for copy number recalibration in each tumor-normal pair. CNVs were defined as regions with a log2 copy number ratio less than -0.3 (deletions) or greater than 0.3 (amplifications). Genic copy number alterations were extracted using FuncotateSegments (GATK). Focal CNV events (both deletions and amplifications) were selected if their size was less than one million base pairs. Comparison between ADC and SCC segments was based on matching start and stop gene loci. The proportion of shared segments was tested between ADC and SCC components from the same patient. Stringent criteria were applied to defined gene amplification with a log2 fold change <-1  and deletion with a fold change >0.6.   

#### Build the Panel of Normals

```bash
java -jar cromwell-47.jar run \
    ./gatk4-somatic-cnvs/cnv_somatic_panel_workflow.wdl \
    --inputs ./gatk4-somatic-cnvs/cnv_somatic_panel_workflow.b37.inputs.json
```

#### Generate inputs for paired analyses

```bash
while read line; do
    CTR="$(echo $line | cut -f1)"
    TEST="$(echo $line | cut -f2)"
    PAIR="$(echo $line | cut -f3)"
    TYPE="$(echo $TEST | perl -pe "s/\d\d//g")"

    tpage --define CTR=$CTR --define TEST=$TEST \
        cnv_somatic_pair_workflow.inputs.tt > RUN_${PAIR}_${TYPE}_cnv_somatic_pair_workflow.b37.inputs
done < Refs_samples2.txt
```

#### Run CNV workflow

```bash
for i in ./gatk4-somatic-cnvs/RUN_Paire*.inputs; do
    java -jar cromwell-47.jar run \
        ./gatk4-somatic-cnvs/cnv_somatic_pair_workflow.wdl \
        --inputs $i
done
```
#### CNV merging and annotation
Merging of denoised CR and called CNV segments.

```bash
while read line;do   

    CTR="$( echo $line | cut -f1 )";
    TEST="$( echo $line | cut -f2 )";
    PAIR="$( echo $line | cut -f3 )";
    TYPE="$( echo $TEST | perl -pe 's/\d\d//g' )"
	echo $TYPE

	if [ ! -f Point_CR_${PAIR}.txt ]; then
		grep "CONTIG" Local_markdup_${CTR}__ALN_final.denoisedCR.tsv|awk -F '\t', '{print $0, "\t","Type" ,"\t","Paire" }' > Point_CR_${PAIR}.txt
	fi
	
	if [ ! -f Point_CR_COMBINED.txt ]; then\
		grep "CONTIG" Local_markdup_${CTR}__ALN_final.denoisedCR.tsv | awk -F '\t', '{print $0, "\t","Type" ,"\t","Paire" }' > Point_CR_COMBINED.txt
	fi
   

	grep -v "@" Local_markdup_${CTR}__ALN_final.denoisedCR.tsv|grep -v "CONTIG" |  awk -F '\t' -v pair="$PAIR" '{print $0, "\tCTR\t" pair}' >> Point_CR_${PAIR}.txt
	grep -v "@" Local_markdup_${TEST}__ALN_final.denoisedCR.tsv|grep -v "CONTIG" |  awk -F '\t' -v pair="$PAIR" -v type="$TYPE" '{print $0, "\t" type "\t" pair}'  >> Point_CR_${PAIR}.txt
	
	grep -v "@" Local_markdup_${CTR}__ALN_final.denoisedCR.tsv|grep -v "CONTIG" |  awk -F '\t' -v pair="$PAIR" '{print $0, "\tCTR\t" pair}' >> Point_CR_COMBINED.txt
	grep -v "@" Local_markdup_${TEST}__ALN_final.denoisedCR.tsv|grep -v "CONTIG" |  awk -F '\t' -v pair="$PAIR" -v type="$TYPE" '{print $0, "\t" type "\t" pair}'  >> Point_CR_COMBINED.txt

done < Refs_samples2.txt

grep "CONTIG" Point_CR_COMBINED.txt > tmp_Point_CR_COMBINED.txt
grep -v "CONTIG" Point_CR_COMBINED.txt | sort | uniq >> tmp_Point_CR_COMBINED.txt
mv tmp_Point_CR_COMBINED.txt Point_CR_COMBINED.txt




while read line;do   
    
    CTR="$( echo $line | cut -f1 )";
    TEST="$( echo $line | cut -f2 )";
    PAIR="$( echo $line | cut -f3 )";
    TYPE="$( echo $TEST | perl -pe "s/\d\d//g" )"
	echo $TEST

	if [ ! -f Test_${PAIR}.txt ]; then

		grep "CONTIG" Local_markdup_${CTR}__ALN_final.called.seg|awk -F '\t', '{print $0, "\t","Type" ,"\t","Paire" }' > Test_${PAIR}.txt
	fi
	
	if [ ! -f Test_COMBINED.txt ]; then
		grep "CONTIG" Local_markdup_${CTR}__ALN_final.called.seg|awk -F '\t', '{print $0, "\t","Type" ,"\t","Paire" }' > Test_COMBINED.txt

    fi
	grep -v "@" Local_markdup_${CTR}__ALN_final.called.seg|grep -v "CONTIG" |  awk -F '\t' -v pair="$PAIR" '{print $0, "\tCTR\t" pair}' >> Test_${PAIR}.txt
	grep -v "@" Local_markdup_${TEST}__ALN_final.called.seg|grep -v "CONTIG" |  awk -F '\t' -v pair="$PAIR" -v type="$TYPE" '{print $0, "\t" type "\t" pair}'  >> Test_${PAIR}.txt
	
    grep -v "@" Local_markdup_${CTR}__ALN_final.called.seg|grep -v "CONTIG" |  awk -F '\t' -v pair="$PAIR" '{print $0, "\tCTR\t" pair}' >> Test_COMBINED.txt
	grep -v "@" Local_markdup_${TEST}__ALN_final.called.seg|grep -v "CONTIG" |  awk -F '\t' -v pair="$PAIR" -v type="$TYPE" '{print $0, "\t" type "\t" pair}'  >> Test_COMBINED.txt
        
done < Refs_samples2.txt

grep  "alt_allele" Local_markdup_33ADC__ALN_final.called.seg.funcotated.tsv | awk -F '\t', '{print $0, "\t","TYPE" ,"\t","PAIRE" }' > combined_files.txt
rm combined_files.txt

while read -r line; do
    CTR="$(echo "$line" | cut -f1)"
    TEST="$(echo "$line" | cut -f2)"
    PAIR="$(echo "$line" | cut -f3)"
    TYPE="$(echo "$TEST" | perl -pe 's/\d\d//g')"

    for SAMPLE in "$CTR" "$TEST"; do
        grep -v "alt_allele" "Local_markdup_${SAMPLE}__ALN_final.called.seg.funcotated.tsv" \
        
            | awk -F '\t' -v type="$TYPE" -v pair="$PAIR" '{
                split($6, genes, ",");  # split column 6 by comma
                $6 = length(genes);     # replace column 6 with number of genes
                print $0 "\t" type "\t" pair
            }' >> combined_files.txt
    done
done < ./Refs_samples2.txt

```
---

## 📈 Figure Generation

Figures were produced using the following R scripts:
```
Circo_plot.R
Figures_papier_MAFs_18.R
Analyze_CNV8.2.R
```

