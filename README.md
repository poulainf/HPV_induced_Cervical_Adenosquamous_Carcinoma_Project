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

	tpage --define CTR=$CTR --define TEST=$TEST --define PAIR=$PAIR --define TYPE=$TYPE Run_mutect2_cluster.tt  > Mutect2_RUN_${PAIR}_${TYPE}_CLUSTER.sh ;
done < Refs_samples2.txt

```

VCF files produced by mutect2 SNV and InDels calling has been next analyze to filtrate mutations, annotated, produce maf, correct HUGO genes symbole. 


```bash
for i in ` ls Paire_*GATK_somatic_filtered.vcf ` ; do

    Trie_vcf_by_gnomad4.pl Selected_chr_0.01_gnomad.exomes.v4.1.vcf $i 0.05 mutect2 no snps

done

for i in ` ls Filtred_0.05*.vcf | sed -e s/".vcf"// ` ; do

	tpage --define CTR=$i ./annovar.tt > Annovar_RUN_${i}.sh;
 
done

chmod a+x *.sh

find . -type f | grep "Annovar" | parallel --tmuxpane '{}'

for i in ` ls Filtred_0.05_*_annotated.hg38_multianno.txt |sed -e s/'.txt'// `; do

	MyID="$( echo $i | perl -pe 's/.+(Paire_.+)_annotated.+/$1/g' )";
	echo $MyID ; annovar2maf.py -t ${MyID} ${i}.txt > ${i}.maf ;

done

for i in `ls *.maf ` ; do

	./maftools.R $i ;

done

cat Annovar_corrected_Filtred_0.05_._GATK_somatic_filtered_annotated.hg38_multianno_maftools.maf | head -n1 > MUTECT2_MAF.txt  

for i in `ls *.hg38_multianno_maftools.maf ` ; do

        SAMPLE="$( echo $i | perl -pe 's/.+Paire\_(\d+)\_(\w+)\_GATK.+/$1$2/g' )";
        NUM="$( echo $i | perl -pe 's/.+Paire\_(\d+)\_\w+\_GATK.+/$1/g' )";
        echo $SAMPLE
        echo $NUM
        tail -n +2 $i  >> MUTECT2_MAF.txt

done



```

Raw vcf files has been combined to ensure fishing of paired mutations. 

```bash

echo -n > Raw_VCF.vcf

for i in Paire_*_GATK_somatic_filtered.vcf; do

    SAMPLE="$(echo $i | perl -pe 's/(Paire_\d+_\w+)_GATK.+/$1/g')"
    grep -v "^#" "$i" | awk -v sample="$SAMPLE" '{ print sample "\t" $0 }' >> Raw_VCF.vcf

done

```


```bash

wget "https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz"

tail -n+2 gnomad.genomes.r4.0.indel.tsv >> whole_genome_SNVs.tsv

Annotate_MAF_CACC.pl Indels_vaf0.8.maf gnomad.genomes.r4.0.indel.tsv INDELs




FILE="whole_genome_SNVs.tsv"

# Nombre de splits
PARTS=35

# Nom du dossier temporaire
mkdir -p splits

# Extraire l'en-tête
head -n 1 "$FILE" > header.txt

# Compter le nombre de lignes (sans l �en-tête)
TOTAL=$(($(wc -l < "$FILE") - 1))
PER_FILE=$(( (TOTAL + PARTS - 1) / PARTS ))

# Split (en omettant la première ligne)
tail -n +2 "$FILE" | split -l "$PER_FILE" - splits/split_

# Ajouter l �en-tête Ã� chaque fichier
for f in splits/split_*; do
    cat header.txt "$f" > "$f.tsv"
    rm "$f"
done


INFILE="Mutect2_VCF0.8.maf"
SCRIPT="./Annotate_MAF_CACC.pl"
mkdir -p bash_jobs

i=0
for f in splits/*.tsv; do
    script_name="bash_jobs/run_$i.sh"
    echo "# > "$script_name"
    echo "$SCRIPT \"$INFILE\" \"$f\" \"_$i\"" >> "$script_name"
    chmod +x "$script_name"
    ((i++))
done

INFILE="Mutect2_VCF0.8.maf"
SCRIPT="./Annotate_MAF_CACC.pl"
mkdir -p bash_jobs

i=0
for f in splits/*.tsv; do
    script_name="bash_jobs/run_$i.sh"
--
rm tmp_fastqc_data.txt

# Calcul de la moyenne
if [ $sample_count -gt 0 ]; then
    avg_reads=$((total_reads / sample_count))
    echo "������� Moyenne de reads par échantillon : $avg_reads"
else
    echo "Aucun fichier FastQC trouvé."
fi


wget "https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz"
./Annotate_MAF_CACC.pl Mutect2_VCF0.8.maf gnomad.genomes.r4.0.indel.tsv INDELs



head -n1 Mutect2_VCF0.8.maf_cadd_ANNOTED_22.txt > Mutect2_VCF0.8.maf_cadd_FULL.txt
grep -vh "GenoCanyon_score" Mutect2_VCF0.8.maf_cadd_ANNOTED_*.txt >> Mutect2_VCF0.8.maf_cadd_FULL.txt
grep -v "GenoCanyon_score" Indels_vaf0.8.maf_cadd_ANNOTEDINDELs.txt >> Mutect2_VCF0.8.maf_cadd_FULL.txt


```


## CNVs calling
Copy variation has been achieved based on [GATK somatic copy number variation calling pipeline](https://gatk.broadinstitute.org/hc/en-us/articles/360035535892-Somatic-copy-number-variant-discovery-CNVs)

A Panel of normal as been first built by the use of noraml tissue BAM read files by following script commande: 
```bash

java -jar cromwell-47.jar run ./gatk4-somatic-cnvs/cnv_somatic_panel_workflow.wdl --inputs ./gatk4-somatic-cnvs/cnv_somatic_panel_workflow.b37.inputs.json


while read line;do\
\
	CTR="$( echo $line | cut -f1 )";\
	TEST="$( echo $line | cut -f2 )";\
	PAIR="$( echo $line | cut -f3 )";\
	TYPE="$( echo $TEST | perl -pe "s/\d\d//g" )"\
	echo $TYPE\
	tpage --define CTR=$CTR --define TEST=$TEST cnv_somatic_pair_workflow.inputs.tt > RUN_${PAIR}_${TYPE}_cnv_somatic_pair_workflow.b37.inputs;\
\
done < ./Refs_samples2.txt


for i in ` ls ./gatk4-somatic-cnvs/RUN_Paire*.inputs ` ; do  java -jar cromwell-47.jar run ./gatk4-somatic-cnvs/cnv_somatic_pair_workflow.wdl --inputs $i ; done
```

The produced 


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


## Figure generation

Based on 

```R

Circo_plot.R
Figures_papier_MAFs_18.R
Analyze_CNV8.2.R

```

