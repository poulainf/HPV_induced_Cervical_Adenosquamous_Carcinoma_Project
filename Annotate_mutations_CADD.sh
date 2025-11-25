#! /bin/bash

########################################
# CADD SNV + Indel annotation pipeline #
########################################

##############################
# 1. Download CADD resources #
##############################

# CADD whole-genome SNV scores (GRCh38)
wget "https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz"

# Decompress SNV file
gunzip whole_genome_SNVs.tsv.gz

# (Optionnel) Télécharger les scores CADD / gnomAD pour les indels
wget "https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz"
gunzip gnomad.genomes.r4.0.indel.tsv.gz

#########################################
# 2. Prepare CADD SNV file for splitting #
#########################################

FILE="whole_genome_SNVs.tsv"

# Number of splits
PARTS=35

# Temporary directory for split files
mkdir -p splits

# Extract header
head -n 1 "$FILE" > header.txt

# Count total number of lines (excluding header)
TOTAL_LINES=$(($(wc -l < "$FILE") - 1))
# Compute number of lines per split (ceiling division)
PER_FILE=$(( (TOTAL_LINES + PARTS - 1) / PARTS ))

# Split file (skip header)
tail -n +2 "$FILE" | split -l "$PER_FILE" - splits/split_

# Add header to each split file
for f in splits/split_*; do
    cat header.txt "$f" > "${f}.tsv"
    rm "$f"
done

##########################################
# 3. Generate bash jobs for SNV annotation
##########################################

INFILE="Mutect2_VCF0.8.maf"         # MAF file from Mutect2
SCRIPT="./Annotate_MAF_CACC.pl"     # Annotation script

mkdir -p bash_jobs

i=0
for f in splits/*.tsv; do
    script_name="bash_jobs/run_${i}.sh"
    {
        echo "#!/usr/bin/env bash"
        echo "# CADD annotation chunk $i"
        echo "$SCRIPT \"$INFILE\" \"$f\" \"_${i}\""
    } > "$script_name"
    chmod +x "$script_name"
    ((i++))
done

# Example: run all jobs in parallel (if GNU parallel is available)
# parallel ::: bash_jobs/run_*.sh

#####################################
# 4. Indel annotation with CADD data
#####################################

# Annotate indels using the CADD/gnomAD indel file
# (Indels_vaf0.8.maf should contain high-confidence indels)
./Annotate_MAF_CACC.pl Indels_vaf0.8.maf gnomad.genomes.r4.0.indel.tsv INDELs

###########################################################
# 5. Merge SNV and Indel CADD-annotated output into one MAF
###########################################################

# Start merged file with header from one SNV-annotated chunk
head -n 1 Mutect2_VCF0.8.maf_cadd_ANNOTED_22.txt > Mutect2_VCF0.8.maf_cadd_FULL.txt

# Append all SNV-annotated chunks, excluding repeated header lines
grep -vh "GenoCanyon_score" Mutect2_VCF0.8.maf_cadd_ANNOTED_*.txt >> Mutect2_VCF0.8.maf_cadd_FULL.txt

# Append CADD-annotated indels (also excluding header if present)
grep -v "GenoCanyon_score" Indels_vaf0.8.maf_cadd_ANNOTEDINDELs.txt >> Mutect2_VCF0.8.maf_cadd_FULL.txt
