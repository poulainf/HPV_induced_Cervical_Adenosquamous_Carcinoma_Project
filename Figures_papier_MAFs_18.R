#!/usr/bin/Rscript
rm(list = ls()) # clean up the environment
setwd("/media/florian2/T7Shield1/Projet_Liege/")
setwd("/media/florian/T7Shield/Projet_Liege/")
#setwd("~/Projet_Liege/")
# library(reshape)
# library(ggpubr)
library("ggplot2")
library("reshape2")
library("grid")
library(ggplot2)
library(randomcoloR)
library(magrittr)
library(ggpubr)
# library(plotly)
# library(ggrepel)
# library(wesanderson)
library(dplyr)
library(RColorBrewer)
library(Formula)
library(lattice)
# library(survival)
library(Hmisc)
library(ggsignif)
require(dplyr)
require(forcats)
library(lsa)
library(scales)
library(gglogo)
library(car)
library(ggplot2)
library(ggpmisc)
library(tidyr)
library(stringr)
library("patchwork")
library(ggh4x)
library(data.table)
library(pheatmap)
library(tidyr)
library(patchwork)
library(maftools)
library(gridExtra)
library(clusterProfiler)
library(org.Hs.eg.db)
library(KEGGREST)
library(HGNChelper)
library(patchwork)


# 
# 

raw_files2 <- fread("/media/florian2/T7/To_send/VCF_Muect2/Raw_VCF.vcf", sep = "\t", header = FALSE, nThread = 4)
raw_files2$Paire<-as.numeric(gsub("Paire_(\\d+)_.+","\\1",raw_files2$V1 ))

##########################################  SNPs ######################################################
Loas_files <- read.delim(file = "/media/florian2/T7/To_send/VCF_Muect2/test_SNPs_MAF.txt", header = TRUE, sep = "\t")
Loas_files$Tumor_Sample_Barcode<-gsub("CIN3","SCC",Loas_files$Tumor_Sample_Barcode)
Loas_files$Tumor_Sample_Barcode<-gsub("SCC2","SCC",Loas_files$Tumor_Sample_Barcode)
Loas_files$Tumor_Sample_Barcode<-gsub("ADC2","ADC",Loas_files$Tumor_Sample_Barcode)
Loas_files$Paire<-gsub("Paire_(\\d+)_.+_GATK_somatic_filtered","\\1",Loas_files$Tumor_Sample_Barcode)
Loas_files$Cancer<-gsub("Paire_\\d+_(.+)_GATK_somatic_filtered","\\1",Loas_files$Tumor_Sample_Barcode)
Loas_files$Tumor_Sample_Barcode<-gsub("_GATK_somatic_filtered","",Loas_files$Tumor_Sample_Barcode)
Loas_files$Cancer[which(Loas_files$Cancer=="SCC3")]<-"SCC"
Loas_files$Cancer[which(Loas_files$Cancer=="CIN3")]<-"SCC"
Loas_files$Cancer[which(Loas_files$Cancer=="SCC2")]<-"SCC"
Loas_files$Cancer[which(Loas_files$Cancer=="ADC2")]<-"ADC"



Loas_files <- Loas_files[!(Loas_files$Paire %in% c(5, 7, 9,10)), ]
datas_cliniques<-read.delim(file = "./Datas_cliniques.txt", header = T, sep = "\t")

Loas_files_introns<-Loas_files[which(Loas_files$Variant_Classification==''|Loas_files$Variant_Classification=="Unknown"),]

Exon_files<-read.delim(file = "Twist_Comprehensive_Exome_Covered_Targets_hg38.bed", header = F, sep = "\t")
Exon_files$up<-Exon_files$V2-120
Exon_files$down<-Exon_files$V2+120
#Exon_files$chr<-paste0("chr",Exon_files$V1)
Exon_files$chr<-Exon_files$V1

# Convert to data.tables
snps <- as.data.table(Loas_files_introns)
exons <- as.data.table(Exon_files[,c(4,5,6)])

# Convert to data.table
setDT(snps)
setDT(exons)

# Rename exon columns for overlap matching
setnames(exons, c("chr", "up", "down"), c("Chromosome", "start", "end"))
# SNPs as point intervals
snps[, start := Start_Position]
snps[, end := Start_Position]
# Set keys for fast overlap
setkey(exons, Chromosome, start, end)
setkey(snps, Chromosome, start, end)
# Find overlaps (SNPs inside any exon)
snps2 <- data.frame(foverlaps(snps, exons, nomatch = 0))
snps2 <- snps2[, c(intersect(colnames(Loas_files), colnames(snps2)))]
snps2$Variant_Classification<-"Non exonic"
Loas_filesDD<-Loas_files[-which(Loas_files$Variant_Classification==''|Loas_files$Variant_Classification=="Unknown"),]
Loas_files<-rbind(Loas_filesDD,snps2)
##########################################  SNPs ######################################################


Loas_files$Depht_tumor<-as.numeric(gsub("[^:]+:[^:]+:[^:]+:([^:,]+):.+","\\1",Loas_files$Otherinfo13))
Loas_files$Depht_normal<-as.numeric(gsub("[^:]+:[^:]+:[^:]+:([^:,]+):.+","\\1",Loas_files$Otherinfo14))
Loas_files$TLOD<-as.numeric(gsub(".+TLOD\\=(\\d*\\.*\\d*).*","\\1",Loas_files$Otherinfo11))
Loas_files <- Loas_files[Loas_files$Otherinfo10 %in% c("PASS", "clustered_events", "haplotype"), ]
Loas_files2<-Loas_files[-which(Loas_files$Depht_tumor<30),]
Loas_files2<-Loas_files2[-which(Loas_files2$Depht_normal<20),]
Loas_files3<-Loas_files2[-which(Loas_files2$TLOD<10),]
Loas_files3$VAF_T<-as.numeric(gsub("[^\\:]+\\:[^\\:\\,]+[^\\:]*\\:([^\\:\\,]+).*\\:.+","\\1",Loas_files3$Otherinfo13))
Loas_files3$VAF_N<-as.numeric(gsub("[^\\:]+\\:[^\\:\\,]+[^\\:]*\\:([^\\:\\,]+).*\\:.+","\\1",Loas_files3$Otherinfo14))
Loas_files3<-Loas_files3[-which(Loas_files3$VAF_T<0.08),]
Loas_files3 <- Loas_files3[which(Loas_files3$VAF_N < 0.04),]

Loas_files3$AD_T_ALT <- as.numeric(gsub("^[^:]+:[^,]+,([^:]+):.*", "\\1", Loas_files3$Otherinfo13))
Loas_files3 <- Loas_files3[which(Loas_files3$AD_T_ALT >= 4),]


raw_files2$V1<-gsub("CIN3","SCC",raw_files2$V1)
raw_files2$V1[raw_files2$V1=="Paire_2_ADC"]<-"Paire_18_ADC"
raw_files2$V1<-gsub("CIN3","SCC",raw_files2$V1)
raw_files2$V1<-gsub("SCC2","SCC",raw_files2$V1)
raw_files2$V1<-gsub("ADC2","ADC",raw_files2$V1)
raw_files2$Cancer<-gsub("Paire_\\d+_(.+)","\\1",raw_files2$V1 )
raw_files2$combi <- paste(raw_files2$V1,raw_files2$V2, raw_files2$V3, raw_files2$V5, raw_files2$V6)
Loas_files3$combi <- paste(paste(sep = "_","Paire",Loas_files3$Paire,Loas_files3$Cancer), Loas_files3$Chromosome, Loas_files3$Start_Position, Loas_files3$Reference_Allele, Loas_files3$Tumor_Seq_Allele2)
Loas_files3$combi2 <- paste(paste(sep = "_","Paire",Loas_files3$Paire), Loas_files3$Chromosome, Loas_files3$Start_Position, Loas_files3$Reference_Allele, Loas_files3$Tumor_Seq_Allele2)


Loas_files3 <- Loas_files3[(Loas_files3$Paire %in% c("11","12","13","14","16","15","18","1","3","4","6","8" )), ]

Loas_files3_repeche <- data.frame()
for (i in 1:length(unique(Loas_files3$Paire))){
  my_sample1a<-Loas_files3[which(Loas_files3$Paire==unique(Loas_files3$Paire)[i]&Loas_files3$Cancer=="ADC"),]
  my_sample2a<-Loas_files3[which(Loas_files3$Paire==unique(Loas_files3$Paire)[i]&Loas_files3$Cancer=="SCC"),]
  
  my_sample1a$research<-gsub("ADC","SCC",my_sample1a$combi)
  my_sample2a$research<-gsub("SCC","ADC",my_sample2a$combi)
  
  my_repeche1<-my_sample1a[my_sample1a$research%in% raw_files2$combi,]
  my_repeche2<-my_sample2a[my_sample2a$research%in% raw_files2$combi,]
  
  if (nrow(my_repeche1) > 0) {
    match_idx1 <- match(my_repeche1$research, raw_files2$combi)
    my_repeche1$Tumor_Sample_Barcode <- gsub("ADC", "SCC", my_repeche1$Tumor_Sample_Barcode)
    my_repeche1$Cancer <- gsub("ADC", "SCC", my_repeche1$Cancer)
    
    my_repeche1$Otherinfo11 <- raw_files2[[9]][match_idx1]
    my_repeche1$Otherinfo13 <- raw_files2[[11]][match_idx1]
    my_repeche1$Otherinfo14 <- raw_files2[[12]][match_idx1]
    my_repeche1$Paire <- raw_files2[[13]][match_idx1]
    my_repeche1$Tumor_Sample_Barcode<- raw_files2[[1]][match_idx1]
    my_repeche1$research <- NULL
  }
  
  
  if (nrow(my_repeche2) > 0) {
    match_idx2 <- match(my_repeche2$research, raw_files2$combi)
    
    my_repeche2$Tumor_Sample_Barcode <- gsub("SCC", "ADC", my_repeche2$Tumor_Sample_Barcode)
    my_repeche2$Cancer <- gsub("SCC", "ADC", my_repeche2$Cancer)
    
    my_repeche2$Otherinfo11 <- raw_files2[[9]][match_idx2]
    my_repeche2$Otherinfo13 <- raw_files2[[11]][match_idx2]
    my_repeche2$Otherinfo14 <- raw_files2[[12]][match_idx2]
    my_repeche2$Paire <- raw_files2[[13]][match_idx2]
    my_repeche2$Tumor_Sample_Barcode<- raw_files2[[1]][match_idx2]
    
    my_repeche2$research <- NULL
  }
  
  Loas_files3_repeche <- rbind( Loas_files3_repeche,my_repeche1, my_repeche2)
}

Loas_files3_repeche<-unique(Loas_files3_repeche)
Loas_files3_repeche$Fishing<-"Yes"
Loas_files3$Fishing<-"No"

Loas_files3<-rbind(Loas_files3,Loas_files3_repeche)
Loas_files3<-unique(Loas_files3)


Loas_files3$Depht_tumor<-as.numeric(gsub("[^:]+:[^:]+:[^:]+:([^:,]+):.+","\\1",Loas_files3$Otherinfo13))
Loas_files3$Depht_normal<-as.numeric(gsub("[^:]+:[^:]+:[^:]+:([^:,]+):.+","\\1",Loas_files3$Otherinfo14))
Loas_files3$TLOD<-as.numeric(gsub(".+TLOD\\=(\\d*\\.*\\d*).*","\\1",Loas_files3$Otherinfo11))
Loas_files3$VAF_T<-as.numeric(gsub("[^\\:]+\\:[^\\:\\,]+[^\\:]*\\:([^\\:\\,]+).*\\:.+","\\1",Loas_files3$Otherinfo13))
Loas_files3$VAF_N<-as.numeric(gsub("[^\\:]+\\:[^\\:\\,]+[^\\:]*\\:([^\\:\\,]+).*\\:.+","\\1",Loas_files3$Otherinfo14))
Loas_files3$AD_T_ALT <- as.numeric(gsub("^[^:]+:[^,]+,([^:]+):.*", "\\1", Loas_files3$Otherinfo13))


Loas_files3$combi3 <- paste(Loas_files3$Chromosome, Loas_files3$Start_Position)
Loas_files3<-Loas_files3[!Loas_files3$Reference_Allele=="GCGGCCGCCGCCGCCGCCGCTGCGGGCGGCGCGCACCAGAACTCGGCCGTGGCGGCGGCGGCGGCGGCG",]
datas_cliniques<-read.delim(file = "./Datas_cliniques.txt", header = T, sep = "\t")
Loas_files3<-merge(x=datas_cliniques, y= Loas_files3, by = "Paire" )

Loas_files3<-unique(Loas_files3)

#

Loas_files3 <- read.delim(file = "./Mutect2_VCF0.8.maf", header = TRUE, sep = "\t")
Loas_files_CADD <- read.delim(file = "./Mutect2_VCF0.8.maf_cadd_FULL.txt", header = TRUE, sep = "\t")

Loas_files_CADD$Paire<-gsub("Paire_(\\d+)_.+_GATK_somatic_filtered","\\1",Loas_files_CADD$Tumor_Sample_Barcode)

Loas_files3$full_ID<-paste(sep = "_",Loas_files3$Start_Position,Loas_files3$Chromosome,Loas_files3$End_Position,
                           Loas_files3$Reference_Allele,Loas_files3$Tumor_Seq_Allele2)

Loas_files_CADD$full_ID<-paste(sep = "_",Loas_files_CADD$Start_Position,Loas_files_CADD$Chromosome,Loas_files_CADD$End_Position,
                               Loas_files_CADD$Reference_Allele,Loas_files_CADD$Tumor_Seq_Allele2)

Loas_files_CADD_simple<-Loas_files_CADD[,c(208,209)]

Loas_files3_completed <- Loas_files3 %>%
  left_join(Loas_files_CADD_simple, by = "full_ID") %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "-", .)))

Loas_files3_completed<-unique(Loas_files3_completed)

write.table(Loas_files3_completed,file = "Mutect2_VCF0.8.maf",quote = F,sep="\t",col.names = T,row.names = F)



Loas_files_INDELs_mutect2<-Loas_files3_completed[Loas_files3_completed$Variant_Type=="INS"|Loas_files3_completed$Variant_Type=="DEL",]
Loas_files_SNPs<-Loas_files3_completed[Loas_files3_completed$Variant_Type=="SNP"|Loas_files3_completed$Variant_Type=="TNP"|Loas_files3_completed$Variant_Type=="DNP",]


####################
## Snps 
Loas_files_SNPs$MUT_TYPE<-paste0(Loas_files_SNPs$Reference_Allele,">",Loas_files_SNPs$Tumor_Seq_Allele2)

summarize_SNPs1 <- Loas_files_SNPs %>%
  group_by(New, Cancer, Variant_Classification) %>%
  mutate(Count = n())%>%
  ungroup() %>%
  dplyr::select(New, Cancer, Variant_Classification, Count) %>%
  distinct() %>%   
  
  complete(
    New = 1:12,
    Cancer  = c("ADC", "SCC"),
    Variant_Classification = c(
      "Missense_Mutation", "Silent", "Non exonic",
      "Nonsense_Mutation", "Nonstop_Mutation", "Translation_Start_Site"
    ),
    fill = list(Count = 0)
  )


summarize_SNPs2 <- Loas_files_SNPs %>%
  group_by(New, Chromosome, Start_Position, MUT_TYPE) %>%
  mutate(commun = n_distinct(Cancer)) %>%
  ungroup() %>%
  mutate(Cancer = ifelse(commun > 1, "Commun", Cancer)) %>%
  group_by(New, Cancer, Variant_Classification) %>%
  mutate(Count = (n()/2))%>%
  ungroup() %>%
  filter(Cancer=="Commun")%>%
  dplyr::select(New, Cancer, Variant_Classification, Count) %>%
  distinct() %>%   
  
  complete(
    New = 1:12,
    Cancer  = "Commun",
    Variant_Classification = c(
      "Missense_Mutation", "Silent", "Non exonic",
      "Nonsense_Mutation", "Nonstop_Mutation", "Translation_Start_Site"
    ),
    fill = list(Count = 0)
  )

summarize_SNPs<-rbind(summarize_SNPs1,summarize_SNPs2)

summarize_SNPs$ordre<-1
summarize_SNPs$ordre[summarize_SNPs$Cancer=="SCC"]<-2
summarize_SNPs$ordre[summarize_SNPs$Cancer=="Commun"]<-3

summarize_SNPs <- summarize_SNPs %>%
  arrange(New, ordre, Variant_Classification)
summarize_SNPs$ordre<-NULL

colnames(summarize_SNPs)<-c("Paire","Cancer","Variant Classification","SNPs")

write.table(x = summarize_SNPs,file = "15SNPs_summary.csv",quote = F,sep = ",",row.names = F,col.names = T)

#########

Full_summarize_SNPs <- Loas_files_SNPs %>%
  group_by(New, Chromosome, Start_Position, MUT_TYPE, Variant_Classification) %>%
  mutate(commun = n_distinct(Cancer)) %>%
  ungroup() %>%
  mutate(Clonal = ifelse(commun > 1, "Yes", "No")) %>%
  dplyr::select(New, Cancer,Clonal,Variant_Type,Variant_Classification,Chromosome,Start_Position,End_Position,Reference_Allele,Tumor_Seq_Allele2,Hugo_Symbol,Transcript_ID,
                Depht_normal,Depht_tumor,VAF_N,VAF_T,TLOD,CADD_PHRED_1.7,Fishing)

test<-Loas_files_SNPs[as.numeric(Loas_files_SNPs$Depht_tumor)<14,]

colnames(Full_summarize_SNPs)<-c("Paire","Cancer","Clonal","Variant Type","Variant Classification","Chromosome","Start","End","Reference Allele","Tumor Allele","Hugo Symbol","Transcript ID","Depht normal",
                                 "Depht tumor","VAF Normal",'VAF Tumor',"TLOD","CADD_PHRED_1.7","Fishing")
nrow(Full_summarize_SNPs)
Full_summarize_SNPs<-unique(Full_summarize_SNPs)
write.table(x = Full_summarize_SNPs,file = "15Full_SNPs_summary.csv",quote = F,sep = ";",row.names = F,col.names = T)

##########

Loas_files_SNPs$MUT_TYPE<-paste0(Loas_files_SNPs$Reference_Allele,">",Loas_files_SNPs$Tumor_Seq_Allele2)

Loas_files_SNPs$MUT_TYPE<-paste0(Loas_files_SNPs$Reference_Allele,">",Loas_files_SNPs$Tumor_Seq_Allele2)

Loas_files_SNPs$MUT_TYPE[grepl("...>...",Loas_files_SNPs$MUT_TYPE)]



Loas_files_SNPs$MUT_TYPE[grepl("^..>..$",Loas_files_SNPs$MUT_TYPE)]


colnames(Loas_files_SNPs)
summarize_SNPs_type <- Loas_files_SNPs %>%
  group_by(New, Chromosome, Start_Position, MUT_TYPE) %>%
  mutate(commun = n_distinct(Cancer)) %>%
  ungroup() %>%
  mutate(Cancer = ifelse(commun > 1, "Commun", Cancer)) %>%
  group_by(Cancer, MUT_TYPE) %>%
  summarise(Count = n())%>%
  ungroup() %>%
  filter(grepl("^.>.$",MUT_TYPE))%>%
  group_by(Cancer) %>%
  mutate(tot = sum(Count))%>%
  ungroup()%>%
  group_by(Cancer, MUT_TYPE) %>%
  mutate(Prop = (sum(Count)/tot)*100)%>%
  ungroup()%>%
  dplyr::select(Cancer, MUT_TYPE,Prop) %>%
  complete(
    Cancer  = c("ADC", "SCC","Commun"),
    MUT_TYPE = c("A>C", "A>T", "A>G","C>A","C>T","C>G","G>A","G>C","G>T","T>A","T>C","T>G"),
    fill = list(Count = 0,Prop=0)
  )



summarize_SNPs_type$ordre<-1
summarize_SNPs_type$ordre[summarize_SNPs_type$Cancer=="SCC"]<-2
summarize_SNPs_type$ordre[summarize_SNPs_type$Cancer=="Commun"]<-3

summarize_SNPs_type <- summarize_SNPs_type %>%
  arrange(ordre, MUT_TYPE)


# Plot
x1<-ggplot(data = summarize_SNPs_type, aes(y = MUT_TYPE, x = Prop, fill = MUT_TYPE)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "Proportion", y = "Mutation Type") +
  geom_vline(xintercept = 0, size = 2) +
  scale_x_continuous(expand = c(0, 0)) + 
  theme_classic() +   # ✅ garde les axes et ticks visibles
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.text.y = element_text(size = 18, face = "bold", color = "black",
                               margin = margin(r = 0)),
    axis.title.y = element_blank(),
    
    # ✅ Axe X et petits ticks
    axis.line.x = element_line(color = "black", linewidth = 1),
    axis.ticks.x = element_line(color = "black", linewidth = 1),
    axis.ticks.length.x = unit(0.3, "cm"),
    
    # ✅ Taille du texte des facettes
    strip.background = element_blank(),
    strip.text = element_text(size = 20, face = "bold", color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.45, size = 20)
  ) +
  facet_grid(. ~ reorder(Cancer,ordre))   # attention: reorder Cancer en amont si besoin

# facet_wrap(reorder(Cancer,ordre)~.,ncol=1, scale="free")

ggsave("15SNPs_compo.pdf", plot = x1, width = 5, height = 5, limitsize = FALSE, device = 'pdf', dpi = 300)

dev.off()
####################
## InDels 
Loas_files_INDELs__INSE<-Loas_files_INDELs_mutect2[Loas_files_INDELs_mutect2$Variant_Type=="INS",]
Loas_files_INDELs__INSE$LEN<-nchar(Loas_files_INDELs__INSE$Tumor_Seq_Allele2)
Loas_files_INDELs__INSE$pos<-"+"
Loas_files_INDELs__DEL<-Loas_files_INDELs_mutect2[Loas_files_INDELs_mutect2$Variant_Type=="DEL",]
Loas_files_INDELs__DEL$LEN<-nchar(Loas_files_INDELs__DEL$Reference_Allele)
Loas_files_INDELs__DEL$pos<-"-"
Loas_files_INDELs__combi<-rbind(Loas_files_INDELs__INSE,Loas_files_INDELs__DEL)
Loas_files_INDELs__combi$Variant_Type

summarize_INDELs1 <- Loas_files_INDELs__combi %>%
  group_by(New, Cancer, Variant_Type) %>%
  mutate(Count = n())%>%
  ungroup() %>%
  dplyr::select(New, Cancer, Variant_Type, Count) %>%
  distinct() %>%   
  
  complete(
    New = 1:12,
    Cancer  = c("ADC", "SCC"),
    Variant_Type = c(
      "INS", "DEL"
    ),
    fill = list(Count = 0)
  )


summarize_INDELs2 <- Loas_files_INDELs__combi %>%
  group_by(New, Chromosome, Start_Position,Reference_Allele, Tumor_Seq_Allele2) %>%
  mutate(commun = n_distinct(Cancer)) %>%
  ungroup() %>%
  mutate(Cancer = ifelse(commun > 1, "Commun", Cancer)) %>%
  group_by(New, Cancer, Variant_Type) %>%
  mutate(Count = (n()/2))%>%
  ungroup() %>%
  filter(Cancer=="Commun")%>%
  dplyr::select(New, Cancer, Variant_Type, Count) %>%
  distinct() %>%   
  
  complete(
    New = 1:12,
    Cancer  = "Commun",
    Variant_Type = c(
      "INS", "DEL"
    ),
    fill = list(Count = 0)
  )

summarize_INDELs<-rbind(summarize_INDELs1,summarize_INDELs2)

summarize_INDELs$ordre<-1
summarize_INDELs$ordre[summarize_INDELs$Cancer=="SCC"]<-2
summarize_INDELs$ordre[summarize_INDELs$Cancer=="Commun"]<-3

summarize_INDELs <- summarize_INDELs %>%
  arrange(New, ordre, Variant_Type)
summarize_INDELs$ordre<-NULL
colnames(summarize_INDELs)<-c("Paire","Cancer","Variant Type","SNPs")

write.table(x = summarize_INDELs,file = "15INDELs_summary.csv",quote = F,sep = ",",row.names = F,col.names = T)



#################

## InDels 

Loas_files_INDELs_mutect2<-Loas_files3_completed[Loas_files3_completed$Variant_Type=="INS"|Loas_files3_completed$Variant_Type=="DEL",]

Loas_files_INDELs__INSE<-Loas_files_INDELs_mutect2[Loas_files_INDELs_mutect2$Variant_Type=="INS",]
Loas_files_INDELs__INSE$LEN<-nchar(Loas_files_INDELs__INSE$Tumor_Seq_Allele2)
Loas_files_INDELs__INSE$pos<-"+"
Loas_files_INDELs__DEL<-Loas_files_INDELs_mutect2[Loas_files_INDELs_mutect2$Variant_Type=="DEL",]
Loas_files_INDELs__DEL$LEN<-nchar(Loas_files_INDELs__DEL$Reference_Allele)
Loas_files_INDELs__DEL$pos<-"-"
Loas_files_INDELs__combi<-rbind(Loas_files_INDELs__INSE,Loas_files_INDELs__DEL)

Loas_files_INDELs__combi$Variant_Classification[Loas_files_INDELs__combi$Variant_Classification==("Frame_Shift_Del")|
                                                  Loas_files_INDELs__combi$Variant_Classification==("Frame_Shift_Ins")]<-"Frame_Shift"

Loas_files_INDELs__combi$Variant_Classification[Loas_files_INDELs__combi$Variant_Classification==("Inframe_INDEL")|
                                                  Loas_files_INDELs__combi$Variant_Classification==("In_Frame_Del")|
                                                  Loas_files_INDELs__combi$Variant_Classification==("In_Frame_Ins")]<-"Inframe"

summarize_INDELs1 <- Loas_files_INDELs__combi %>%
  group_by(New, Cancer, Variant_Classification) %>%
  mutate(Count = n())%>%
  ungroup() %>%
  dplyr::select(New, Cancer, Variant_Classification, Count) %>%
  distinct() %>%   
  complete(
    New = 1:12,
    Cancer = c("ADC","SCC"),
    Variant_Classification  = c("Non exonic","Frame_Shift","Nonsense_Mutation","Inframe","Translation_Start_Site","Nonstop_Mutation"),
    fill = list(Count = 0)
  )


summarize_INDELs2 <- Loas_files_INDELs__combi %>%
  group_by(New, Chromosome, Start_Position,Reference_Allele, Tumor_Seq_Allele2) %>%
  mutate(commun = n_distinct(Cancer)) %>%
  ungroup() %>%
  mutate(Cancer = ifelse(commun > 1, "Commun", Cancer)) %>%
  group_by(New, Cancer, Variant_Classification) %>%
  mutate(Count = (n()/2))%>%
  ungroup() %>%
  filter(Cancer=="Commun")%>%
  dplyr::select(New, Cancer, Variant_Classification, Count) %>%
  distinct() %>%   
  
  complete(
    New = 1:12,
    Cancer  = "Commun",
    Variant_Classification  = c("Non exonic","Frame_Shift","Nonsense_Mutation","Inframe","Translation_Start_Site","Nonstop_Mutation"),
    fill = list(Count = 0)
  )

summarize_INDELs<-rbind(summarize_INDELs1,summarize_INDELs2)

summarize_INDELs$ordre<-1
summarize_INDELs$ordre[summarize_INDELs$Cancer=="SCC"]<-2
summarize_INDELs$ordre[summarize_INDELs$Cancer=="Commun"]<-3

summarize_INDELs <- summarize_INDELs %>%
  arrange(New, ordre, Variant_Classification)
summarize_INDELs$ordre<-NULL
colnames(summarize_INDELs)<-c("Paire","Cancer","Variant Classification","SNPs")

write.table(x = summarize_INDELs,file = "15INDELs_summary_variant.csv",quote = F,sep = ",",row.names = F,col.names = T)


#########

Loas_files_INDELs_mutect2<-Loas_files3_completed[Loas_files3_completed$Variant_Type=="INS"|Loas_files3_completed$Variant_Type=="DEL",]
Full_summarize_INDELs <- Loas_files_INDELs__combi %>%
  group_by(New, Chromosome, Start_Position,Reference_Allele, Tumor_Seq_Allele2) %>%
  mutate(commun = n_distinct(Cancer)) %>%
  ungroup() %>%
  mutate(Clonal = ifelse(commun > 1, "Yes", "No")) %>%
  dplyr::select(New, Cancer,Clonal,Variant_Type,Variant_Classification,Chromosome,Start_Position,End_Position,Reference_Allele,Tumor_Seq_Allele2,Hugo_Symbol,Transcript_ID,
                Depht_normal,Depht_tumor,VAF_N,VAF_T,TLOD,LEN,CADD_PHRED_1.7,Fishing)



colnames(Full_summarize_INDELs)<-c("Paire","Cancer","Clonal","Variant Type","Variant Classification","Chromosome","Start","End","Reference Allele","Tumor Allele","Hugo Symbol","Transcript ID","Depht normal",
                                   "Depht tumor","VAF Normal",'VAF Tumor',"TLOD","InDels size","CADD_PHRED_1.7","Fishing")

write.table(x = Full_summarize_INDELs,file = "15Full_INDELs_summary.csv",quote = F,sep = ";",row.names = F,col.names = T)

###########

Loas_files_INDELs__INSE<-Loas_files_INDELs_mutect2[Loas_files_INDELs_mutect2$Variant_Type=="INS",]
Loas_files_INDELs__INSE$LEN<-nchar(Loas_files_INDELs__INSE$Tumor_Seq_Allele2)
Loas_files_INDELs__INSE$pos<-"+"
Loas_files_INDELs__DEL<-Loas_files_INDELs_mutect2[Loas_files_INDELs_mutect2$Variant_Type=="DEL",]
Loas_files_INDELs__DEL$LEN<-nchar(Loas_files_INDELs__DEL$Reference_Allele)
Loas_files_INDELs__DEL$pos<-"-"
Loas_files_INDELs__combi2<-rbind(Loas_files_INDELs__INSE,Loas_files_INDELs__DEL)
#Loas_files_INDELs__combi2<-Loas_files_INDELs__combi[,c(10,11,196,197)]

Loas_files_INDELs__combi2$LEN[Loas_files_INDELs__combi2$LEN>6]<-6

colnames(Loas_files_INDELs__combi2)

summarize_INDELs <- Loas_files_INDELs__combi2 %>%
  group_by(New, Chromosome, Start_Position,Reference_Allele, Tumor_Seq_Allele2) %>%
  mutate(commun = n_distinct(Cancer)) %>%
  ungroup() %>%
  mutate(Cancer = ifelse(commun > 1, "Commun", Cancer)) %>%
  group_by(New, Cancer, Variant_Classification) %>%
  mutate(Count = (n()/2))%>%
  
  group_by(LEN, pos,Cancer) %>%
  summarise(count = n()) %>%
  ungroup()%>%
  group_by(Cancer) %>%
  mutate(Prop = (count / sum(count)) * 100)

colnames(summarize_INDELs)
summarize_INDELs$y<-paste0(summarize_INDELs$pos,summarize_INDELs$LEN)




summarize_INDELs$ordre<-1
summarize_INDELs$ordre[summarize_INDELs$Cancer=="SCC"]<-2
summarize_INDELs$ordre[summarize_INDELs$Cancer=="Commun"]<-3

summarize_INDELs <- summarize_INDELs %>%
  arrange(ordre, MUT_TYPE)



# Plot
x1 <- ggplot(data = summarize_INDELs, aes(y = Prop, x = reorder(y,as.numeric(y)), fill = pos)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  
  scale_fill_manual(
    values = c(
      "+" = "#8AC78F",
      "-" = "#C389C9"
    )
  ) + 
  
  labs( x = "Proportion", y = "Mutation Type") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) + 
  geom_hline(yintercept = 0,col="black",size=1.5)+
  
  theme(legend.position = "none")+
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 16, face = "bold", color = "black"),
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.text.y = element_text(size = 16, face = "bold", color = "black",
                               margin = margin(r = 0)),
    axis.title.y = element_blank(),
    
    # ✅ Axe X et petits ticks
    axis.line.x = element_blank(),
    axis.ticks.x = element_line(color = "black", linewidth = 1),
    #axis.ticks.length.x = unit(0.3, "cm"),
    
    # ✅ Taille du texte des facettes
    strip.background = element_blank(),
    strip.text = element_text(size = 20, face = "bold", color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.45, size = 20)
    
    
  )+
  facet_wrap(reorder(Cancer,ordre)~.,ncol=1, scale="free_y")


summarize_INDELs2 <- Loas_files_INDELs__combi2 %>%
  group_by(New, Chromosome, Start_Position,Reference_Allele, Tumor_Seq_Allele2) %>%
  mutate(commun = n_distinct(Cancer)) %>%
  ungroup() %>%
  mutate(Cancer = ifelse(commun > 1, "Commun", Cancer)) %>%
  group_by(New, Cancer, Variant_Classification) %>%
  mutate(Count = (n()/2))%>%
  
  group_by(Variant_Type,Cancer) %>%
  summarise(count = n()) %>%
  ungroup()%>%
  group_by(Cancer) %>%
  mutate(Prop = (count / sum(count)) * 100)

summarize_INDELs2$ordre<-1
summarize_INDELs2$ordre[summarize_INDELs2$Cancer=="SCC"]<-2
summarize_INDELs2$ordre[summarize_INDELs2$Cancer=="Commun"]<-3

summarize_INDELs2 <- summarize_INDELs2 %>%
  arrange(ordre)



# Plot
x2 <- ggplot(data = summarize_INDELs2, aes(y = Prop, x = reorder(Cancer,as.numeric(ordre)), fill = Variant_Type)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  scale_fill_manual(
    values = c(
      "INS" = "#9ACBB3",
      "DEL" = "#C89CCB"
    )
  ) + 
  labs( x = "Proportion", y = "Mutation Type") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) + 
  theme(legend.position = "none")+
  geom_hline(yintercept = 0,col="black",size=1.5)+
  # theme(
  #   axis.text.x = element_text(size = 18, face = "bold", color = "black", vjust = 1),
  #   axis.title.x = element_text(size = 20, face = "bold", color = "black", vjust = 1),
  #   axis.text.y   = element_text(size = 18, face = "bold", color = "black"),
  #   axis.title.y = element_blank(),
  #   panel.background = element_blank(),
  #   # axis.title = element_blank(),
  #   # axis.ticks.x = element_blank(),
  #   # axis.ticks.y = element_blank(),
  #   # plot.margin = unit(c(0, 1, 1, 1), "cm"),
  #   # strip.text  = element_blank(),
  #   plot.title = element_text(face = "bold", hjust = 0.45, size = 20)
  
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title.x =  element_blank(),
    axis.text.y = element_text(size = 16, face = "bold", color = "black",
                               margin = margin(r = 0)),
    axis.title.y = element_blank(),
    
    
    # ✅ Axe X et petits ticks
    axis.line.x = element_blank(),
    axis.ticks.x =  element_blank(),
    #axis.ticks.length.x = unit(0.3, "cm"),
    
    # ✅ Taille du texte des facettes
    strip.background = element_blank(),
    strip.text = element_text(size = 20, face = "bold", color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.45, size = 20)
    
    
  )+
  facet_wrap(reorder(Cancer,ordre)~.,ncol=1, scale="free")






bottom_row <- x2+x1+
  plot_layout(ncol = 2, guides = "collect")+ plot_layout(widths = c(1, 10))


bottom_row


ggsave("15INDELs_compo.pdf", plot = bottom_row, width = 12, height = 8, limitsize = FALSE, device = 'pdf', dpi = 300)

dev.off()



##################
dirs <- c("MAFs_SCC", "MAFs_ADC", "MAFs_COMMUN")

# Optional: Check existence before deleting
for (d in dirs) {
  if (dir.exists(d)) {
    message("Deleting directory: ", d)
    system(paste("rm -rf", d))
  }
}

# Recreate directories
system("mkdir MAFs_SCC MAFs_ADC MAFs_COMMUN")




##################################################################

variant_colors <- c(
  "Frame_Shift" = "#F1E417",
  "Missense_Mutation" = "#F1E417",
  "Translation_Start_Site" = "#C29C33",
  "Nonsense_Mutation" = "#C85EC7",
  "Nonstop_Mutation" = "#ED2024",
  "Inframe" = "#32C581",
  "Silent" = "#32C581",
  "Non exonic" = "#F05C3B"
)



Result_SNPs<- data.frame(Paire = integer(),
                         ADC = integer(),
                         SCC = integer(),
                         class = numeric(),
                         New = numeric(),
                         Subtype = integer(),
                         Integration = integer(),
                         Cancer = integer(),
                         Variant_Classification = integer(),
                         prop_commun = numeric(),
                         count_variant = numeric(),
                         statut = integer(),
                         compare = integer(),
                         stringsAsFactors = FALSE)

#Loas_files_SNPs<-Loas_files_SNPs[-which(Loas_files_SNPs$Variant_Classification=="",)]
#           Loas_files_SNPs$Variant_Classification=="Unknown"|Loas_files_SNPs$Variant_Classification=="Silent")]<-"Non-coding"
# Loas_files_SNPs$Variant_Classification
unique(Loas_files_INDELs_mutect2$Variant_Classification)

datas_cliniques<-read.delim(file = "./Datas_cliniques.txt", header = T, sep = "\t")



#Loas_files_SNPs$Cancer<-gsub("_crossed","",Loas_files_SNPs$Cancer)
unique(Loas_files_SNPs$Paire)





# 
# 
# lamlCommun <- maftools::read.maf(maf = Loas_files_SNPs[Loas_files_SNPs$New==12,], vc_nonSyn = unique(Loas_files_SNPs$Variant_Classification))
# 
# rainfallPlot(maf = lamlCommun, detectChangePoints = TRUE, pointSize = 0.4)
# 
# 
# 
# 
# lamlCommun <- maftools::read.maf(maf = Loas_files_SNPs[Loas_files_SNPs$New==4&Loas_files_SNPs$Cancer=="ADC",], 
#                                  vc_nonSyn = unique(Loas_files_SNPs$Variant_Classification))
# rainfallPlot(maf = lamlCommun, detectChangePoints = TRUE, pointSize = 0.4)
# lamlCommun <- maftools::read.maf(maf = Loas_files_SNPs[Loas_files_SNPs$New==1&Loas_files_SNPs$Cancer=="SCC",], 
#                                  vc_nonSyn = unique(Loas_files_SNPs$Variant_Classification))
# rainfallPlot(maf = lamlCommun, detectChangePoints = TRUE, pointSize = 0.4)
# 
# 
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
# # 
# 
# library(BSgenome.Hsapiens.UCSC.hg38)
# 
# 
# lamlCommun <- maftools::read.maf(maf = Loas_files_SNPs[Loas_files_SNPs$New==6&Loas_files_SNPs$Cancer=="ADC",], 
#                                  vc_nonSyn = unique(Loas_files_SNPs$Variant_Classification))
# laml.tnm = trinucleotideMatrix(maf = lamlCommun, prefix = '', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
# plotApobecDiff(tnm = laml.tnm, maf = lamlCommun, pVal = 0.9)
# 
# 
# 
# 
# 
# 
# 
# 
unique(Loas_files_SNPs$Tumor_Sample_Barcode)
# 
# 
#Loas_files_SNPs$Cancer<-gsub("_crossed","",Loas_files_SNPs$Cancer)
unique(Loas_files_SNPs$Paire)[10]

unique(Loas_files_SNPs$Cancer[Loas_files_SNPs$Tumor_Sample_Barcode=="Paire_15_ADC_GATK_somatic_filtered"])
unique(Loas_files_SNPs$Cancer[Loas_files_SNPs$Tumor_Sample_Barcode=="Paire_15_SCC_GATK_somatic_filtered"])

for (i in 1:length(unique(Loas_files_SNPs$Paire))){
  
  
  
  my_sample1<-Loas_files_SNPs[which(Loas_files_SNPs$Paire==unique(Loas_files_SNPs$Paire)[i]&Loas_files_SNPs$Cancer=="SCC"),]
  my_sample3<-Loas_files_SNPs[which(Loas_files_SNPs$Paire==unique(Loas_files_SNPs$Paire)[i]&Loas_files_SNPs$Cancer=="ADC"),]
  
  unique(my_sample3$Tumor_Sample_Barcode[my_sample3$Cancer=="ADC"])
  unique(my_sample3$Cancer[my_sample3$Tumor_Sample_Barcode=="Paire_15_SCC_GATK_somatic_filtered"])
  
  # #################
  # 
  Maf_communs<-my_sample3[my_sample3$combi%in%my_sample1$combi,]
  Maf_communs$Tumor_Sample_Barcode
  Maf_communs$Variant_Classification[which(Maf_communs$Variant_Classification=="Non exonic")]<-""
  if(nrow(Maf_communs>0)){
    lamlCommun <- maftools::read.maf(maf = Maf_communs, vc_nonSyn = unique(Maf_communs$Variant_Classification))
    maf_df_Commun <- lamlCommun@data[, c("Tumor_Sample_Barcode", "Chromosome", "Start_Position", 
                                         "Reference_Allele", "Tumor_Seq_Allele2")]
    # maf_df_Commun
    sampCommun<-unique(lamlCommun@data$Tumor_Sample_Barcode)
    sampCommun
    maf_df_Commun$Tumor_Sample_Barcode
    
    spemaf_df_Commun<-maf_df_Commun[maf_df_Commun$Tumor_Sample_Barcode==sampCommun,]
    # 
    vcf_df_Commun <- data.frame(
      CHROM = spemaf_df_Commun$Chromosome,
      POS = spemaf_df_Commun$Start_Position,
      ID = ".",
      REF = spemaf_df_Commun$Reference_Allele,
      ALT = spemaf_df_Commun$Tumor_Seq_Allele2,
      QUAL = ".",
      FILTER = "PASS",
      INFO = paste0("SAMPLE=", spemaf_df_Commun$Tumor_Sample_Barcode)
    )
    # 
    # Write VCF with a header
    writeLines(c(
      "##fileformat=VCFv4.2",
      paste0("##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Tumor sample ID\">"),
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    ), con = "output.vcf")
    
    write.table(vcf_df_Commun, file = paste0("./MAFs_COMMUN/",sampCommun,".vcf"), append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  ###
  
  Maf_no_communsADC<-my_sample3[!my_sample3$combi%in%my_sample1$combi,]
  Maf_no_communsSCC<-my_sample1[!my_sample1$combi%in%my_sample3$combi,]
  Maf_no_communsTOT<-unique(rbind(Maf_no_communsADC,Maf_no_communsSCC))
  
  Maf_no_communsADC$Variant_Classification[which(Maf_no_communsADC$Variant_Classification=="Non exonic")]<-""
  Maf_no_communsSCC$Variant_Classification[which(Maf_no_communsSCC$Variant_Classification=="Non exonic")]<-""
  
  ###
  if(nrow(Maf_no_communsADC>0)){
    lamlADC = maftools::read.maf(maf = Maf_no_communsADC,vc_nonSyn = unique(Maf_no_communsADC$Variant_Classification))
    maf_df_ADC <- lamlADC@data[, c("Tumor_Sample_Barcode", "Chromosome", "Start_Position", 
                                   "Reference_Allele", "Tumor_Seq_Allele2")]
    
    sampADC<-unique(lamlADC@data$Tumor_Sample_Barcode)
    spemaf_df_ADC<-maf_df_ADC[maf_df_ADC$Tumor_Sample_Barcode==sampADC,]
    
    vcf_df_ADC <- data.frame(
      CHROM = spemaf_df_ADC$Chromosome,
      POS = spemaf_df_ADC$Start_Position,
      ID = ".",
      REF = spemaf_df_ADC$Reference_Allele,
      ALT = spemaf_df_ADC$Tumor_Seq_Allele2,
      QUAL = ".",
      FILTER = "PASS",
      INFO = paste0("SAMPLE=", spemaf_df_ADC$Tumor_Sample_Barcode)
    )
    
    # Write VCF with a header
    writeLines(c(
      "##fileformat=VCFv4.2",
      paste0("##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Tumor sample ID\">"),
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    ), con = "output.vcf")
    
    write.table(vcf_df_ADC, file = paste0("./MAFs_ADC/",sampADC,".vcf"), append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  ###
  if(nrow(Maf_no_communsSCC>0)){
    lamlSCC = maftools::read.maf(maf = Maf_no_communsSCC,vc_nonSyn = unique(Maf_no_communsSCC$Variant_Classification))
    maf_df_SCC <- lamlSCC@data[, c("Tumor_Sample_Barcode", "Chromosome", "Start_Position", 
                                   "Reference_Allele", "Tumor_Seq_Allele2")]
    
    sampSCC<-unique(lamlSCC@data$Tumor_Sample_Barcode)
    spemaf_df_SCC<-maf_df_SCC[maf_df_SCC$Tumor_Sample_Barcode==sampSCC,]
    # Create a simple VCF-like table
    vcf_df_SCC <- data.frame(
      CHROM = spemaf_df_SCC$Chromosome,
      POS = spemaf_df_SCC$Start_Position,
      ID = ".",
      REF = spemaf_df_SCC$Reference_Allele,
      ALT = spemaf_df_SCC$Tumor_Seq_Allele2,
      QUAL = ".",
      FILTER = "PASS",
      INFO = paste0("SAMPLE=", spemaf_df_SCC$Tumor_Sample_Barcode)
    )
    
    # Write VCF with a header
    writeLines(c(
      "##fileformat=VCFv4.2",
      paste0("##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Tumor sample ID\">"),
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    ), con = "output.vcf")
    
    write.table(vcf_df_SCC, file = paste0("./MAFs_SCC/",sampSCC,".vcf"), append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  
  
  
  ################
  
  
  for (y in 1:length(unique(Loas_files_SNPs$Paire))){
    
    #y<-1
    # i<-12
    # y<-12
    my_sample2<-Loas_files_SNPs[which(Loas_files_SNPs$Paire==unique(Loas_files_SNPs$Paire)[y]&Loas_files_SNPs$Cancer=="ADC"),]
    
    my_sample<-rbind(my_sample1,my_sample2)
    
    SUM_my_sample<-my_sample%>%
      mutate(total_mut=n())%>%
      group_by(Chromosome,Start_Position) %>%
      mutate(iter= n_distinct(Cancer))%>%
      ungroup()%>%
      mutate(n_com=sum(iter == 2),prop_commun=((sum(iter == 2)/2)/(total_mut-(sum(iter == 2)/2)))*100)%>%
      group_by(Paire,Cancer, Variant_Classification,prop_commun,total_mut,n_com) %>%
      summarise(count_variant = n())
    
    
    SUM_my_sample_n<- merge(x=datas_cliniques, y= SUM_my_sample, by = "Paire" )
    
    if(nrow(SUM_my_sample_n)>0){
      
      SUM_my_sample_n$statut<-"CTR"
      SUM_my_sample_n$compare<-paste(i,y,sep = "_")
      
      if(i==y){
        
        SUM_my_sample_n$statut<-as.character(i)
        
      }  
      
      Result_SNPs<-rbind(Result_SNPs,SUM_my_sample_n)
      
    }
  }
}





Result2<-Result_SNPs[-which(Result_SNPs$statut=="CTR"),]
Result3<-Result_SNPs[-which(Result_SNPs$statut=="CTR"),]
Result2$Cancer<-NULL
Result2$Variant_Classification<-NULL
Result2$count_variant<-NULL
Result2<-unique(Result2)

ResultCTR<-Result_SNPs[which(Result_SNPs$statut=="CTR"),]

Result2$ctr<-max(ResultCTR$prop_commun)
max(ResultCTR$prop_commun)
sd(ResultCTR$prop_commun)


colnames(SUM_my_sample)

Result2a<-Result2[which(Result2$FIGO=="1a1"),]
Result2b<-Result2[which(Result2$FIGO=="1b1"),]
Result2c<-Result2[which(Result2$FIGO=="1b2"),]



max_prop<-max(Result_SNPs$prop_commun)

Result2a$Figo_score

Result2a$ctr


df <- data.frame(x = factor(), y = numeric())

# Start the plot
x10<-ggplot(df, aes(x, y)) +
  geom_bar(stat = "identity") +  # Won't actually draw anything because df is empty
  theme_void() +
  labs(y = "Proportion of common INDELs") +
  ylim(0,max_prop)+# Removes almost everything
  ggtitle("FIGO") +
  theme(
    axis.title.y = element_text(size = 18, face = "bold", color = "Grey",angle = 90,hjust = 0.25,vjust = 1.5),
    axis.text.y = element_text(size = 16, face = "bold", color = "black",  hjust = 1, vjust = 0.5  ),
    # axis.ticks.y = element_line(),   # Keep y-axis ticks
    # axis.line.y = element_line(),    # Keep y-axis line
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.45, size = 16)
  )





x1a <-ggplot(data = Result2a, aes(x =  as.character(New), 
                                  y = prop_commun)) +
  geom_bar(stat = "identity", fill = "gray") +
  ylim(0,max_prop)+
  geom_hline(aes(yintercept = ctr), linetype = "dashed") +
  geom_hline(aes(yintercept = 0)) +
  #geom_hline(yintercept = max(ResultCTR$prop_commun), linetype = "dashed")+
  facet_wrap(~ reorder(as.character(New),as.numeric(New)), nrow = 1,scale="free_x") +
  # scale_y_continuous(expand = c(0, 0)) +
  ggtitle("1a1") +
  
  theme(
    axis.text.x = element_text(size = 20, face = "bold", color = "gray", vjust = 1),
    axis.title.x = element_blank(),
    axis.text.y =  element_blank(),
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = unit(c(0, 1, 1, 1), "cm"),
    strip.text  = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.45, size = 16)
  )


x1b <- ggplot(data = Result2b, aes(x =  as.character(Result2b$New), 
                                   y = prop_commun)) +
  geom_bar(stat = "identity", fill = "gray") +
  ylim(0,max_prop)+
  geom_hline(aes(yintercept = 0)) +
  #geom_hline(data = ResultCTR, aes(yintercept = max(ResultCTR$prop_commun)), linetype = "dashed") +
  geom_hline(aes(yintercept = ctr), linetype = "dashed") +
  facet_wrap(~ reorder(as.character(New),as.numeric(New)), nrow = 1,scale="free_x") +
  # scale_y_continuous(expand = c(0, 0)) +
  ggtitle("1b1") +
  theme(
    axis.text.x = element_text(size = 20, face = "bold", color = "gray", vjust = 1),
    axis.title.x = element_blank(),
    axis.text.y =  element_blank(),
    axis.ticks.y = element_blank(),
    panel.background = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    plot.margin = unit(c(0, 1, 1, 1), "cm"),
    strip.text = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.45, size = 16)
  )

x1c <- ggplot(data = Result2c, aes(x = as.character(Result2c$New), 
                                   y = prop_commun)) +
  geom_bar(stat = "identity", fill = "gray") +
  ylim(0,max_prop)+
  geom_hline(aes(yintercept = ctr), linetype = "dashed") +
  geom_hline(aes(yintercept = 0)) +
  facet_wrap(~ reorder(as.character(New),as.numeric(New)), nrow = 1,scale="free_x") +
  # scale_y_continuous(expand = c(0, 0)) +
  ggtitle("1b2") +
  theme(
    axis.text.x = element_text(size = 20, face = "bold", color = "gray", vjust = 1),
    axis.title.x = element_blank(),
    axis.text.y =  element_blank(),
    panel.background = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    plot.margin = unit(c(0, 1, 1, 1), "cm"),
    strip.text = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.45, size = 16)
  )


####
Result3 <- Result_SNPs[Result_SNPs$statut != "CTR",]
Result3a<-Result3[which(Result3$FIGO=="1a1"),]
Result3b<-Result3[which(Result3$FIGO=="1b1"),]
Result3c<-Result3[which(Result3$FIGO=="1b2"),]


# Step 2: Get all unique Variant_Classification values
all_variants <- unique(Result3$Variant_Classification)

# Step 3: Create a named color palette
colors <- RColorBrewer::brewer.pal(n = length(all_variants), name = "Set2")  # or any palette you like
names(colors) <- all_variants

#####"


Result3<-Result3%>%
  group_by(Paire,Cancer)%>%
  mutate(tot_sum=sum(count_variant))

max_count<-max(Result3$tot_sum)
######"


max_count<-max(Result3$tot_sum)
# Start the plot
x20<-ggplot(df, aes(x, y)) +
  geom_bar(stat = "identity") +  # Won't actually draw anything because df is empty
  theme_void() +
  labs(y = "Number of INDELs") +
  ylim(as.numeric(max_count),0)+
  theme(
    axis.title.y = element_text(size = 18, face = "bold", color = "Grey",angle = 90,hjust = 0.25,vjust = 1.5),
    axis.text.y = element_text(size = 16, face = "bold", color = "black",  hjust = 1, vjust = 0.5  ),
    # axis.ticks.y = element_line(),   # Keep y-axis ticks
    # axis.line.y = element_line(),    # Keep y-axis line
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  )+
  guides(x = ggh4x::guide_axis_manual(position = "top"))

x2a <- ggplot(data = Result3a,
              aes(x = Cancer, 
                  y = count_variant, 
                  fill = Variant_Classification)) +
  geom_bar(stat = "identity", color = "black") +
  geom_hline(aes(yintercept = 0)) +
  ylim(as.numeric(max_count),0)+
  scale_fill_manual(values = variant_colors, drop = FALSE) +
  # scale_fill_brewer(palette = "Set2") + 
  facet_wrap(~ reorder(as.character(New),as.numeric(New)), nrow = 1,scale="free_x") +
  #scale_y_continuous(trans = "reverse") +
  theme(
    axis.text.x = element_text(size = 16, face = "bold", color = "gray", vjust = 0.5,angle = 45,hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.y =  element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    plot.margin = unit(c(0, 1, 1, 1), "cm"),
    strip.text = element_blank()
  )+
  guides(x = ggh4x::guide_axis_manual(position = "top"))


x2b <- ggplot(data = Result3b,
              aes(x = Cancer, 
                  y = count_variant, 
                  fill = Variant_Classification)) +
  geom_bar(stat = "identity", color = "black") +
  geom_hline(aes(yintercept = 0)) +
  #scale_fill_brewer(palette = "Set2") + 
  scale_fill_manual(values = variant_colors, drop = FALSE) +
  ylim(as.numeric(max_count),0)+
  facet_wrap(~ reorder(as.character(New),as.numeric(New)), nrow = 1,scale="free_x") +
  # scale_y_continuous(trans = "reverse") +
  theme(
    axis.text.x = element_text(size = 16, face = "bold", color = "gray", vjust = 0.5,angle = 45,hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.y =  element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    plot.margin = unit(c(0, 1, 1, 1), "cm"),
    strip.text = element_blank()
  )+
  guides(x = ggh4x::guide_axis_manual(position = "top"))


x2c <- ggplot(data = Result3c,
              aes(x = Cancer, 
                  y = count_variant, 
                  fill = Variant_Classification)) +
  geom_bar(stat = "identity", color = "black") +
  geom_hline(aes(yintercept = 0)) +
  #  scale_fill_brewer(palette = "Set2") + 
  scale_fill_manual(values = variant_colors, drop = FALSE) +
  ylim(as.numeric(max_count),0)+
  facet_wrap(~ reorder(as.character(New),as.numeric(New)), nrow = 1,scale="free_x") +
  #scale_y_continuous(trans = "reverse") +
  theme(
    axis.text.x =  element_text(size = 16, face = "bold", color = "gray", vjust = 0.5,angle = 45,hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.y =  element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    plot.margin = unit(c(0, 1, 1, 1), "cm"),
    strip.text = element_blank()
  )+
  guides(x = ggh4x::guide_axis_manual(position = "top"))


widths_top <- c(
  0.5,  # For x10 (y-axis-only)
  length(unique(Result2a$New)),
  length(unique(Result2b$New)),
  length(unique(Result2c$New))
)

widths_bottom <- c(
  0.5,  # For x20 (y-axis-only)
  length(unique(Result3a$New)),
  length(unique(Result3b$New)),
  length(unique(Result3c$New))
)



top_row <- x10 + x1a + x1b + x1c + 
  plot_layout(ncol = 4, guides = "collect", widths = widths_top)+ theme(plot.margin = margin(0, 1, 0, 1, "cm"))

bottom_row <- x20 + x2a + x2b + x2c + 
  plot_layout(ncol = 4, guides = "collect", widths = widths_bottom)+ theme(plot.margin = margin(0, 1, 0, 1, "cm"))


final_plot <- top_row / bottom_row + plot_layout(heights = c(1, 1),tag_level = "keep",guides = "collect")&theme(plot.margin = margin(0, 0, 0, 0)) 

final_plot

ggsave("Part1_filtr15_commun_newnewnewFig2_SNV_commun.pdf", plot = final_plot, width = 12, height = 8, limitsize = FALSE, device = 'pdf', dpi = 300)
#ggsave("filtr8_newnewnewFig2_SNV_commun.png", plot = final_plot, width = 12, height = 8, limitsize = FALSE, device = 'png', dpi = 300)

dev.off()


###########


Result4<-Result3%>%
  group_by(Paire,FIGO,Cancer,Integration,type)%>%
  summarise(counting=sum(count_variant))

Result4$Integration[Result4$Integration=="INTE_VS"]<-"INTE"
Result4$Integration[Result4$Integration=="VS"]<-"EPI"


TMB_plot1 <- ggplot(Result4, aes(x = type, y = counting), fill = "grey") +
  geom_violin(alpha = 0.6, fill = "gray") +  # Violin plot for distributione
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(shape = 16, size = 2, width = 0.15, alpha = 0.7) +  # Add points for individual values
  scale_y_log10() +  # Use log scale for better visualization
  scale_fill_brewer(palette = "Set2") +  # Nice color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text( size = 14,angle = 45,hjust = 1,vjust = 1),
    axis.text.y =  element_text( size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "top"
  ) +
  
  
  labs(y = "Number of mutation", title = "HPV types")+
  stat_compare_means(method = "anova", label.y = 1,label="p.format",label.x.npc = "right") 


colnames(Result4)
TMB_plot2 <- ggplot(Result4, aes(x = FIGO, y = counting), fill = "grey") +
  geom_violin(alpha = 0.6, fill = "gray") +  # Violin plot for distributione
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(shape = 16, size = 2, width = 0.15, alpha = 0.7) +  # Add points for individual values
  scale_y_log10() +  # Use log scale for better visualization
  scale_fill_brewer(palette = "Set2") +  # Nice color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text( size = 14,angle = 45,hjust = 1,vjust = 1),
    axis.text.y =  element_text( size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "top"
  ) +
  labs(y = "Number of mutation", title = "FIGO cancer classification")+
  stat_compare_means(method = "anova", label.y = 1,label="p.format",label.x.npc = "right") 




TMB_plot3 <- ggplot(Result4, aes(x = Integration, y = counting)) +
  geom_violin(alpha = 0.6, fill = "gray") +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(shape = 16, size = 2, width = 0.15, alpha = 0.7) +
  scale_y_log10() +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "top"
  ) +
  labs(y = "Number of mutation", title = "Integration Statut") +
  
  # 🟡 Add statistical test brackets
  stat_compare_means(method = "anova", label.y = 1,label="p.format",label.x.npc = "right") 
# stat_compare_means(
#   comparisons = list(
#     c("group1", "group2"),  # Replace with your actual group names
#     c("group2", "group3")
#     # Add more if needed
#   ),
#   method = "t.test",
#   label = "p.signif",  # You can use "p.format" for actual p-value
#   bracket.size = 0.5,
#   tip.length = 0.03,
#   size = 5
# )


TMB_plot4 <- ggplot(Result4, aes(x = Cancer, y = counting)) +
  geom_violin(alpha = 0.6, fill = "gray") +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(shape = 16, size = 2, width = 0.15, alpha = 0.7) +
  scale_y_log10() +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "top"
  ) +
  labs(y = "Number of mutation", title = "Integration Statut") +
  
  # 🟡 Add statistical test brackets
  stat_compare_means(method = "anova", label.y = 1,label="p.format",label.x.npc = "right") 

combi_row <- TMB_plot1 + TMB_plot2 + TMB_plot3 + TMB_plot4+
  plot_layout(ncol = 4, guides = "collect")


combi_row

# Save the plot
ggsave("filtr15_newnewTMB_plot_by_Strain2.pdf", width = 12, height = 6, dpi = 300)


Result3snps<-Result3

#########################################################################################################


####

################################################################################
##########################################  InDels ######################################################

min(Loas_files_INDELs_mutect2$TLOD)

unique(Loas_files_INDELs_mutect2$Tumor_Sample_Barcode[Loas_files_INDELs_mutect2$Paire==15&Loas_files_INDELs_mutect2$Cancer=="ADC"])

unique(Loas_files_INDELs_mutect2$Otherinfo10[Loas_files_INDELs_mutect2$Paire==15&Loas_files_INDELs_mutect2$Cancer=="ADC"])

min(Loas_files_INDELs_mutect2$Otherinfo10[Loas_files_INDELs_mutect2$Paire==15&Loas_files_INDELs_mutect2$Cancer=="ADC"])

xx<-table(Loas_files_INDELs_mutect2$Otherinfo8[Loas_files_INDELs_mutect2$Paire == 15 & Loas_files_INDELs_mutect2$Cancer == "ADC"])
table(Loas_files_INDELs_mutect2$Otherinfo10)

df_counts <- as.data.frame(table(
  Loas_files_INDELs_mutect2$Otherinfo8[
    Loas_files_INDELs_mutect2$Paire == 15 &
      Loas_files_INDELs_mutect2$Cancer == "ADC"
  ]
))


df_counts_snps <- as.data.frame(table(
  Loas_files_SNPs$Otherinfo10[
    Loas_files_SNPs$Paire == 15 &
      Loas_files_SNPs$Cancer == "ADC"
  ]
))



grep("(A{1,}|T{1,}|C{1,}|G{1,})", df$Otherinfo10, value = TRUE)
grepl("(A{3,}|T{3,}|C{3,}|G{3,})", df_counts$Otherinfo10, perl = TRUE)

Result_SNPs<- data.frame(Paire = integer(),
                         ADC = integer(),
                         SCC = integer(),
                         class = numeric(),
                         New = numeric(),
                         Subtype = integer(),
                         Integration = integer(),
                         Cancer = integer(),
                         Variant_Classification = integer(),
                         prop_commun = numeric(),
                         count_variant = numeric(),
                         statut = integer(),
                         compare = integer(),
                         stringsAsFactors = FALSE)



Loas_files_INDELs_mutect2$Variant_Classification[which(Loas_files_INDELs_mutect2$Variant_Classification==""| Loas_files_INDELs_mutect2$Variant_Classification=="Unknown")]<-"Non exonic"

#Loas_files_INDELs_mutect2<-Loas_files_INDELs_mutect2[!Loas_files_INDELs_mutect2$Variant_Classification=="Non exonic",]
unique(Loas_files_INDELs_mutect2$Paire)


for (i in 1:length(unique(Loas_files_INDELs_mutect2$Paire))){
  my_sample1<-Loas_files_INDELs_mutect2[which(Loas_files_INDELs_mutect2$Paire==unique(Loas_files_INDELs_mutect2$Paire)[i]&Loas_files_INDELs_mutect2$Cancer=="SCC"),]
  my_sample3<-Loas_files_INDELs_mutect2[which(Loas_files_INDELs_mutect2$Paire==unique(Loas_files_INDELs_mutect2$Paire)[i]&Loas_files_INDELs_mutect2$Cancer=="ADC"),]
  # #################
  # 
  Maf_communs<-my_sample3[my_sample3$combi%in%my_sample1$combi,]
  Maf_communs$Variant_Classification[which(Maf_communs$Variant_Classification=="Non exonic")]<-""
  if(nrow(Maf_communs>0)){
    lamlCommun <- maftools::read.maf(maf = Maf_communs, vc_nonSyn = unique(Maf_communs$Variant_Classification))
    maf_df_Commun <- lamlCommun@data[, c("Tumor_Sample_Barcode", "Chromosome", "Start_Position", 
                                         "Reference_Allele", "Tumor_Seq_Allele2")]
    # maf_df_Commun
    sampCommun<-unique(lamlCommun@data$Tumor_Sample_Barcode)
    spemaf_df_Commun<-maf_df_Commun[maf_df_Commun$Tumor_Sample_Barcode==sampCommun,]
    # 
    vcf_df_Commun <- data.frame(
      CHROM = spemaf_df_Commun$Chromosome,
      POS = spemaf_df_Commun$Start_Position,
      ID = ".",
      REF = spemaf_df_Commun$Reference_Allele,
      ALT = spemaf_df_Commun$Tumor_Seq_Allele2,
      QUAL = ".",
      FILTER = "PASS",
      INFO = paste0("SAMPLE=", spemaf_df_Commun$Tumor_Sample_Barcode)
    )
    # 
    # Write VCF with a header
    writeLines(c(
      "##fileformat=VCFv4.2",
      paste0("##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Tumor sample ID\">"),
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    ), con = "output.vcf")
    
    write.table(vcf_df_Commun, file = paste0("./MAFs_COMMUN/",sampCommun,"INDELs.vcf"), append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  ###
  
  Maf_no_communsADC<-my_sample3[!my_sample3$combi%in%my_sample1$combi,]
  Maf_no_communsSCC<-my_sample1[!my_sample1$combi%in%my_sample3$combi,]
  Maf_no_communsTOT<-unique(rbind(Maf_no_communsADC,Maf_no_communsSCC))
  
  Maf_no_communsADC$Variant_Classification[which(Maf_no_communsADC$Variant_Classification=="Non exonic")]<-""
  Maf_no_communsSCC$Variant_Classification[which(Maf_no_communsSCC$Variant_Classification=="Non exonic")]<-""
  
  ###
  if(nrow(Maf_no_communsADC>0)){
    lamlADC = maftools::read.maf(maf = Maf_no_communsADC,vc_nonSyn = unique(Maf_no_communsADC$Variant_Classification))
    maf_df_ADC <- lamlADC@data[, c("Tumor_Sample_Barcode", "Chromosome", "Start_Position", 
                                   "Reference_Allele", "Tumor_Seq_Allele2")]
    
    sampADC<-unique(lamlADC@data$Tumor_Sample_Barcode)
    spemaf_df_ADC<-maf_df_ADC[maf_df_ADC$Tumor_Sample_Barcode==sampADC,]
    
    vcf_df_ADC <- data.frame(
      CHROM = spemaf_df_ADC$Chromosome,
      POS = spemaf_df_ADC$Start_Position,
      ID = ".",
      REF = spemaf_df_ADC$Reference_Allele,
      ALT = spemaf_df_ADC$Tumor_Seq_Allele2,
      QUAL = ".",
      FILTER = "PASS",
      INFO = paste0("SAMPLE=", spemaf_df_ADC$Tumor_Sample_Barcode)
    )
    
    # Write VCF with a header
    writeLines(c(
      "##fileformat=VCFv4.2",
      paste0("##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Tumor sample ID\">"),
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    ), con = "output.vcf")
    
    write.table(vcf_df_ADC, file = paste0("./MAFs_ADC/",sampADC,"INDELs.vcf"), append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  ###
  if(nrow(Maf_no_communsSCC>0)){
    lamlSCC = maftools::read.maf(maf = Maf_no_communsSCC,vc_nonSyn = unique(Maf_no_communsSCC$Variant_Classification))
    maf_df_SCC <- lamlSCC@data[, c("Tumor_Sample_Barcode", "Chromosome", "Start_Position", 
                                   "Reference_Allele", "Tumor_Seq_Allele2")]
    
    sampSCC<-unique(lamlSCC@data$Tumor_Sample_Barcode)
    spemaf_df_SCC<-maf_df_SCC[maf_df_SCC$Tumor_Sample_Barcode==sampSCC,]
    # Create a simple VCF-like table
    vcf_df_SCC <- data.frame(
      CHROM = spemaf_df_SCC$Chromosome,
      POS = spemaf_df_SCC$Start_Position,
      ID = ".",
      REF = spemaf_df_SCC$Reference_Allele,
      ALT = spemaf_df_SCC$Tumor_Seq_Allele2,
      QUAL = ".",
      FILTER = "PASS",
      INFO = paste0("SAMPLE=", spemaf_df_SCC$Tumor_Sample_Barcode)
    )
    
    # Write VCF with a header
    writeLines(c(
      "##fileformat=VCFv4.2",
      paste0("##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Tumor sample ID\">"),
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    ), con = "output.vcf")
    
    write.table(vcf_df_SCC, file = paste0("./MAFs_SCC/",sampADC,"INDELs.vcf"), append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  
  
  ################
  for (y in 1:length(unique(Loas_files_INDELs_mutect2$Paire))){
    Loas_files_INDELs_mutect2$Reference_Allele
    my_sample2<-Loas_files_INDELs_mutect2[which(Loas_files_INDELs_mutect2$Paire==unique(Loas_files_INDELs_mutect2$Paire)[y]&Loas_files_INDELs_mutect2$Cancer=="ADC"),]
    my_sample<-rbind(my_sample1,my_sample2)
    # 
    # my_sample1$combi3
    # my_sample2$combi3
    # my_sample_commun<-my_sample1[my_sample1$combi3%in%my_sample2$combi3,]
    # my_sample_commun2<-my_sample2[my_sample2$combi3%in%my_sample1$combi3,]

    
    # Loas_files_INDELs_mutect2$Variant_Classification
    SUM_my_sample<-my_sample%>%
      mutate(total_mut=n())%>%
      group_by(Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2,End_Position) %>%
      mutate(iter= n_distinct(Cancer))%>%
      ungroup()%>%
      mutate(n_com=sum(iter == 2),prop_commun=((sum(iter == 2)/2)/(total_mut-(sum(iter == 2)/2)))*100)%>%
      group_by(Paire,Cancer, Variant_Classification,Variant_Type,prop_commun,total_mut,n_com) %>%
      summarise(count_variant = n())
    
    SUM_my_sample_n<- merge(x=datas_cliniques, y= SUM_my_sample, by = "Paire" )
    
    if(nrow(SUM_my_sample_n)>0){
      
      SUM_my_sample_n$statut<-"CTR"
      SUM_my_sample_n$compare<-paste(i,y,sep = "_")
      
      if(i==y){
        
        SUM_my_sample_n$statut<-as.character(i)
        
      }  
      
      Result_SNPs<-rbind(Result_SNPs,SUM_my_sample_n)
      
    }
  }
}



Result2<-Result_SNPs[-which(Result_SNPs$statut=="CTR"),]
Result3<-Result_SNPs[-which(Result_SNPs$statut=="CTR"),]
Result2$Cancer<-NULL
Result2$Variant_Classification<-NULL
Result2$Variant_Type<-NULL
Result2$count_variant<-NULL
Result2<-unique(Result2)

ResultCTR<-Result_SNPs[which(Result_SNPs$statut=="CTR"),]
ResultCTR
Result2$ctr<-max(ResultCTR$prop_commun)

Result2a<-Result2[which(Result2$FIGO=="1a1"),]
Result2b<-Result2[which(Result2$FIGO=="1b1"),]
Result2c<-Result2[which(Result2$FIGO=="1b2"),]



max_prop<-max(Result_SNPs$prop_commun)
max_count<-max(Result_SNPs$count_variant)
Result2a$Figo_score



#####"


Result3<-Result3%>%
  group_by(Paire,Cancer)%>%
  mutate(tot_sum=sum(count_variant))

max_count<-max(Result3$tot_sum)
######"


df <- data.frame(x = factor(), y = numeric())

# Start the plot
x10<-ggplot(df, aes(x, y)) +
  geom_bar(stat = "identity") +  # Won't actually draw anything because df is empty
  theme_void() +
  labs(y = "Proportion of common INDELs") +
  ylim(0,max_prop)+# Removes almost everything
  ggtitle("FIGO") +
  theme(
    axis.title.y = element_text(size = 18, face = "bold", color = "Grey",angle = 90,hjust = 0.25,vjust = 1.5),
    axis.text.y = element_text(size = 16, face = "bold", color = "black",  hjust = 1, vjust = 0.5  ),
    # axis.ticks.y = element_line(),   # Keep y-axis ticks
    # axis.line.y = element_line(),    # Keep y-axis line
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.45, size = 16)
  )





x1a <-ggplot(data = Result2a, aes(x =  as.character(New), 
                                  y = prop_commun)) +
  geom_bar(stat = "identity", fill = "gray") +
  ylim(0,max_prop)+
  geom_hline(aes(yintercept = ctr), linetype = "dashed") +
  geom_hline(aes(yintercept = 0)) +
  #geom_hline(yintercept = max(ResultCTR$prop_commun), linetype = "dashed")+
  facet_wrap(~ reorder(as.character(New),as.numeric(New)), nrow = 1,scale="free_x") +
  # scale_y_continuous(expand = c(0, 0)) +
  ggtitle("1a1") +
  
  theme(
    axis.text.x = element_text(size = 20, face = "bold", color = "gray", vjust = 1),
    axis.title.x = element_blank(),
    axis.text.y =  element_blank(),
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = unit(c(0, 1, 1, 1), "cm"),
    strip.text  = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.45, size = 16)
  )


x1b <- ggplot(data = Result2b, aes(x =  as.character(Result2b$New), 
                                   y = prop_commun)) +
  geom_bar(stat = "identity", fill = "gray") +
  ylim(0,max_prop)+
  geom_hline(aes(yintercept = 0)) +
  #geom_hline(data = ResultCTR, aes(yintercept = max(ResultCTR$prop_commun)), linetype = "dashed") +
  geom_hline(aes(yintercept = ctr), linetype = "dashed") +
  facet_wrap(~ reorder(as.character(New),as.numeric(New)), nrow = 1,scale="free_x") +
  # scale_y_continuous(expand = c(0, 0)) +
  ggtitle("1b1") +
  theme(
    axis.text.x = element_text(size = 20, face = "bold", color = "gray", vjust = 1),
    axis.title.x = element_blank(),
    axis.text.y =  element_blank(),
    axis.ticks.y = element_blank(),
    panel.background = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    plot.margin = unit(c(0, 1, 1, 1), "cm"),
    strip.text = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.45, size = 16)
  )
x1b
x1c <- ggplot(data = Result2c, aes(x = as.character(Result2c$New), 
                                   y = prop_commun)) +
  geom_bar(stat = "identity", fill = "gray") +
  ylim(0,max_prop)+
  geom_hline(aes(yintercept = ctr), linetype = "dashed") +
  geom_hline(aes(yintercept = 0)) +
  facet_wrap(~ reorder(as.character(New),as.numeric(New)), nrow = 1,scale="free_x") +
  # scale_y_continuous(expand = c(0, 0)) +
  ggtitle("1b2") +
  theme(
    axis.text.x = element_text(size = 20, face = "bold", color = "gray", vjust = 1),
    axis.title.x = element_blank(),
    axis.text.y =  element_blank(),
    panel.background = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    plot.margin = unit(c(0, 1, 1, 1), "cm"),
    strip.text = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.45, size = 16)
  )


####
Result3 <- Result_SNPs[Result_SNPs$statut != "CTR",]

unique(Result3$Variant_Type)

Result3_complet<-Result3%>%
  dplyr::select(New,Variant_Classification,Cancer,count_variant)%>%
  complete(
    New=c(1:12),
    Variant_Classification=c("Frame_Shift_Del","Inframe_INDEL","Non exonic","Nonsense_Mutation","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Translation_Start_Site","Nonstop_Mutation"),
    Cancer = c("ADC","SCC"),
    fill = list(count_variant = 0)
  )

Result3_complet$Variant_Classification[Result3_complet$Variant_Classification==("Frame_Shift_Del")|Result3_complet$Variant_Classification==("Frame_Shift_Ins")]<-"Frame_Shift"

Result3_complet$Variant_Classification[Result3_complet$Variant_Classification==("Inframe_INDEL")|
                                         Result3_complet$Variant_Classification==("In_Frame_Del")|
                                         Result3_complet$Variant_Classification==("In_Frame_Ins")]<-"Inframe"

library(dplyr)
library(tidyr)
Result3_complet2<-Result3_complet%>%
  group_by(New,Variant_Classification,Cancer)%>%
  mutate(count_variant=sum(count_variant))%>%
  ungroup()%>%
  dplyr::select(New,Variant_Classification,Cancer,count_variant)%>%
  distinct()%>%
  ungroup()

Result3_complet2<-merge(datas_cliniques, Result3_complet2, by.x = "New", by.y = "New")

Result3a<-Result3_complet2[which(Result3_complet2$FIGO=="1a1"),]
Result3b<-Result3_complet2[which(Result3_complet2$FIGO=="1b1"),]
Result3c<-Result3_complet2[which(Result3_complet2$FIGO=="1b2"),]

# Step 2: Get all unique Variant_Classification values
all_variants <- unique(Result3_complet2$Variant_Classification)
all_variants
# Step 3: Create a named color palette
colors <- RColorBrewer::brewer.pal(n = length(all_variants), name = "Set3")  # or any palette you like
names(colors) <- all_variants


# Start the plot
x20<-ggplot(df, aes(x, y)) +
  geom_bar(stat = "identity") +  # Won't actually draw anything because df is empty
  theme_void() +
  labs(y = "Number of INDELs") +
  ylim(as.numeric(max_count),0)+
  theme(
    axis.title.y = element_text(size = 18, face = "bold", color = "Grey",angle = 90,hjust = 0.25,vjust = 1.5),
    axis.text.y = element_text(size = 16, face = "bold", color = "black",  hjust = 1, vjust = 0.5  ),
    # axis.ticks.y = element_line(),   # Keep y-axis ticks
    # axis.line.y = element_line(),    # Keep y-axis line
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  )+
  guides(x = ggh4x::guide_axis_manual(position = "top"))

x2a <- ggplot(data = Result3a,
              aes(x = Cancer, 
                  y = count_variant, 
                  fill = Variant_Classification)) +
  geom_bar(stat = "identity", color = "black") +
  geom_hline(aes(yintercept = 0)) +
  ylim(as.numeric(max_count),0)+
  scale_fill_manual(values = variant_colors, drop = FALSE) +
  #scale_fill_brewer(palette = "Set2") + 
  facet_wrap(~ reorder(as.character(New),as.numeric(New)), nrow = 1,scale="free_x") +
  #scale_y_continuous(trans = "reverse") +
  theme(
    axis.text.x = element_text(size = 16, face = "bold", color = "gray", vjust = 0.5,angle = 45,hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.y =  element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    plot.margin = unit(c(0, 1, 1, 1), "cm"),
    strip.text = element_blank()
  )+
  guides(x = ggh4x::guide_axis_manual(position = "top"))


x2b <- ggplot(data = Result3b,
              aes(x = Cancer, 
                  y = count_variant, 
                  fill = Variant_Classification)) +
  geom_bar(stat = "identity", color = "black") +
  geom_hline(aes(yintercept = 0)) +
  # scale_fill_brewer(palette = "Set2") + 
  scale_fill_manual(values = variant_colors, drop = FALSE) +
  ylim(as.numeric(max_count),0)+
  facet_wrap(~ reorder(as.character(New),as.numeric(New)), nrow = 1,scale="free_x") +
  # scale_y_continuous(trans = "reverse") +
  theme(
    axis.text.x = element_text(size = 16, face = "bold", color = "gray", vjust = 0.5,angle = 45,hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.y =  element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    plot.margin = unit(c(0, 1, 1, 1), "cm"),
    strip.text = element_blank()
  )+
  guides(x = ggh4x::guide_axis_manual(position = "top"))


x2c <- ggplot(data = Result3c,
              aes(x = Cancer, 
                  y = count_variant, 
                  fill = Variant_Classification)) +
  geom_bar(stat = "identity", color = "black") +
  geom_hline(aes(yintercept = 0)) +
  # scale_fill_brewer(palette = "Set2") + 
  scale_fill_manual(values = variant_colors, drop = FALSE) +
  ylim(as.numeric(max_count),0)+
  facet_wrap(~ reorder(as.character(New),as.numeric(New)), nrow = 1,scale="free_x") +
  #scale_y_continuous(trans = "reverse") +
  theme(
    axis.text.x =  element_text(size = 16, face = "bold", color = "gray", vjust = 0.5,angle = 45,hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.y =  element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    plot.margin = unit(c(0, 1, 1, 1), "cm"),
    strip.text = element_blank()
  )+
  guides(x = ggh4x::guide_axis_manual(position = "top"))


widths_top <- c(
  0.5,  # For x10 (y-axis-only)
  length(unique(Result2a$New)),
  length(unique(Result2b$New)),
  length(unique(Result2c$New))
)

widths_bottom <- c(
  0.5,  # For x20 (y-axis-only)
  length(unique(Result3a$New)),
  length(unique(Result3b$New)),
  length(unique(Result3c$New))
)


top_row <- x10 + x1a + x1b + x1c + 
  plot_layout(ncol = 4, guides = "collect", widths = widths_top)+ theme(plot.margin = margin(0, 1, 0, 1, "cm"))

bottom_row <- x20 + x2a + x2b + x2c + 
  plot_layout(ncol = 4, guides = "collect", widths = widths_bottom)+ theme(plot.margin = margin(0, 1, 0, 1, "cm"))


x2a

final_plot <- top_row / bottom_row + plot_layout(heights = c(1, 1),tag_level = "keep",guides = "collect")&theme(plot.margin = margin(0, 0, 0, 0)) 
top_row
final_plot



ggsave("15_newnewFig2_INDELs_commun_dplyr.pdf", plot = final_plot, width = 12, height = 8, limitsize = FALSE, device = 'pdf', dpi = 300)
ggsave("8_newnewFig2_INDELs_commun_select.png", plot = final_plot, width = 12, height = 8, limitsize = FALSE, device = 'png', dpi = 300)

colnames(Result3)

###################
#############
Result4<-Result3%>%
  group_by(Paire,FIGO,Cancer,Integration,type)%>%
  summarise(counting=sum(count_variant))


Result4$Integration[Result4$Integration=="INTE_VS"]<-"INTE"
Result4$Integration[Result4$Integration=="VS"]<-"EPI"

TMB_plot1 <- ggplot(Result4, aes(x = type, y = counting), fill = "grey") +
  geom_violin(alpha = 0.6, fill = "gray") +  # Violin plot for distributione
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(shape = 16, size = 2, width = 0.15, alpha = 0.7) +  # Add points for individual values
  scale_y_log10() +  # Use log scale for better visualization
  scale_fill_brewer(palette = "Set2") +  # Nice color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text( size = 14,angle = 45,hjust = 1,vjust = 1),
    axis.text.y =  element_text( size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "top"
  ) +
  
  labs(y = "Number of mutation", title = "HPV subtypes")+
  stat_compare_means(method = "anova", label.y = 1,label="p.format",label.x.npc = "right") 


colnames(Result4)
TMB_plot2 <- ggplot(Result4, aes(x = FIGO, y = counting), fill = "grey") +
  geom_violin(alpha = 0.6, fill = "gray") +  # Violin plot for distributione
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(shape = 16, size = 2, width = 0.15, alpha = 0.7) +  # Add points for individual values
  scale_y_log10() +  # Use log scale for better visualization
  scale_fill_brewer(palette = "Set2") +  # Nice color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text( size = 14,angle = 45,hjust = 1,vjust = 1),
    axis.text.y =  element_text( size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "top"
  ) +
  labs(y = "Number of mutation", title = "FIGO cancer classification")+
  stat_compare_means(method = "anova", label.y = 1,label="p.format",label.x.npc = "right") 



TMB_plot3 <- ggplot(Result4, aes(x = Integration, y = counting)) +
  geom_violin(alpha = 0.6, fill = "gray") +  # Violin plot for distributione
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(shape = 16, size = 2, width = 0.15, alpha = 0.7) +  # Add points for individual values
  scale_y_log10() +  # Use log scale for better visualization
  scale_fill_brewer(palette = "Set2") +  # Nice color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text( size = 14,angle = 45,hjust = 1,vjust = 1),
    axis.text.y =  element_text( size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "top"
  ) +
  labs(y = "Number of mutation", title = "Integration Statut")+
  stat_compare_means(method = "anova", label.y = 1,label="p.format",label.x.npc = "right") 
TMB_plot3




TMB_plot4 <- ggplot(Result4, aes(x = Cancer, y = counting)) +
  geom_violin(alpha = 0.6, fill = "gray") +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(shape = 16, size = 2, width = 0.15, alpha = 0.7) +
  scale_y_log10() +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "top"
  ) +
  labs(y = "Number of mutation", title = "Integration Statut") +
  
  # 🟡 Add statistical test brackets
  stat_compare_means(method = "anova", label.y = 1,label="p.format",label.x.npc = "right") 

combi_row <- TMB_plot1 + TMB_plot2 + TMB_plot3 + TMB_plot4+
  plot_layout(ncol = 4, guides = "collect")

combi_row
Result2c$ctr
# Save the plot
ggsave("15_TMB_plot_by_Strain2_INDELs_select.pdf", width = 12, height = 6, dpi = 300)

Result3indels<-Result3

dev.off()

#########################################################################################################
##########################################  fill tmb ######################################################


common_cols <- intersect(colnames(Result3snps), colnames(Result3indels))

Result3snps55<-Result3snps%>%
  dplyr::select(common_cols)

Result3indels55<-Result3indels%>%
  dplyr::select(common_cols)

Result32<-rbind(Result3snps55,Result3indels55)

Result4<-Result32%>%
  group_by(Paire,FIGO,Cancer,Integration,type)%>%
  summarise(counting=sum(count_variant))

Result41<-Result32%>%
  group_by(Paire,FIGO,Integration,type)%>%
  summarise(counting=sum(count_variant))


Result4$Integration[Result4$Integration=="INTE_VS"]<-"INTE"
Result4$Integration[Result4$Integration=="VS"]<-"EPI"



TMB_plot1 <- ggplot(Result4, aes(x = Cancer, y = counting)) +
  geom_violin(alpha = 1, fill = "white",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  geom_violin(alpha = 0.4, aes(fill = Cancer),col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +  # Violin plot for distributione
  #geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(shape = 16, size = 3, width = 0.15, alpha = 0.7) +
  scale_y_log10() +
  #scale_fill_brewer(palette = "Set2") +
  
  scale_fill_manual(
    values = c(
      "ADC" = "#FF6666",
      "SCC" = "#3399ff",
      "Commun" = "gray"
    )
  )+
  

  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  labs(y = "Number of mutation", title = "Integration Statut") +
  
  # 🟡 Add statistical test brackets
  stat_compare_means(method = "t.test", label.y = 1,label="p.format",label.x.npc = "right") 

TMB_plot1

Result4$type[Result4$type=="HPV58"]<-"Other HPVs"
Result4$type[Result4$type=="HPV45"]<-"Other HPVs"


TMB_plot2 <- ggplot(Result4, aes(x = type, y = counting)) +
  geom_violin(alpha = 1, fill = "white",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  geom_violin(alpha = 0.4, fill = "gray",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  #geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(aes(col=Cancer),shape = 16, size = 3, width = 0.15, alpha = 0.7) +  # Add points for individual values
  #scale_y_log10() +  # Use log scale for better visualization
  scale_color_manual(
    values = c(
      "ADC" = "#FF6666",
      "SCC" = "#3399ff",
      "Commun" = "grey"
    )
  )+
  scale_y_log10() +  # Use log scale for better visualization
  scale_fill_brewer(palette = "Set2") +  # Nice color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text( size = 14,angle = 45,hjust = 1,vjust = 1),
    axis.text.y =  element_text( size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  
  labs(y = "Number of mutation", title = "HPV subtypes")+
  stat_compare_means(method = "anova", label.y = 1,label="p.format",label.x.npc = "right") 


TMB_plot3 <- ggplot(Result4, aes(x = FIGO, y = counting)) +
  geom_violin(alpha = 1, fill = "white",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  geom_violin(alpha = 0.4, fill = "gray",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  #geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(aes(col=Cancer),shape = 16, size = 3, width = 0.15, alpha = 0.7) +  # Add points for individual values
  #scale_y_log10() +  # Use log scale for better visualization
  scale_color_manual(
    values = c(
      "ADC" = "#FF6666",
      "SCC" = "#3399ff",
      "Commun" = "gray"
    )
  )+
  scale_y_log10() +  # Use log scale for better visualization
  scale_fill_brewer(palette = "Set2") +  # Nice color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text( size = 14,angle = 45,hjust = 1,vjust = 1),
    axis.text.y =  element_text( size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  labs(y = "Number of mutation", title = "FIGO cancer classification")+
  stat_compare_means(method = "anova", label.y = 1,label="p.format",label.x.npc = "right") 
?stat_compare_means


TMB_plot4 <- ggplot(Result4, aes(x = Integration, y = counting)) +
  geom_violin(alpha = 1, fill = "white",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  geom_violin(alpha = 0.4, fill = "gray",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  #geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(aes(col=Cancer),shape = 16, size = 3, width = 0.15, alpha = 0.7) +  # Add points for individual values
  #scale_y_log10() +  # Use log scale for better visualization
  scale_color_manual(
    values = c(
      "ADC" = "#FF6666",
      "SCC" = "#3399ff",
      "Commun" = "gray"
    )
  )+
  scale_y_log10() +  # Use log scale for better visualization
  scale_fill_brewer(palette = "Set2") +  # Nice color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text( size = 14,angle = 45,hjust = 1,vjust = 1),
    axis.text.y =  element_text( size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  labs(y = "Number of mutation", title = "Integration Statut")+
  stat_compare_means(method = "t.test", label.y = 1,label="p.format",label.x.npc = "right") 



combi_row <- TMB_plot1 + TMB_plot2 + TMB_plot3 + TMB_plot4+
  plot_layout(ncol = 4)

combi_row

# Save the plot
ggsave("16_TMB_plot_by_Strain2_FULLs_select.pdf", width = 12, height = 6, dpi = 300)



##########################################  fill tmb ######################################################


common_cols <- intersect(colnames(Result3snps), colnames(Result3indels))

Result3snps55<-Result3snps%>%
  dplyr::select(common_cols)

Result3indels55<-Result3indels%>%
  dplyr::select(common_cols)

Result32<-rbind(Result3snps55,Result3indels55)

Result4<-Result32%>%
  group_by(Paire,FIGO,Cancer,Integration,type)%>%
  summarise(counting=sum(count_variant))

Result41<-Result32%>%
  group_by(Paire,FIGO,Integration,type)%>%
  summarise(counting=sum(count_variant))


Result4$Integration[Result4$Integration=="INTE_VS"]<-"INTE"
Result4$Integration[Result4$Integration=="VS"]<-"EPI"



TMB_plot1 <- ggplot(Result4, aes(x = Cancer, y = counting)) +
  geom_violin(alpha = 1, fill = "white",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  geom_violin(alpha = 0.4, aes(fill = Cancer),col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +  # Violin plot for distributione
  #geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(shape = 16, size = 3, width = 0.15, alpha = 0.7) +
  scale_y_log10() +
  #scale_fill_brewer(palette = "Set2") +
  
  scale_fill_manual(
    values = c(
      "ADC" = "#FF6666",
      "SCC" = "#3399ff",
      "Commun" = "gray"
    )
  )+
  
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  labs(y = "Number of mutation", title = "Integration Statut") +
  
  # 🟡 Add statistical test brackets
  stat_compare_means(method = "t.test", label.y = 1,label="p.format",label.x.npc = "right") 

TMB_plot1

Result4$type[Result4$type=="HPV58"]<-"Other HPVs"
Result4$type[Result4$type=="HPV45"]<-"Other HPVs"


TMB_plot2 <- ggplot(Result4, aes(x = type, y = counting)) +
  geom_violin(alpha = 1, fill = "white",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  geom_violin(alpha = 0.4, fill = "gray",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  #geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(aes(col=Cancer),shape = 16, size = 3, width = 0.15, alpha = 0.7) +  # Add points for individual values
  #scale_y_log10() +  # Use log scale for better visualization
  scale_color_manual(
    values = c(
      "ADC" = "#FF6666",
      "SCC" = "#3399ff",
      "Commun" = "grey"
    )
  )+
  scale_y_log10() +  # Use log scale for better visualization
  scale_fill_brewer(palette = "Set2") +  # Nice color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text( size = 14,angle = 45,hjust = 1,vjust = 1),
    axis.text.y =  element_text( size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  
  labs(y = "Number of mutation", title = "HPV subtypes")+
  stat_compare_means(method = "anova", label.y = 1,label="p.format",label.x.npc = "right") 


TMB_plot3 <- ggplot(Result4, aes(x = FIGO, y = counting)) +
  geom_violin(alpha = 1, fill = "white",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  geom_violin(alpha = 0.4, fill = "gray",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  #geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(aes(col=Cancer),shape = 16, size = 3, width = 0.15, alpha = 0.7) +  # Add points for individual values
  #scale_y_log10() +  # Use log scale for better visualization
  scale_color_manual(
    values = c(
      "ADC" = "#FF6666",
      "SCC" = "#3399ff",
      "Commun" = "gray"
    )
  )+
  scale_y_log10() +  # Use log scale for better visualization
  scale_fill_brewer(palette = "Set2") +  # Nice color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text( size = 14,angle = 45,hjust = 1,vjust = 1),
    axis.text.y =  element_text( size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  labs(y = "Number of mutation", title = "FIGO cancer classification")+
  stat_compare_means(method = "anova", label.y = 1,label="p.format",label.x.npc = "right") 
?stat_compare_means


TMB_plot4 <- ggplot(Result4, aes(x = Integration, y = counting)) +
  geom_violin(alpha = 1, fill = "white",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  geom_violin(alpha = 0.4, fill = "gray",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  #geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(aes(col=Cancer),shape = 16, size = 3, width = 0.15, alpha = 0.7) +  # Add points for individual values
  #scale_y_log10() +  # Use log scale for better visualization
  scale_color_manual(
    values = c(
      "ADC" = "#FF6666",
      "SCC" = "#3399ff",
      "Commun" = "gray"
    )
  )+
  scale_y_log10() +  # Use log scale for better visualization
  scale_fill_brewer(palette = "Set2") +  # Nice color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text( size = 14,angle = 45,hjust = 1,vjust = 1),
    axis.text.y =  element_text( size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  labs(y = "Number of mutation", title = "Integration Statut")+
  stat_compare_means(method = "t.test", label.y = 1,label="p.format",label.x.npc = "right") 



combi_row <- TMB_plot1 + TMB_plot2 + TMB_plot3 + TMB_plot4+
  plot_layout(ncol = 4)

combi_row

# Save the plot
ggsave("16_TMB_plot_by_Strain2_FULLs_select.pdf", width = 12, height = 6, dpi = 300)









############################################# Mutational signature ##################################################
#################################

library(reshape2)
library(RColorBrewer)

ADC_Sample_correspondances_indels<-read.delim(file = "./MAFs_ADC/Assignement_SNPs/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",header = T)
SCC_Sample_correspondances_indels<-read.delim(file = "./MAFs_SCC/Assignement_SNPs/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",header = T)
COM_Sample_correspondances_indels<-read.delim(file = "./MAFs_COMMUN/Assignement_SNPs/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",header = T)

# ADC_Sample_correspondances_indels<-read.delim(file = "./back_upaOUT/MAFs_ADC/Assignement_SNPs/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",header = T)
# SCC_Sample_correspondances_indels<-read.delim(file = "./back_upaOUT/MAFs_SCC/Assignement_SNPs/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",header = T)
# COM_Sample_correspondances_indels<-read.delim(file = "./back_upaOUT/MAFs_COMMUN/Assignement_SNPs/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",header = T)

ADC_melt_sample<-reshape2::melt(ADC_Sample_correspondances_indels,id = c("Samples"))
SCC_melt_sample<-reshape2::melt(SCC_Sample_correspondances_indels,id = c("Samples"))
COM_melt_sample<-reshape2::melt(COM_Sample_correspondances_indels,id = c("Samples"))

ADC_melt_sample$cancer<-"ADC"
SCC_melt_sample$cancer<-"SCC"
COM_melt_sample$cancer<-"COMMUN"

melt_sample<-rbind(ADC_melt_sample,SCC_melt_sample,COM_melt_sample)

melt_sample2<-melt_sample[-grep("INDELs",melt_sample$Samples),]

melt_sample3<-melt_sample2%>%
  group_by(cancer,Samples)%>%
  mutate(full=sum(value))%>%
  filter(full>500)

melt_sample3$prop<-(melt_sample3$value/melt_sample3$full)*100
melt_sample3$Paire<-gsub("Paire\\_(\\d+)\\_.+","\\1",melt_sample3$Samples)
datas_cliniques<-read.delim(file = "./Datas_cliniques.txt", header = T, sep = "\t")
melt_sample4<- merge(x=datas_cliniques, y= melt_sample3, by = "Paire" )

melt_sample5<-melt_sample4%>%
  filter(prop>0.5)%>%
  group_by(cancer,Samples)%>%
  mutate(full2=sum(value))

values = brewer.pal(n = length(unique(melt_sample4$variable)), name = "Set3")

melt_sample5$prop2<-(melt_sample5$value/melt_sample5$full2)*100
melt_sample5$prop2
ggplot(data = melt_sample5, aes(x= reorder(New,New), y= prop2,fill = variable))+
  geom_bar(color="black",stat="identity")+
  # scale_x_discrete(position="top")+
  
  
  scale_fill_manual(  values = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB", "#E78AC3","#A6D854", "#FFD92F")) +
  
  
  theme(axis.text.x = element_text(angle = 45,hjust= 0.5, vjust = 0.5, size=12),
        panel.background = element_blank(),axis.title =  element_blank())+
  facet_grid(.~cancer,scale="free_x",space="free_x")  


myOutFile1 <- paste0("propRAWFig_SIg_SNPs_nmut.pdf")
ggsave(myOutFile1, width=15, height=10, limitsize = FALSE, device='pdf', dpi=300)

ggplot(data = melt_sample5, aes(x= reorder(New,New), y= value,fill = variable))+
  geom_bar(color="black",stat="identity")+
  # scale_x_discrete(position="top")+
 # scale_fill_manual(  values = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB", "#E78AC3","#A6D854", "#FFD92F")) +
  
  
  theme(axis.text.x = element_text(angle = 45,hjust= 0.5, vjust = 0.5, size=12),
        panel.background = element_blank(),axis.title =  element_blank())+
  facet_grid(.~cancer,scale="free_x",space="free_x")  


myOutFile1 <- paste0("RAWFig_SIg_SNPs_nmut.pdf")
ggsave(myOutFile1, width=15, height=10, limitsize = FALSE, device='pdf', dpi=300)

#####


ADC_Sample_correspondances_indels<-read.delim(file = "./MAFs_ADC/Assignement_INDELS/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",header = T)
SCC_Sample_correspondances_indels<-read.delim(file = "./MAFs_SCC/Assignement_INDELS/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",header = T)
COM_Sample_correspondances_indels<-read.delim(file = "./MAFs_COMMUN/Assignement_INDELS/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",header = T)

# ADC_Sample_correspondances_indels<-read.delim(file = "./back_upaOUT/MAFs_ADC/Assignement_INDELS/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",header = T)
# SCC_Sample_correspondances_indels<-read.delim(file = "./back_upaOUT/MAFs_SCC/Assignement_INDELS/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",header = T)
# COM_Sample_correspondances_indels<-read.delim(file = "./back_upaOUT/MAFs_COMMUN/Assignement_INDELS/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",header = T)

ADC_melt_sample<-reshape2::melt(ADC_Sample_correspondances_indels,id = c("Samples"))
SCC_melt_sample<-reshape2::melt(SCC_Sample_correspondances_indels,id = c("Samples"))
COM_melt_sample<-reshape2::melt(COM_Sample_correspondances_indels,id = c("Samples"))

ADC_melt_sample$cancer<-"ADC"
SCC_melt_sample$cancer<-"SCC"
COM_melt_sample$cancer<-"COMMUN"

melt_sample<-rbind(ADC_melt_sample,SCC_melt_sample,COM_melt_sample)

#melt_sample2<-melt_sample[-grep("INDELs",melt_sample$Samples),]

melt_sample3<-melt_sample%>%
  group_by(cancer,Samples)%>%
  mutate(full=sum(value))%>%
  filter(full>50)

melt_sample3$prop<-(melt_sample3$value/melt_sample3$full)*100
melt_sample3$Paire<-gsub("Paire\\_(\\d+)\\_.+","\\1",melt_sample3$Samples)
datas_cliniques<-read.delim(file = "./Datas_cliniques.txt", header = T, sep = "\t")
melt_sample4<- merge(x=datas_cliniques, y= melt_sample3, by = "Paire" )

melt_sample5<-melt_sample4%>%
  filter(prop>5)%>%
  group_by(cancer,Samples)%>%
  mutate(full2=sum(value))
melt_sample5$prop2<-(melt_sample5$value/melt_sample5$full2)*100

ggplot(data = melt_sample5, aes(x= reorder(New,New), y= value,fill = variable))+
  geom_bar(color="black",stat="identity")+
  # scale_x_discrete(position="top")+
  scale_fill_manual(values = brewer.pal(n = length(unique(melt_sample4$variable)), name = "Set2")) +
  
  
  theme(axis.text.x = element_text(angle = 45,hjust= 0.5, vjust = 0.5, size=12),
        panel.background = element_blank(),axis.title =  element_blank())+
  facet_grid(.~cancer,scale="free_x",space="free_x")  


myOutFile1 <- paste0("RAWFig_SIg_id_nmut.pdf")
ggsave(myOutFile1, width=15, height=10, limitsize = FALSE, device='pdf', dpi=300)


ggplot(data = melt_sample5, aes(x= reorder(New,New), y= prop2,fill = variable))+
  geom_bar(color="black",stat="identity")+
  # scale_x_discrete(position="top")+
  scale_fill_manual(values = brewer.pal(n = length(unique(melt_sample4$variable)), name = "Set2")) +
  
  
  theme(axis.text.x = element_text(angle = 45,hjust= 0.5, vjust = 0.5, size=12),
        panel.background = element_blank(),axis.title =  element_blank())+
  facet_grid(.~cancer,scale="free_x",space="free_x")  


myOutFile1 <- paste0("propRAWFig_SIg_id_nmut.pdf")
ggsave(myOutFile1, width=15, height=10, limitsize = FALSE, device='pdf', dpi=300)















# ==================== 4. Tous les gènes KEGG humains par pathway ====================
message("Téléchargement de tous les pathways KEGG humains...")

kegg_pathways <- keggList("pathway", "hsa")
pathway_df <- data.frame(
  PathwayID    = names(kegg_pathways),
  Description  = gsub(" - Homo sapiens \\(human\\)", "", kegg_pathways),
  stringsAsFactors = FALSE
)

message("Téléchargement des gènes par pathway... (cela peut prendre plusieurs minutes)")

get_entrez_genes <- function(pid) {
  entry <- tryCatch(keggGet(pid)[[1]], error = function(e) return(NULL))
  if (is.null(entry) || is.null(entry$GENE)) return(NA)
  gene_ids <- entry$GENE[seq(1, length(entry$GENE), by = 2)]
  return(gene_ids)
}

entrez_list <- lapply(pathway_df$PathwayID, get_entrez_genes)
names(entrez_list) <- pathway_df$Description

# Conversion ENTREZ → SYMBOL
message("Conversion ENTREZ → SYMBOL pour chaque pathway...")

symbol_list <- lapply(entrez_list, function(entrez_ids) {
  if (is.na(entrez_ids[1])) return(NA)
  mapIds(
    org.Hs.eg.db,
    keys     = entrez_ids,
    column   = "SYMBOL",
    keytype  = "ENTREZID",
    multiVals = "first"
  )
})

# Construction du data.frame final
pathway_gene_df <- do.call(rbind, lapply(names(symbol_list), function(pw) {
  genes <- symbol_list[[pw]]
  if (is.na(genes[1])) return(NULL)
  data.frame(
    Pathway     = pw,
    GeneSymbol  = genes,
    stringsAsFactors = FALSE
  )
}))

# Nettoyage des doublons
pathway_gene_df <- unique(pathway_gene_df)



######################################################  ######################################################  ######################################################  
######################################################  FULL mutations  ##################################################################
######################################################  ######################################################  ###################################################### 

Loas_files3<-read.delim(file = "Mutect2_VCF0.8.maf", header = T, sep = "\t")
# Loas_files3<-subset(Loas_files3,
#                     CADD_PHRED_1.7 > 20 &
#                       # REVEL_score > 0.5 &
#                       MutationTaster_pred == "D" &
#                       Variant_Classification %in% c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site"))

Loas_files_INDELs_mutect2<-Loas_files3[Loas_files3$Variant_Type=="INS"|Loas_files3$Variant_Type=="DEL",]
Loas_files_SNPs<-Loas_files3[Loas_files3$Variant_Type=="SNP"|Loas_files3$Variant_Type=="TNP"|Loas_files3$Variant_Type=="DNP",]



datas_cliniques<-read.delim(file = "./Datas_cliniques.txt", header = T, sep = "\t")
CODING_snps<-Loas_files_SNPs[Loas_files_SNPs$Variant_Classification=="Missense_Mutation"|Loas_files_SNPs$Variant_Classification=="Nonsense_Mutation"|Loas_files_SNPs$Variant_Classification=="Nonstop_Mutation"|Loas_files_SNPs$Variant_Classification=="Translation_Start_Site",]
#CODING_snps<- merge(x=datas_cliniques, y= CODING_snps, by = "Paire" )
COMBI_snps<-CODING_snps[,c('New','Cancer','Hugo_Symbol','Variant_Classification','Chromosome','Start_Position','Reference_Allele','End_Position','Tumor_Seq_Allele2')]

CODING_INDELs<-Loas_files_INDELs_mutect2[!Loas_files_INDELs_mutect2$Variant_Classification=="Non exonic",]
CODING_INDELs
#CODING_INDELs<- merge(x=datas_cliniques, y= CODING_INDELs, by = "Paire" )
COMBI_indels<-CODING_INDELs[,c('New','Cancer','Hugo_Symbol','Variant_Type','Chromosome','Start_Position','Reference_Allele','End_Position','Tumor_Seq_Allele2')]

########################################################"
gene_refs<-read.delim("../CNV_analyse/output_CDS4_Refs_genes_HG38.full.txt",header = F)

datas_cliniques <- read.delim(file = "./Datas_cliniques.txt", header = TRUE, sep = "\t")


Ful_datas<-read.delim(file = "./CNVs/full_gene_focal_event2.txt", header = T, sep = "\t")
Ful_datas$Paire_num<-gsub("Paire_(\\d+)","\\1",Ful_datas$Paire_ID)
Ful_datas$Paire<-as.numeric(Ful_datas$Paire_num)
Ful_datas$TYPE<-Ful_datas$Cancer
Ful_datas$gene<-Ful_datas$Gene
Ful_datas$sig<-"AMP"
Ful_datas$sig[which(Ful_datas$Seg_mean<0)]<-"DEL"
Ful_datas$Paire<-gsub("17","15",Ful_datas$Paire)
Ful_datas <- merge(datas_cliniques, Ful_datas, by.x = "Paire", by.y = "Paire")

Ful_datas2<-Ful_datas[Ful_datas$Seg_mean > 0.6 | Ful_datas$Seg_mean < -1,]


COMBI_cnvs<-Ful_datas2[,c('New','TYPE','gene','sig')]
COMBI_cnvs$Chromosome<-"."
COMBI_cnvs$Start_Position<-""
COMBI_cnvs$Reference_Allele<-"."
COMBI_cnvs$End_Position<-"."
COMBI_cnvs$Tumor_Seq_Allele2<-"."

colnames(COMBI_snps)<-colnames(COMBI_cnvs)
colnames(COMBI_indels)<-colnames(COMBI_cnvs)

COMBI_snps$Mut<-"SNPs"
COMBI_cnvs$Mut<-"CNVs"
COMBI_indels$Mut<-"INDELs"

COMBI<-rbind(COMBI_cnvs,COMBI_snps,COMBI_indels)

colnames(COMBI)<-c("New","Cancer","Gene","sig","Chromosome","Start_Position","Reference_Allele","End_Position","Tumor_Seq_Allele2","Mut")

#install.packages("HGNChelper")   # si pas encore installé
library(HGNChelper)

corrected <- checkGeneSymbols(COMBI$Gene)
COMBI$Gene<-corrected$Suggested.Symbo

####################################

COMBI_ADC<-COMBI%>%
  group_by(New,Gene,sig,Chromosome,Start_Position,Reference_Allele,End_Position)%>%
  mutate(test=n_distinct(Cancer))%>%
  filter(test<2)%>%
  filter(Cancer=="ADC")

COMBI_SCC<-COMBI%>%
  group_by(New,Gene,sig,Chromosome,Start_Position,Reference_Allele,End_Position)%>%
  mutate(test=n_distinct(Cancer))%>%
  filter(test<2)%>%
  filter(Cancer=="SCC")

COMBI_COM<-COMBI%>%
  group_by(New,Gene,sig,Chromosome,Start_Position,Reference_Allele,End_Position)%>%
  mutate(test=n_distinct(Cancer))%>%
  filter(test==2)%>%
  filter(Cancer=="ADC")

####################################
COMBI_analyze<-COMBI%>%
  group_by(New,Gene,Cancer)%>%
  mutate(multi=n(),cross2=n_distinct(Mut))%>%
  ungroup()%>%
  group_by(Gene)%>%
  mutate(cross=n_distinct(New))%>%
  ungroup()%>%
  group_by(Gene,Mut)%>%
  mutate(cross3=n_distinct(New))%>%
  ungroup()%>%
  mutate(FillColor = ifelse(cross2 > 1, "Multi", Mut))%>%
  arrange(cross)

COMBI_COM_analyze<-COMBI_COM%>%
  group_by(New,Gene,Cancer)%>%
  mutate(multi=n(),cross2=n_distinct(Mut))%>%
  ungroup()%>%
  group_by(Gene)%>%
  mutate(cross=n_distinct(New))%>%
  ungroup()%>%
  group_by(Gene,Mut)%>%
  mutate(cross3=n_distinct(New))%>%
  ungroup()%>%
  mutate(FillColor = ifelse(cross2 > 1, "Multi", Mut))%>%
  arrange(cross)


COMBI_ADC_analyze<-COMBI_ADC%>%
  group_by(New,Gene,Cancer)%>%
  mutate(multi=n(),cross2=n_distinct(Mut))%>%
  ungroup()%>%
  group_by(Gene)%>%
  mutate(cross=n_distinct(New))%>%
  ungroup()%>%
  group_by(Gene,Mut)%>%
  mutate(cross3=n_distinct(New))%>%
  ungroup()%>%
  mutate(FillColor = ifelse(cross2 > 1, "Multi", Mut))%>%
  arrange(cross)

COMBI_SCC_analyze<-COMBI_SCC%>%
  group_by(New,Gene,Cancer)%>%
  mutate(multi=n(),cross2=n_distinct(Mut))%>%
  ungroup()%>%
  group_by(Gene)%>%
  mutate(cross=n_distinct(New))%>%
  ungroup()%>%
  group_by(Gene,Mut)%>%
  mutate(cross3=n_distinct(New))%>%
  ungroup()%>%
  mutate(FillColor = ifelse(cross2 > 1, "Multi", Mut))%>%
  arrange(cross)

####################


COMBI_analyze$PROP<-(COMBI_analyze$cross/12)*100
COMBI_COM_analyze$PROP<-(COMBI_COM_analyze$cross/12)*100
COMBI_ADC_analyze$PROP<-(COMBI_ADC_analyze$cross/12)*100
COMBI_SCC_analyze$PROP<-(COMBI_SCC_analyze$cross/12)*100

#####

COMBI_plot <- COMBI_analyze %>%
  mutate(FillColor = ifelse(cross2 > 1, "Multi", Mut))%>%
  arrange(PROP)

COMBI_plot2<-unique(COMBI_plot[,c("New","Cancer","Gene","PROP","FillColor")])
COMBI_plot3<-unique(COMBI_plot[,c("Gene","PROP")])
COMBI_analyze2<-COMBI_analyze%>%
  dplyr::select(Cancer,Mut)

#####

COMBI_plot <- COMBI_analyze %>%
  mutate(FillColor = ifelse(cross2 > 1, "Multi", Mut))%>%
  arrange(desc(cross))

COMBI_plot2<-unique(COMBI_plot[,c("New","Cancer","Gene","PROP","FillColor","cross")])
COMBI_plot3<-unique(COMBI_plot[,c("Gene","cross3","Mut","cross")])
COMBI_plot_selection<-unique(COMBI_plot[,c("Gene","cross")])

x0 <- ggplot(COMBI_analyze, aes(y = factor(Cancer), fill = Mut)) +
  geom_bar(position = "stack",col="white") +
  theme_minimal() +
  scale_fill_manual(
    values = c("CNVs" = "#A6C8F4", "SNPs" = "#90EE90", "Multi" = "black", "INDELs" = "#F4A6A6"), 
    na.value = "white"
  ) +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank()
  ) +
  facet_grid(New ~ ., switch = "x") +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "none" 
  )

x0
# # 1) Build a single canonical gene ordering from COMBI_plot3 (adjust asc/desc below if you prefer)
# gene_levels <- COMBI_plot3 %>%
#   distinct(Gene, .keep_all = TRUE) %>%
#   arrange(cross) %>%    # ascending ORDER: smallest PROP -> bottom ; largest -> top
#   pull(Gene)%>%head(50)
COMBI_plot3$cross

gene_levels <- COMBI_plot_selection %>%
  #distinct(Gene, .keep_all = TRUE) %>%
  #arrange(cross) %>%  
  head(50)  %>%# ascending order
  dplyr::select(Gene) 
# dplyr::select only the first 50 genes

# 
# COMBI_plot_selection2SCC<-COMBI_plot_selection[!COMBI_plot_selection$Gene%in%COMBI_ADC_analyze$Gene,]
# gene_levels <- COMBI_plot_selection2SCC %>%
#   filter(cross > 2) %>%
#   dplyr::select(Gene)

gene_levels$Gene
COMBI_plot2<-COMBI_plot2[COMBI_plot2$Gene%in%gene_levels$Gene,]
COMBI_plot3<-COMBI_plot3[COMBI_plot3$Gene%in%gene_levels$Gene,]
# COMBI_plot3
# COMBI_plot2 <- COMBI_plot2 %>%
#   mutate(Gene = factor(Gene, levels = gene_levels))
# 
# COMBI_plot3 <- COMBI_plot3 %>%
#   mutate(Gene = factor(Gene, levels = gene_levels))



# 2) Rebuild x1 and x2 so they use y = Gene (the factor) instead of reorder()
# 2) Rebuild x1 and x2 so they use y = Gene (the factor) instead of reorder()
x1_clean <- ggplot(COMBI_plot3, aes(x = cross, y = reorder(COMBI_plot3$Gene,(as.numeric(COMBI_plot3$cross))))) +
  geom_col(width = 1, fill = "black",col="white") +          # change fill if you want a different color
  #  facet_grid(. ~ New, switch = "x") +             # same faceting as x2 so columns align
  scale_fill_manual(
    values = c("CNVs" = "#A6C8F4", "SNPs" = "#90EE90", "Multi" = "black", "INDELs" = "#F4A6A6"), 
    na.value = "white"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_blank()                 # hide duplicate gene labels (keep them in x2 if desired)
  )


# Création de la grille complète
complete_data <- expand_grid(
  Cancer = unique(COMBI_plot2$Cancer),
  Gene = unique(COMBI_plot2$Gene),
  New = unique(COMBI_plot2$New)
)

# Fusion avec les vraies données
COMBI_plot2_filled <- complete_data %>%
  left_join(COMBI_plot2, by = c("Cancer", "Gene", "New")) %>%
  mutate(x_label = paste(New, Cancer, sep = "_"))

# Heatmap avec quadrillage continu (sans facette)
x2_clean <- ggplot(COMBI_plot2_filled, aes(
  x = factor(reorder(x_label,New)),
  y = reorder(Gene, as.numeric(as.factor(Gene))),
  fill = FillColor
)) +
  geom_tile(color = "grey", size = 0.3, na.rm = FALSE) +
  scale_fill_manual(
    values = c("CNVs" = "#A6C8F4", "SNPs" = "#90EE90", "Multi" = "black", "INDELs" = "#F4A6A6"),
    na.value = "white"
  ) +
  theme_minimal() +
  labs(x = "Patient ID", y = "Gene", fill = "Mutation Type") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5,vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank()
  )


x2_clean


# 3) Prepare x0: remove its legend (we'll collect a single legend from x2_clean)
x0_clean <- x0 + theme(legend.position = "none")

# 4) Layout:
# top row = x0 | spacer
# bottom row = x2 | x1
top_row    <- x0_clean  | plot_spacer()   # spacer prevents x1 from spanning up
bottom_row <- x2_clean | x1_clean
bottom_row2<-bottom_row+plot_layout( widths = c(15, 3), guides = "collect") 
final_plot <- top_row / bottom_row2 +
  plot_layout(heights = c(1, 4), guides = "collect") &
  theme(legend.position = "left")

# show
final_plot

ggsave("FULL_Full_total_genes.pdf", plot = final_plot, width = 12, height = 12, limitsize = FALSE, device = 'pdf', dpi = 300)



gene_df <- bitr(COMBI_analyze$Gene,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)



kegg_enrich <- enrichKEGG(
  gene         = gene_df$ENTREZID,
  organism     = 'hsa',      # hsa = Homo sapiens
  pvalueCutoff = 0.05
)


kegg_annotated <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
head(as.data.frame(kegg_annotated))

pdf("FULL_mutations_barplot_kegg_enrich.pdf", width = 10, height = 28)  # Ouvre le PDF
barplot(kegg_enrich, showCategory = 331)                # Crée le plot
dev.off()                                               # Ferme le PDF

# Figure 2 : dotplot
pdf("FULL_mutations_dotplot_kegg_enrich.pdf", width = 10, height = 28)
dotplot(kegg_enrich, showCategory = 331)
dev.off()


# 
# 
# # Extraire les résultats sous forme de data.frame
# kegg_df <- as.data.frame(kegg_enrich)
# 
# # Récupérer les noms (colonne "Description")
# kegg_pathways <- kegg_df$Description
# 
# # Voir les N premiers (ici 331 si tu veux exactement ceux du barplot)
# head(kegg_pathways, 331)
# 
# 
# hpv_genes <- keggGet("hsa05165")[[1]]$GENE
# 
# 
# sel<-COMBI_SCC_analyze[COMBI_SCC_analyze$Gene%in%hpv_genes_symbols,]
# 
# 
# # Les gènes KEGG sont alternés : EntrezID, Description
# hpv_genes_symbols <- hpv_genes[seq(2, length(hpv_genes), by = 2)]
# 
# # Extraire juste le nom du gène en enlevant la description
# hpv_genes_symbols <- sapply(strsplit(hpv_genes_symbols, ";"), `[`, 1)
# 
# # Résultat final
# hpv_genes_symbols <- unique(trimws(hpv_genes_symbols))
# print(hpv_genes_symbols)
# 
# 
# hpv_genes_symbols=="COL4A6"




############################################################################################################
#####################################      ADC    ##########################################################
############################################################################################################
############################################################################################################
#####

COMBI_plot <- COMBI_ADC_analyze
COMBI_plot2<-unique(COMBI_plot[,c("New","Cancer","Gene","FillColor","cross")])
COMBI_plot3<-unique(COMBI_plot[,c("Gene","cross")])
COMBI_plot_selection<-unique(COMBI_plot[,c("Gene","cross")])

COMBI_plot4<-data.frame(New=c(12:1),VC=c(12:1))


COMBI_plot4$New <- factor(COMBI_plot4$New, levels = 12:1)
COMBI_ADC_analyze$New <- factor(COMBI_ADC_analyze$New, levels = 12:1)
COMBI_plot$New

x0 <- ggplot(COMBI_plot, aes(x = reorder(New, -New), fill = Mut)) +
  geom_bar(position = "stack", col = "white") +
  coord_flip() +  # Barres horizontales
  theme_minimal() +
  scale_x_discrete(position = "top") +  # 👈 Affiche les étiquettes de New à droite après flip
  scale_y_reverse() +  # Inverser l'axe des valeurs : 0 à droite
  scale_fill_manual(
    values = c("CNVs" = "#A6C8F4", "SNPs" = "#90EE90", "Multi" = "black", "INDELs" = "#F4A6A6"),
    na.value = "white"
  ) +
  #labs(x = "Number of mutation") +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    
    axis.text.y = element_text(angle = 0, size = 12,hjust = 0.5 ),     # Valeurs numériques (horizontal)
    axis.ticks.y = element_blank(),
    
    # Supprimer les étiquettes de New à gauche (en bas de x après flip)
    axis.text.x.bottom = element_text(angle = 0, size = 12),
    axis.ticks.x.bottom = element_line(),
    
    # Afficher les étiquettes de New à droite (en haut de x après flip)
    axis.text.x.top = element_text(angle = 0, size = 12),
    axis.ticks.x.top = element_line(),
    
    legend.position = "none"
  )

x0


# # 1) Build a single canonical gene ordering from COMBI_plot3 (adjust asc/desc below if you prefer)
# gene_levels <- COMBI_plot3 %>%
#   distinct(Gene, .keep_all = TRUE) %>%
#   arrange(cross) %>%    # ascending ORDER: smallest PROP -> bottom ; largest -> top
#   pull(Gene)%>%head(50)
COMBI_plot3$cross

library(dplyr)

#COMBI_plot_selection2ADC<-COMBI_plot_selection[!COMBI_plot_selection$Gene%in%COMBI_SCC_analyze$Gene,]
gene_levels <- COMBI_plot_selection %>%
  filter(cross > 2) %>%
  dplyr::select(Gene)

gene_levels_adc <- gene_levels$Gene

COMBI_plot2<-COMBI_plot2[COMBI_plot2$Gene%in%gene_levels$Gene,]
COMBI_plot3<-COMBI_plot3[COMBI_plot3$Gene%in%gene_levels$Gene,]
# COMBI_plot3
# COMBI_plot2 <- COMBI_plot2 %>%
#   mutate(Gene = factor(Gene, levels = gene_levels))
# 
# COMBI_plot3 <- COMBI_plot3 %>%
#   mutate(Gene = factor(Gene, levels = gene_levels))



# 2) Rebuild x1 and x2 so they use y = Gene (the factor) instead of reorder()
x1_clean <- ggplot(COMBI_plot3, aes(x = cross, y = reorder(COMBI_plot3$Gene,(-as.numeric(COMBI_plot3$cross))))) +
  geom_col(width = 1, fill = "black",col="white") +          # change fill if you want a different color
  #  facet_grid(. ~ New, switch = "x") +             # same faceting as x2 so columns align
  coord_flip() +
  scale_fill_manual(
    values = c("CNVs" = "#A6C8F4", "SNPs" = "#90EE90", "Multi" = "black", "INDELs" = "#F4A6A6"), 
    na.value = "white"
  ) +
  theme_minimal() +
  labs(y = "n Pair") +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(angle = 0, size = 12)# hide duplicate gene labels (keep them in x2 if desired)
  )

x0_clean
# Force factor levels 1–12 in both data frames
COMBI_plot4$New <- factor(COMBI_plot4$New, levels = 1:12)
COMBI_plot2$New <- factor(COMBI_plot2$New, levels = 1:12)


# Création de la grille complète
complete_data <- expand_grid(
  Cancer = unique(COMBI_plot2$Cancer),
  Gene = unique(COMBI_plot2$Gene),
  New = c("1","2","3","4","5","6","7","8","9","10","11","12")
)

#View(COMBI_plot2_filled)
# Fusion avec les vraies données
COMBI_plot2_filled <- complete_data %>%
  left_join(COMBI_plot2, by = c("Cancer", "Gene", "New")) %>%
  mutate(x_label = paste(New, Cancer, sep = "_"))

as.numeric(as.factor(COMBI_plot2_filled$Gene))

# Heatmap avec quadrillage continu (sans facette)

COMBI_plot2_filled <- COMBI_plot2_filled %>%
  mutate(
    Gene = factor(Gene, levels = unique(Gene[order(-cross,Gene)]))
  )

COMBI_plot2_filled

x2_clean <- ggplot(COMBI_plot2_filled, aes(
  x = factor(reorder(x_label,-as.numeric(New))),
  y = reorder(Gene, -cross),
  fill = FillColor
)) +
  coord_flip() +
  geom_tile(color = "grey", size = 0.3, na.rm = FALSE) +
  scale_fill_manual(
    values = c("CNVs" = "#A6C8F4", "SNPs" = "#90EE90", "Multi" = "black", "INDELs" = "#F4A6A6"),
    na.value = "white"
  ) +
  theme_minimal() +
  labs(x = "Pair", y = "Gene", fill = "Mutation Type") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,vjust = 0.95, size = 10),
    axis.text.y = element_blank(),
    panel.grid = element_blank()
  )


x2_clean

dev.off()



x2_clean
x0_clean

# 3) Prepare x0: remove its legend (we'll collect a single legend from x2_clean)
x0_clean <- x0 + theme(legend.position = "none")

# 4) Layout:
# top row = x0 | spacer
# bottom row = x2 | x1
top_row    <- plot_spacer()|x1_clean   # spacer prevents x1 from spanning up
top_row2<-top_row+plot_layout( widths = c(4, 17), guides = "collect") 
top_row2
bottom_row <- x0_clean | x2_clean
bottom_row2<-bottom_row+plot_layout( widths = c(3, 15), guides = "collect") 
bottom_row
top_row
final_plot <- top_row2 / bottom_row2 +
  plot_layout(heights = c(1, 4), guides = "collect") &
  theme(legend.position = "none")

# show
final_plot

ggsave("FULL_ADC_total_genes.pdf", plot = final_plot, width = 14, height = 6, limitsize = FALSE, device = 'pdf', dpi = 300)


#############################################
#############################################
#############################################



############################################################################################################
#####################################      SCC    ##########################################################
############################################################################################################
############################################################################################################
#####

COMBI_plot <- COMBI_SCC_analyze
COMBI_plot2<-unique(COMBI_plot[,c("New","Cancer","Gene","FillColor","cross")])
COMBI_plot3<-unique(COMBI_plot[,c("Gene","cross")])
COMBI_plot_selection<-unique(COMBI_plot[,c("Gene","cross")])

COMBI_plot4<-data.frame(New=c(1:12),VC=c(1:12))


COMBI_plot4$New <- factor(COMBI_plot4$New, levels = 1:12)
COMBI_SCC_analyze$New <- factor(COMBI_SCC_analyze$New, levels = 1:12)
COMBI_plot$Mut
library(dplyr)
library(tidyr)
library(ggplot2)

# On force New comme numérique de 1 à 12 et on complète
COMBI_plot <- COMBI_plot %>%
  mutate(New = as.numeric(New)) %>%
  group_by(New, Mut) %>%
  summarise(Mut_num = n(), .groups = "drop") %>%
  complete(
    New = 1:12,
    fill = list(Mut_num = 0)
  )

# Graphique
x0 <- ggplot(COMBI_plot, aes(x = reorder(factor(New),-New), y = Mut_num, fill = Mut)) +
  geom_col(col = "white") +
  coord_flip() +  # barres horizontales
  theme_minimal() +
  scale_x_discrete(position = "top", breaks = 1:12, labels = 1:12) +  
  scale_y_reverse() +  
  scale_fill_manual(
    values = c("CNVs" = "#A6C8F4", 
               "SNPs" = "#90EE90", 
               "Multi" = "black", 
               "INDELs" = "#F4A6A6"),
    na.value = "white"
  ) +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12, hjust = 0.5),
    axis.ticks.y = element_blank(),
    axis.text.x.bottom = element_text(size = 12),
    axis.ticks.x.bottom = element_line(),
    axis.text.x.top = element_text(size = 12),
    axis.ticks.x.top = element_line(),
    legend.position = "none"
  )

x0



# # 1) Build a single canonical gene ordering from COMBI_plot3 (adjust asc/desc below if you prefer)
# gene_levels <- COMBI_plot3 %>%
#   distinct(Gene, .keep_all = TRUE) %>%
#   arrange(cross) %>%    # ascending ORDER: smallest PROP -> bottom ; largest -> top
#   pull(Gene)%>%head(50)
COMBI_plot3$cross

library(dplyr)
COMBI_SCC_analyze

#COMBI_plot_selection2SCC<-COMBI_plot_selection[!COMBI_plot_selection$Gene%in%COMBI_SCC_analyze$Gene,]



gene_levels <- COMBI_SCC_analyze %>%
  filter(cross > 2) %>%
  dplyr::select(Gene)

gene_levels_SCC <- gene_levels$Gene

COMBI_plot2<-COMBI_plot2[COMBI_plot2$Gene%in%gene_levels$Gene,]
COMBI_plot3<-COMBI_plot3[COMBI_plot3$Gene%in%gene_levels$Gene,]
# COMBI_plot3
# COMBI_plot2 <- COMBI_plot2 %>%
#   mutate(Gene = factor(Gene, levels = gene_levels))
# 
# COMBI_plot3 <- COMBI_plot3 %>%
#   mutate(Gene = factor(Gene, levels = gene_levels))



# 2) Rebuild x1 and x2 so they use y = Gene (the factor) instead of reorder()
x1_clean <- ggplot(COMBI_plot3, aes(x = cross, y = reorder(COMBI_plot3$Gene,(-as.numeric(COMBI_plot3$cross))))) +
  geom_col(width = 1, fill = "black",col="white") +          # change fill if you want a different color
  #  facet_grid(. ~ New, switch = "x") +             # same faceting as x2 so columns align
  coord_flip() +
  scale_fill_manual(
    values = c("CNVs" = "#A6C8F4", "SNPs" = "#90EE90", "Multi" = "black", "INDELs" = "#F4A6A6"), 
    na.value = "white"
  ) +
  theme_minimal() +
  labs(y = "n Pair") +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(angle = 0, size = 12)# hide duplicate gene labels (keep them in x2 if desired)
  )

x1_clean
# Force factor levels 1–12 in both data frames
COMBI_plot4$New <- factor(COMBI_plot4$New, levels = 1:12)
COMBI_plot2$New <- factor(COMBI_plot2$New, levels = 1:12)


# Création de la grille complète
complete_data <- expand_grid(
  Cancer = unique(COMBI_plot2$Cancer),
  Gene = unique(COMBI_plot2$Gene),
  New = c("1","2","3","4","5","6","7","8","9","10","11","12")
)

#View(COMBI_plot2_filled)
# Fusion avec les vraies données
COMBI_plot2_filled <- complete_data %>%
  left_join(COMBI_plot2, by = c("Cancer", "Gene", "New")) %>%
  mutate(x_label = paste(New, Cancer, sep = "_"))

as.numeric(as.factor(COMBI_plot2_filled$Gene))

# Heatmap avec quadrillage continu (sans facette)

COMBI_plot2_filled <- COMBI_plot2_filled %>%
  mutate(
    Gene = factor(Gene, levels = unique(Gene[order(-cross,Gene)]))
  )

COMBI_plot2_filled

x2_clean <- ggplot(COMBI_plot2_filled, aes(
  x = factor(reorder(x_label,-as.numeric(New))),
  y = reorder(Gene, -cross),
  fill = FillColor
)) +
  coord_flip() +
  geom_tile(color = "grey", size = 0.3, na.rm = FALSE) +
  scale_fill_manual(
    values = c("CNVs" = "#A6C8F4", "SNPs" = "#90EE90", "Multi" = "black", "INDELs" = "#F4A6A6"),
    na.value = "white"
  ) +
  theme_minimal() +
  labs(x = "Pair", y = "Gene", fill = "Mutation Type") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,vjust = 0.95, size = 12),
    axis.text.y = element_blank(),
    panel.grid = element_blank()
  )


x2_clean

dev.off()






# 3) Prepare x0: remove its legend (we'll collect a single legend from x2_clean)
x0_clean <- x0 + theme(legend.position = "none")

# 4) Layout:
# top row = x0 | spacer
# bottom row = x2 | x1
top_row    <- plot_spacer()|x1_clean   # spacer prevents x1 from spanning up
top_row2<-top_row+plot_layout( widths = c(4, 17), guides = "collect") 
top_row2
bottom_row <- x0_clean | x2_clean
bottom_row2
bottom_row2<-bottom_row+plot_layout( widths = c(3, 15), guides = "collect") 

final_plot <- top_row2 / bottom_row2 +
  plot_layout(heights = c(1, 4), guides = "collect") &
  theme(legend.position = "none")

# show
final_plot

ggsave("FULL_SCC_total_genes.pdf", plot = final_plot, width = 10, height = 6, limitsize = FALSE, device = 'pdf', dpi = 300)
ggsave("FULL_SCC_total_genes2.pdf", plot = final_plot, width = 14, height = 6, limitsize = FALSE, device = 'pdf', dpi = 300)

COMBI_plot3_SCC<-COMBI_plot3




###########################
# ==================== 1. Conversion SYMBOL → ENTREZ ====================
message("Conversion des gènes SYMBOL en ENTREZID...")
gene_df <- bitr(COMBI_SCC_analyze$Gene,
                fromType = "SYMBOL",
                toType   = "ENTREZID",
                OrgDb    = org.Hs.eg.db)

# ==================== 2. Enrichissement KEGG ====================
message("Lancement de l'enrichissement KEGG...")
kegg_enrich <- enrichKEGG(
  gene         = gene_df$ENTREZID,
  organism     = 'hsa',
  pvalueCutoff = 0.05
)

# Rendre les résultats lisibles (ENTREZ → SYMBOL)
kegg_annotated <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_df <- as.data.frame(kegg_annotated)

# ==================== 3. Visualisation KEGG (barplot + dotplot) ====================
pdf("CADD_FULL_mutations_SCC_barplot_kegg_enrich.pdf", width = 10, height = 28)
x1_clean <- barplot(kegg_enrich, showCategory = 331)
dev.off()

pdf("CADD_FULL_mutations_SCC_dotplot_kegg_enrich.pdf", width = 14, height = 6)
x2_clean <- dotplot(kegg_enrich, showCategory = 331, orderBy = "GeneRatio", decreasing = FALSE) +
  coord_flip() +
  scale_x_reverse() +   
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


# ==================== 5. Jointure gènes enrichis/pathways ====================
message("Jointure avec les gènes CADD_COMBI...")

df_merged_SCC <- inner_join(COMBI_SCC, pathway_gene_df, by = c("Gene" = "GeneSymbol"))

# Comptage par pathway
COMBI_SCC_analyze <- df_merged_SCC %>%
  group_by(Pathway) %>%
  mutate(cross = n_distinct(New)) %>%
  ungroup() %>%
  arrange(desc(cross))

# Données pour le plot
COMBI_plot3 <- distinct(COMBI_SCC_analyze[, c("Pathway", "cross")])

# ==================== 6. Récupérer infos du KEGG enrichissement ====================
infos_generatios <- kegg_df[, c("Description", "GeneRatio", "p.adjust")]

COMBI_plot3_join <- inner_join(COMBI_plot3, infos_generatios, by = c("Pathway" = "Description")) %>%
  filter(p.adjust < 0.05) %>%
  mutate(Pathway = factor(Pathway, levels = Pathway[order(GeneRatio, decreasing = TRUE)]))

# ==================== 7. Plot final combiné ====================
x1_clean <- ggplot(COMBI_plot3_join, aes(x = cross, y = Pathway)) +
  geom_col(width = 1, fill = "black", col = "white") +
  coord_flip() +
  theme_minimal() +
  labs(y = "n Pair") +
  theme(
    axis.title   = element_blank(),
    panel.grid   = element_blank(),
    axis.text.x  = element_blank(),
    axis.text.y  = element_text(angle = 0, size = 12)
  )

bottom_row2 <- x1_clean / x2_clean
bottom_row2
ggsave("CADD_FULL_SCC_total_genes_barplot.pdf",
       plot = bottom_row2,
       width = 14, height = 6,
       limitsize = FALSE, device = 'pdf', dpi = 300)



#######################################################################################################

###########################
# ==================== 1. Conversion SYMBOL → ENTREZ ====================
message("Conversion des gènes SYMBOL en ENTREZID...")
gene_df <- bitr(COMBI_ADC_analyze$Gene,
                fromType = "SYMBOL",
                toType   = "ENTREZID",
                OrgDb    = org.Hs.eg.db)

# ==================== 2. Enrichissement KEGG ====================
message("Lancement de l'enrichissement KEGG...")
kegg_enrich <- enrichKEGG(
  gene         = gene_df$ENTREZID,
  organism     = 'hsa',
  pvalueCutoff = 0.05
)

# Rendre les résultats lisibles (ENTREZ → SYMBOL)
kegg_annotated <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_df <- as.data.frame(kegg_annotated)

# ==================== 3. Visualisation KEGG (barplot + dotplot) ====================
pdf("CADD_FULL_mutations_ADC_barplot_kegg_enrich.pdf", width = 10, height = 28)
x1_clean <- barplot(kegg_enrich, showCategory = 331)
dev.off()

pdf("CADD_FULL_mutations_ADC_dotplot_kegg_enrich.pdf", width = 14, height = 6)
x2_clean <- dotplot(kegg_enrich, showCategory = 331, orderBy = "GeneRatio", decreasing = FALSE) +
  coord_flip() +
  scale_x_reverse() +   
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



# ==================== 5. Jointure gènes enrichis/pathways ====================
message("Jointure avec les gènes CADD_COMBI...")

df_merged_ADC <- inner_join(COMBI_ADC, pathway_gene_df, by = c("Gene" = "GeneSymbol"))

# Comptage par pathway
COMBI_ADC_analyze <- df_merged_ADC %>%
  group_by(Pathway) %>%
  mutate(cross = n_distinct(New)) %>%
  ungroup() %>%
  arrange(desc(cross))

# Données pour le plot
COMBI_plot3 <- distinct(COMBI_ADC_analyze[, c("Pathway", "cross")])

# ==================== 6. Récupérer infos du KEGG enrichissement ====================
infos_generatios <- kegg_df[, c("Description", "GeneRatio", "p.adjust")]

COMBI_plot3_join <- inner_join(COMBI_plot3, infos_generatios, by = c("Pathway" = "Description")) %>%
  filter(p.adjust < 0.05) %>%
  mutate(Pathway = factor(Pathway, levels = Pathway[order(GeneRatio, decreasing = TRUE)]))

# ==================== 7. Plot final combiné ====================
x1_clean <- ggplot(COMBI_plot3_join, aes(x = cross, y = Pathway)) +
  geom_col(width = 1, fill = "black", col = "white") +
  coord_flip() +
  theme_minimal() +
  labs(y = "n Pair") +
  theme(
    axis.title   = element_blank(),
    panel.grid   = element_blank(),
    axis.text.x  = element_blank(),
    axis.text.y  = element_text(angle = 0, size = 12)
  )

bottom_row2 <- x1_clean / x2_clean
bottom_row2
ggsave("CADD_FULL_ADC_total_genes_barplot.pdf",
       plot = bottom_row2,
       width = 14, height = 6,
       limitsize = FALSE, device = 'pdf', dpi = 300)

message("✅ Script terminé.")



test_for_fulladc<-COMBI_plot3_join
#############################################################

######################################################  ######################################################  ######################################################  
######################################################  selected  mutations  ##################################################################
######################################################  ######################################################  ###################################################### 

Loas_files3<-read.delim(file = "Mutect2_VCF0.8.maf", header = T, sep = "\t")
Loas_files3<-subset(Loas_files3,
                    CADD_PHRED_1.7 > 20 ) #&
# REVEL_score > 0.5 &
#  MutationTaster_pred == "D" &
# Variant_Classification %in% c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site"))

Loas_files_INDELs_mutect2<-Loas_files3[Loas_files3$Variant_Type=="INS"|Loas_files3$Variant_Type=="DEL",]
Loas_files_SNPs<-Loas_files3[Loas_files3$Variant_Type=="SNP"|Loas_files3$Variant_Type=="TNP"|Loas_files3$Variant_Type=="DNP",]



datas_cliniques<-read.delim(file = "./Datas_cliniques.txt", header = T, sep = "\t")
CODING_snps<-Loas_files_SNPs[Loas_files_SNPs$Variant_Classification=="Missense_Mutation"|Loas_files_SNPs$Variant_Classification=="Nonsense_Mutation"|Loas_files_SNPs$Variant_Classification=="Nonstop_Mutation"|Loas_files_SNPs$Variant_Classification=="Translation_Start_Site",]
#CODING_snps<- merge(x=datas_cliniques, y= CODING_snps, by = "Paire" )
COMBI_snps<-CODING_snps[,c('New','Cancer','Hugo_Symbol','Variant_Classification','Chromosome','Start_Position','Reference_Allele','End_Position','Tumor_Seq_Allele2')]

CODING_INDELs<-Loas_files_INDELs_mutect2[!Loas_files_INDELs_mutect2$Variant_Classification=="Non exonic",]
#CODING_INDELs<- merge(x=datas_cliniques, y= CODING_INDELs, by = "Paire" )
COMBI_indels<-CODING_INDELs[,c('New','Cancer','Hugo_Symbol','Variant_Type','Chromosome','Start_Position','Reference_Allele','End_Position','Tumor_Seq_Allele2')]

########################################################"
gene_refs<-read.delim("../CNV_analyse/output_CDS4_Refs_genes_HG38.full.txt",header = F)

datas_cliniques <- read.delim(file = "./Datas_cliniques.txt", header = TRUE, sep = "\t")


Ful_datas<-read.delim(file = "./CNVs/full_gene_focal_event2.txt", header = T, sep = "\t")
Ful_datas$Paire_num<-gsub("Paire_(\\d+)","\\1",Ful_datas$Paire_ID)
Ful_datas$Paire<-as.numeric(Ful_datas$Paire_num)
Ful_datas$TYPE<-Ful_datas$Cancer
Ful_datas$gene<-Ful_datas$Gene
Ful_datas$sig<-"AMP"
Ful_datas$sig[which(Ful_datas$Seg_mean<0)]<-"DEL"
Ful_datas$Paire<-gsub("17","15",Ful_datas$Paire)
Ful_datas <- merge(datas_cliniques, Ful_datas, by.x = "Paire", by.y = "Paire")

Ful_datas2<-Ful_datas[Ful_datas$Seg_mean > 0.6 | Ful_datas$Seg_mean < -1,]


COMBI_cnvs<-Ful_datas2[,c('New','TYPE','gene','sig')]
COMBI_cnvs$Chromosome<-"."
COMBI_cnvs$Start_Position<-""
COMBI_cnvs$Reference_Allele<-"."
COMBI_cnvs$End_Position<-"."
COMBI_cnvs$Tumor_Seq_Allele2<-"."

colnames(COMBI_snps)<-colnames(COMBI_cnvs)
colnames(COMBI_indels)<-colnames(COMBI_cnvs)

COMBI_snps$Mut<-"SNPs"
COMBI_cnvs$Mut<-"CNVs"
COMBI_indels$Mut<-"INDELs"

COMBI<-rbind(COMBI_cnvs,COMBI_snps,COMBI_indels)

colnames(COMBI)<-c("New","Cancer","Gene","sig","Chromosome","Start_Position","Reference_Allele","End_Position","Tumor_Seq_Allele2","Mut")

#install.packages("HGNChelper")   # si pas encore installé
library(HGNChelper)

corrected <- checkGeneSymbols(COMBI$Gene)
COMBI$Gene<-corrected$Suggested.Symbol
####################################

COMBI_ADC<-COMBI%>%
  group_by(New,Gene,sig,Chromosome,Start_Position,Reference_Allele,End_Position)%>%
  mutate(test=n_distinct(Cancer))%>%
  filter(test<2)%>%
  filter(Cancer=="ADC")

COMBI_SCC<-COMBI%>%
  group_by(New,Gene,sig,Chromosome,Start_Position,Reference_Allele,End_Position)%>%
  mutate(test=n_distinct(Cancer))%>%
  filter(test<2)%>%
  filter(Cancer=="SCC")

COMBI_COM<-COMBI%>%
  group_by(New,Gene,sig,Chromosome,Start_Position,Reference_Allele,End_Position)%>%
  mutate(test=n_distinct(Cancer))%>%
  filter(test==2)%>%
  filter(Cancer=="ADC")

####################################
COMBI_analyze<-COMBI%>%
  group_by(New,Gene,Cancer)%>%
  mutate(multi=n(),cross2=n_distinct(Mut))%>%
  ungroup()%>%
  group_by(Gene)%>%
  mutate(cross=n_distinct(New))%>%
  ungroup()%>%
  group_by(Gene,Mut)%>%
  mutate(cross3=n_distinct(New))%>%
  ungroup()%>%
  mutate(FillColor = ifelse(cross2 > 1, "Multi", Mut))%>%
  arrange(cross)

COMBI_COM_analyze<-COMBI_COM%>%
  group_by(New,Gene,Cancer)%>%
  mutate(multi=n(),cross2=n_distinct(Mut))%>%
  ungroup()%>%
  group_by(Gene)%>%
  mutate(cross=n_distinct(New))%>%
  ungroup()%>%
  group_by(Gene,Mut)%>%
  mutate(cross3=n_distinct(New))%>%
  ungroup()%>%
  mutate(FillColor = ifelse(cross2 > 1, "Multi", Mut))%>%
  arrange(cross)


COMBI_ADC_analyze<-COMBI_ADC%>%
  group_by(New,Gene,Cancer)%>%
  mutate(multi=n(),cross2=n_distinct(Mut))%>%
  ungroup()%>%
  group_by(Gene)%>%
  mutate(cross=n_distinct(New))%>%
  ungroup()%>%
  group_by(Gene,Mut)%>%
  mutate(cross3=n_distinct(New))%>%
  ungroup()%>%
  mutate(FillColor = ifelse(cross2 > 1, "Multi", Mut))%>%
  arrange(cross)

COMBI_SCC_analyze<-COMBI_SCC%>%
  group_by(New,Gene,Cancer)%>%
  mutate(multi=n(),cross2=n_distinct(Mut))%>%
  ungroup()%>%
  group_by(Gene)%>%
  mutate(cross=n_distinct(New))%>%
  ungroup()%>%
  group_by(Gene,Mut)%>%
  mutate(cross3=n_distinct(New))%>%
  ungroup()%>%
  mutate(FillColor = ifelse(cross2 > 1, "Multi", Mut))%>%
  arrange(cross)

####################


COMBI_analyze$PROP<-(COMBI_analyze$cross/12)*100
COMBI_COM_analyze$PROP<-(COMBI_COM_analyze$cross/12)*100
COMBI_ADC_analyze$PROP<-(COMBI_ADC_analyze$cross/12)*100
COMBI_SCC_analyze$PROP<-(COMBI_SCC_analyze$cross/12)*100

#####

COMBI_plot <- COMBI_analyze %>%
  mutate(FillColor = ifelse(cross2 > 1, "Multi", Mut))%>%
  arrange(PROP)

COMBI_plot2<-unique(COMBI_plot[,c("New","Cancer","Gene","PROP","FillColor")])
COMBI_plot3<-unique(COMBI_plot[,c("Gene","PROP")])
COMBI_analyze2<-COMBI_analyze%>%
  dplyr::select(Cancer,Mut)

#####

COMBI_plot <- COMBI_analyze %>%
  mutate(FillColor = ifelse(cross2 > 1, "Multi", Mut))%>%
  arrange(desc(cross))

COMBI_plot2<-unique(COMBI_plot[,c("New","Cancer","Gene","PROP","FillColor","cross")])
COMBI_plot3<-unique(COMBI_plot[,c("Gene","cross3","Mut","cross")])
COMBI_plot_selection<-unique(COMBI_plot[,c("Gene","cross")])

x0 <- ggplot(COMBI_analyze, aes(y = factor(Cancer), fill = Mut)) +
  geom_bar(position = "stack",col="white") +
  theme_minimal() +
  scale_fill_manual(
    values = c("CNVs" = "#A6C8F4", "SNPs" = "#90EE90", "Multi" = "black", "INDELs" = "#F4A6A6"), 
    na.value = "white"
  ) +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank()
  ) +
  facet_grid(New ~ ., switch = "x") +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "none" 
  )

x0
# # 1) Build a single canonical gene ordering from COMBI_plot3 (adjust asc/desc below if you prefer)
# gene_levels <- COMBI_plot3 %>%
#   distinct(Gene, .keep_all = TRUE) %>%
#   arrange(cross) %>%    # ascending ORDER: smallest PROP -> bottom ; largest -> top
#   pull(Gene)%>%head(50)
COMBI_plot3$cross

gene_levels <- COMBI_plot_selection %>%
  #distinct(Gene, .keep_all = TRUE) %>%
  #arrange(cross) %>%  
  head(50)  %>%# ascending order
  dplyr::select(Gene) 
# dplyr::select only the first 50 genes

# 
# COMBI_plot_selection2SCC<-COMBI_plot_selection[!COMBI_plot_selection$Gene%in%COMBI_ADC_analyze$Gene,]
# gene_levels <- COMBI_plot_selection2SCC %>%
#   filter(cross > 2) %>%
#   dplyr::select(Gene)

gene_levels$Gene
COMBI_plot2<-COMBI_plot2[COMBI_plot2$Gene%in%gene_levels$Gene,]
COMBI_plot3<-COMBI_plot3[COMBI_plot3$Gene%in%gene_levels$Gene,]
# COMBI_plot3
# COMBI_plot2 <- COMBI_plot2 %>%
#   mutate(Gene = factor(Gene, levels = gene_levels))
# 
# COMBI_plot3 <- COMBI_plot3 %>%
#   mutate(Gene = factor(Gene, levels = gene_levels))



# 2) Rebuild x1 and x2 so they use y = Gene (the factor) instead of reorder()
# 2) Rebuild x1 and x2 so they use y = Gene (the factor) instead of reorder()
x1_clean <- ggplot(COMBI_plot3, aes(x = cross, y = reorder(COMBI_plot3$Gene,(as.numeric(COMBI_plot3$cross))))) +
  geom_col(width = 1, fill = "black",col="white") +          # change fill if you want a different color
  #  facet_grid(. ~ New, switch = "x") +             # same faceting as x2 so columns align
  scale_fill_manual(
    values = c("CNVs" = "#A6C8F4", "SNPs" = "#90EE90", "Multi" = "black", "INDELs" = "#F4A6A6"), 
    na.value = "white"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_blank()                 # hide duplicate gene labels (keep them in x2 if desired)
  )


# Création de la grille complète
complete_data <- expand_grid(
  Cancer = unique(COMBI_plot2$Cancer),
  Gene = unique(COMBI_plot2$Gene),
  New = unique(COMBI_plot2$New)
)

# Fusion avec les vraies données
COMBI_plot2_filled <- complete_data %>%
  left_join(COMBI_plot2, by = c("Cancer", "Gene", "New")) %>%
  mutate(x_label = paste(New, Cancer, sep = "_"))

# Heatmap avec quadrillage continu (sans facette)
x2_clean <- ggplot(COMBI_plot2_filled, aes(
  x = factor(reorder(x_label,New)),
  y = reorder(Gene, as.numeric(as.factor(Gene))),
  fill = FillColor
)) +
  geom_tile(color = "grey", size = 0.3, na.rm = FALSE) +
  scale_fill_manual(
    values = c("CNVs" = "#A6C8F4", "SNPs" = "#90EE90", "Multi" = "black", "INDELs" = "#F4A6A6"),
    na.value = "white"
  ) +
  theme_minimal() +
  labs(x = "Patient ID", y = "Gene", fill = "Mutation Type") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5,vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank()
  )


x2_clean


# 3) Prepare x0: remove its legend (we'll collect a single legend from x2_clean)
x0_clean <- x0 + theme(legend.position = "none")

# 4) Layout:
# top row = x0 | spacer
# bottom row = x2 | x1
top_row    <- x0_clean  | plot_spacer()   # spacer prevents x1 from spanning up
bottom_row <- x2_clean | x1_clean
bottom_row2<-bottom_row+plot_layout( widths = c(15, 3), guides = "collect") 
final_plot <- top_row / bottom_row2 +
  plot_layout(heights = c(1, 4), guides = "collect") &
  theme(legend.position = "left")

# show
final_plot

ggsave("FULL_Full_total_genes.pdf", plot = final_plot, width = 12, height = 12, limitsize = FALSE, device = 'pdf', dpi = 300)



gene_df <- bitr(COMBI_analyze$Gene,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)



kegg_enrich <- enrichKEGG(
  gene         = gene_df$ENTREZID,
  organism     = 'hsa',      # hsa = Homo sapiens
  pvalueCutoff = 0.05
)


kegg_annotated <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
head(as.data.frame(kegg_annotated))

pdf("FULL_mutations_barplot_kegg_enrich.pdf", width = 10, height = 28)  # Ouvre le PDF
barplot(kegg_enrich, showCategory = 331)                # Crée le plot
dev.off()                                               # Ferme le PDF

# Figure 2 : dotplot
pdf("FULL_mutations_dotplot_kegg_enrich.pdf", width = 10, height = 28)
dotplot(kegg_enrich, showCategory = 331)
dev.off()


# 
# 
# # Extraire les résultats sous forme de data.frame
# kegg_df <- as.data.frame(kegg_enrich)
# 
# # Récupérer les noms (colonne "Description")
# kegg_pathways <- kegg_df$Description
# 
# # Voir les N premiers (ici 331 si tu veux exactement ceux du barplot)
# head(kegg_pathways, 331)
# 
# 
# hpv_genes <- keggGet("hsa05165")[[1]]$GENE
# 
# 
# sel<-COMBI_SCC_analyze[COMBI_SCC_analyze$Gene%in%hpv_genes_symbols,]
# 
# 
# # Les gènes KEGG sont alternés : EntrezID, Description
# hpv_genes_symbols <- hpv_genes[seq(2, length(hpv_genes), by = 2)]
# 
# # Extraire juste le nom du gène en enlevant la description
# hpv_genes_symbols <- sapply(strsplit(hpv_genes_symbols, ";"), `[`, 1)
# 
# # Résultat final
# hpv_genes_symbols <- unique(trimws(hpv_genes_symbols))
# print(hpv_genes_symbols)
# 
# 
# hpv_genes_symbols=="COL4A6"




############################################################################################################
#####################################      ADC    ##########################################################
############################################################################################################
############################################################################################################
#####

COMBI_plot <- COMBI_ADC_analyze
COMBI_plot2<-unique(COMBI_plot[,c("New","Cancer","Gene","FillColor","cross")])
COMBI_plot3<-unique(COMBI_plot[,c("Gene","cross")])
COMBI_plot_selection<-unique(COMBI_plot[,c("Gene","cross")])

COMBI_plot4<-data.frame(New=c(12:1),VC=c(12:1))


COMBI_plot4$New <- factor(COMBI_plot4$New, levels = 12:1)
COMBI_ADC_analyze$New <- factor(COMBI_ADC_analyze$New, levels = 12:1)
COMBI_plot$New

x0 <- ggplot(COMBI_plot, aes(x = reorder(New, -New), fill = Mut)) +
  geom_bar(position = "stack", col = "white") +
  coord_flip() +  # Barres horizontales
  theme_minimal() +
  scale_x_discrete(position = "top") +  # 👈 Affiche les étiquettes de New à droite après flip
  scale_y_reverse() +  # Inverser l'axe des valeurs : 0 à droite
  scale_fill_manual(
    values = c("CNVs" = "#A6C8F4", "SNPs" = "#90EE90", "Multi" = "black", "INDELs" = "#F4A6A6"),
    na.value = "white"
  ) +
  #labs(x = "Number of mutation") +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    
    axis.text.y = element_text(angle = 0, size = 12,hjust = 0.5 ),     # Valeurs numériques (horizontal)
    axis.ticks.y = element_blank(),
    
    # Supprimer les étiquettes de New à gauche (en bas de x après flip)
    axis.text.x.bottom = element_text(angle = 0, size = 12),
    axis.ticks.x.bottom = element_line(),
    
    # Afficher les étiquettes de New à droite (en haut de x après flip)
    axis.text.x.top = element_text(angle = 0, size = 12),
    axis.ticks.x.top = element_line(),
    
    legend.position = "none"
  )

x0


# # 1) Build a single canonical gene ordering from COMBI_plot3 (adjust asc/desc below if you prefer)
# gene_levels <- COMBI_plot3 %>%
#   distinct(Gene, .keep_all = TRUE) %>%
#   arrange(cross) %>%    # ascending ORDER: smallest PROP -> bottom ; largest -> top
#   pull(Gene)%>%head(50)
COMBI_plot3$cross

library(dplyr)

#COMBI_plot_selection2ADC<-COMBI_plot_selection[!COMBI_plot_selection$Gene%in%COMBI_SCC_analyze$Gene,]
gene_levels <- COMBI_plot_selection %>%
  filter(cross > 2) %>%
  dplyr::select(Gene)

gene_levels_adc <- gene_levels$Gene

COMBI_plot2<-COMBI_plot2[COMBI_plot2$Gene%in%gene_levels$Gene,]
COMBI_plot3<-COMBI_plot3[COMBI_plot3$Gene%in%gene_levels$Gene,]
# COMBI_plot3
# COMBI_plot2 <- COMBI_plot2 %>%
#   mutate(Gene = factor(Gene, levels = gene_levels))
# 
# COMBI_plot3 <- COMBI_plot3 %>%
#   mutate(Gene = factor(Gene, levels = gene_levels))



# 2) Rebuild x1 and x2 so they use y = Gene (the factor) instead of reorder()
x1_clean <- ggplot(COMBI_plot3, aes(x = cross, y = reorder(COMBI_plot3$Gene,(-as.numeric(COMBI_plot3$cross))))) +
  geom_col(width = 1, fill = "black",col="white") +          # change fill if you want a different color
  #  facet_grid(. ~ New, switch = "x") +             # same faceting as x2 so columns align
  coord_flip() +
  scale_fill_manual(
    values = c("CNVs" = "#A6C8F4", "SNPs" = "#90EE90", "Multi" = "black", "INDELs" = "#F4A6A6"), 
    na.value = "white"
  ) +
  theme_minimal() +
  labs(y = "n Pair") +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(angle = 0, size = 12)# hide duplicate gene labels (keep them in x2 if desired)
  )

x1_clean
# Force factor levels 1–12 in both data frames
COMBI_plot4$New <- factor(COMBI_plot4$New, levels = 1:12)
COMBI_plot2$New <- factor(COMBI_plot2$New, levels = 1:12)


# Création de la grille complète
complete_data <- expand_grid(
  Cancer = unique(COMBI_plot2$Cancer),
  Gene = unique(COMBI_plot2$Gene),
  New = c("1","2","3","4","5","6","7","8","9","10","11","12")
)

#View(COMBI_plot2_filled)
# Fusion avec les vraies données
COMBI_plot2_filled <- complete_data %>%
  left_join(COMBI_plot2, by = c("Cancer", "Gene", "New")) %>%
  mutate(x_label = paste(New, Cancer, sep = "_"))

as.numeric(as.factor(COMBI_plot2_filled$Gene))

# Heatmap avec quadrillage continu (sans facette)

COMBI_plot2_filled <- COMBI_plot2_filled %>%
  mutate(
    Gene = factor(Gene, levels = unique(Gene[order(-cross,Gene)]))
  )

COMBI_plot2_filled

x2_clean <- ggplot(COMBI_plot2_filled, aes(
  x = factor(reorder(x_label,-as.numeric(New))),
  y = reorder(Gene, -cross),
  fill = FillColor
)) +
  coord_flip() +
  geom_tile(color = "grey", size = 0.3, na.rm = FALSE) +
  scale_fill_manual(
    values = c("CNVs" = "#A6C8F4", "SNPs" = "#90EE90", "Multi" = "black", "INDELs" = "#F4A6A6"),
    na.value = "white"
  ) +
  theme_minimal() +
  labs(x = "Pair", y = "Gene", fill = "Mutation Type") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,vjust = 0.95, size = 10),
    axis.text.y = element_blank(),
    panel.grid = element_blank()
  )


x2_clean

dev.off()






# 3) Prepare x0: remove its legend (we'll collect a single legend from x2_clean)
x0_clean <- x0 + theme(legend.position = "none")

# 4) Layout:
# top row = x0 | spacer
# bottom row = x2 | x1
top_row    <- plot_spacer()|x1_clean   # spacer prevents x1 from spanning up
top_row2<-top_row+plot_layout( widths = c(4, 17), guides = "collect") 
bottom_row <- x0_clean | x2_clean
bottom_row2<-bottom_row+plot_layout( widths = c(3, 15), guides = "collect") 

final_plot <- top_row2 / bottom_row2 +
  plot_layout(heights = c(1, 4), guides = "collect") &
  theme(legend.position = "none")

# show
final_plot

ggsave("CADD_ADC_total_mutations.pdf", plot = final_plot, width = 14, height = 6, limitsize = FALSE, device = 'pdf', dpi = 300)


#############################################
#############################################
#############################################



############################################################################################################
#####################################      SCC    ##########################################################
############################################################################################################
############################################################################################################
#####

COMBI_plot <- COMBI_SCC_analyze
COMBI_plot2<-unique(COMBI_plot[,c("New","Cancer","Gene","FillColor","cross")])
COMBI_plot3<-unique(COMBI_plot[,c("Gene","cross")])
COMBI_plot_selection<-unique(COMBI_plot[,c("Gene","cross")])

COMBI_plot4<-data.frame(New=c(1:12),VC=c(1:12))


COMBI_plot4$New <- factor(COMBI_plot4$New, levels = 1:12)
COMBI_SCC_analyze$New <- factor(COMBI_SCC_analyze$New, levels = 1:12)
COMBI_plot$Mut
library(dplyr)
library(tidyr)
library(ggplot2)

# On force New comme numérique de 1 à 12 et on complète
COMBI_plot <- COMBI_plot %>%
  mutate(New = as.numeric(New)) %>%
  group_by(New, Mut) %>%
  summarise(Mut_num = n(), .groups = "drop") %>%
  complete(
    New = 1:12,
    fill = list(Mut_num = 0)
  )

# Graphique
x0 <- ggplot(COMBI_plot, aes(x = reorder(factor(New),-New), y = Mut_num, fill = Mut)) +
  geom_col(col = "white") +
  coord_flip() +  # barres horizontales
  theme_minimal() +
  scale_x_discrete(position = "top", breaks = 1:12, labels = 1:12) +  
  scale_y_reverse() +  
  scale_fill_manual(
    values = c("CNVs" = "#A6C8F4", 
               "SNPs" = "#90EE90", 
               "Multi" = "black", 
               "INDELs" = "#F4A6A6"),
    na.value = "white"
  ) +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12, hjust = 0.5),
    axis.ticks.y = element_blank(),
    axis.text.x.bottom = element_text(size = 12),
    axis.ticks.x.bottom = element_line(),
    axis.text.x.top = element_text(size = 12),
    axis.ticks.x.top = element_line(),
    legend.position = "none"
  )

x0



# # 1) Build a single canonical gene ordering from COMBI_plot3 (adjust asc/desc below if you prefer)
# gene_levels <- COMBI_plot3 %>%
#   distinct(Gene, .keep_all = TRUE) %>%
#   arrange(cross) %>%    # ascending ORDER: smallest PROP -> bottom ; largest -> top
#   pull(Gene)%>%head(50)
COMBI_plot3$cross

library(dplyr)
COMBI_SCC_analyze

#COMBI_plot_selection2SCC<-COMBI_plot_selection[!COMBI_plot_selection$Gene%in%COMBI_SCC_analyze$Gene,]



gene_levels <- COMBI_SCC_analyze %>%
  filter(cross > 2) %>%
  dplyr::select(Gene)

gene_levels_SCC <- gene_levels$Gene

COMBI_plot2<-COMBI_plot2[COMBI_plot2$Gene%in%gene_levels$Gene,]
COMBI_plot3<-COMBI_plot3[COMBI_plot3$Gene%in%gene_levels$Gene,]
# COMBI_plot3
# COMBI_plot2 <- COMBI_plot2 %>%
#   mutate(Gene = factor(Gene, levels = gene_levels))
# 
# COMBI_plot3 <- COMBI_plot3 %>%
#   mutate(Gene = factor(Gene, levels = gene_levels))



# 2) Rebuild x1 and x2 so they use y = Gene (the factor) instead of reorder()
x1_clean <- ggplot(COMBI_plot3, aes(x = cross, y = reorder(COMBI_plot3$Gene,(-as.numeric(COMBI_plot3$cross))))) +
  geom_col(width = 1, fill = "black",col="white") +          # change fill if you want a different color
  #  facet_grid(. ~ New, switch = "x") +             # same faceting as x2 so columns align
  coord_flip() +
  scale_fill_manual(
    values = c("CNVs" = "#A6C8F4", "SNPs" = "#90EE90", "Multi" = "black", "INDELs" = "#F4A6A6"), 
    na.value = "white"
  ) +
  theme_minimal() +
  labs(y = "n Pair") +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(angle = 0, size = 12)# hide duplicate gene labels (keep them in x2 if desired)
  )

x1_clean
# Force factor levels 1–12 in both data frames
COMBI_plot4$New <- factor(COMBI_plot4$New, levels = 1:12)
COMBI_plot2$New <- factor(COMBI_plot2$New, levels = 1:12)


# Création de la grille complète
complete_data <- expand_grid(
  Cancer = unique(COMBI_plot2$Cancer),
  Gene = unique(COMBI_plot2$Gene),
  New = c("1","2","3","4","5","6","7","8","9","10","11","12")
)

#View(COMBI_plot2_filled)
# Fusion avec les vraies données
COMBI_plot2_filled <- complete_data %>%
  left_join(COMBI_plot2, by = c("Cancer", "Gene", "New")) %>%
  mutate(x_label = paste(New, Cancer, sep = "_"))

as.numeric(as.factor(COMBI_plot2_filled$Gene))

# Heatmap avec quadrillage continu (sans facette)

COMBI_plot2_filled <- COMBI_plot2_filled %>%
  mutate(
    Gene = factor(Gene, levels = unique(Gene[order(-cross,Gene)]))
  )

COMBI_plot2_filled

x2_clean <- ggplot(COMBI_plot2_filled, aes(
  x = factor(reorder(x_label,-as.numeric(New))),
  y = reorder(Gene, -cross),
  fill = FillColor
)) +
  coord_flip() +
  geom_tile(color = "grey", size = 0.3, na.rm = FALSE) +
  scale_fill_manual(
    values = c("CNVs" = "#A6C8F4", "SNPs" = "#90EE90", "Multi" = "black", "INDELs" = "#F4A6A6"),
    na.value = "white"
  ) +
  theme_minimal() +
  labs(x = "Pair", y = "Gene", fill = "Mutation Type") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,vjust = 0.95, size = 12),
    axis.text.y = element_blank(),
    panel.grid = element_blank()
  )


x2_clean

dev.off()






# 3) Prepare x0: remove its legend (we'll collect a single legend from x2_clean)
x0_clean <- x0 + theme(legend.position = "none")

# 4) Layout:
# top row = x0 | spacer
# bottom row = x2 | x1
top_row    <- plot_spacer()|x1_clean   # spacer prevents x1 from spanning up
top_row2<-top_row+plot_layout( widths = c(4, 17), guides = "collect") 
top_row2
bottom_row <- x0_clean | x2_clean
bottom_row2
bottom_row2<-bottom_row+plot_layout( widths = c(3, 15), guides = "collect") 

final_plot <- top_row2 / bottom_row2 +
  plot_layout(heights = c(1, 4), guides = "collect") &
  theme(legend.position = "none")

# show
final_plot

ggsave("CADD_SCC_total_mutations.pdf", plot = final_plot, width = 10, height = 6, limitsize = FALSE, device = 'pdf', dpi = 300)
ggsave("CADD_SCC_total_mutations2.pdf", plot = final_plot, width = 14, height = 6, limitsize = FALSE, device = 'pdf', dpi = 300)

COMBI_plot3_SCC<-COMBI_plot3






###########################
# ==================== 1. Conversion SYMBOL → ENTREZ ====================
message("Conversion des gènes SYMBOL en ENTREZID...")
gene_df <- bitr(COMBI_SCC_analyze$Gene,
                fromType = "SYMBOL",
                toType   = "ENTREZID",
                OrgDb    = org.Hs.eg.db)

# ==================== 2. Enrichissement KEGG ====================
message("Lancement de l'enrichissement KEGG...")
kegg_enrich <- enrichKEGG(
  gene         = gene_df$ENTREZID,
  organism     = 'hsa',
  pvalueCutoff = 0.05
)

# Rendre les résultats lisibles (ENTREZ → SYMBOL)
kegg_annotated <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_df <- as.data.frame(kegg_annotated)

# ==================== 3. Visualisation KEGG (barplot + dotplot) ====================
pdf("CADD_mutations_mutations_SCC_barplot_kegg_enrich.pdf", width = 10, height = 28)
x1_clean <- barplot(kegg_enrich, showCategory = 331)
dev.off()

pdf("CADD_mutations_mutations_SCC_dotplot_kegg_enrich.pdf", width = 14, height = 6)
x2_clean <- dotplot(kegg_enrich, showCategory = 331, orderBy = "GeneRatio", decreasing = FALSE) +
  coord_flip() +
  scale_x_reverse() +   
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


# ==================== 5. Jointure gènes enrichis/pathways ====================
message("Jointure avec les gènes CADD_COMBI...")

df_merged_SCC <- inner_join(COMBI_SCC, pathway_gene_df, by = c("Gene" = "GeneSymbol"))

# Comptage par pathway
COMBI_SCC_analyze <- df_merged_SCC %>%
  group_by(Pathway) %>%
  mutate(cross = n_distinct(New)) %>%
  ungroup() %>%
  arrange(desc(cross))

# Données pour le plot
COMBI_plot3 <- distinct(COMBI_SCC_analyze[, c("Pathway", "cross")])

# ==================== 6. Récupérer infos du KEGG enrichissement ====================
infos_generatios <- kegg_df[, c("Description", "GeneRatio", "p.adjust")]

COMBI_plot3_join <- inner_join(COMBI_plot3, infos_generatios, by = c("Pathway" = "Description")) %>%
  filter(p.adjust < 0.05) %>%
  mutate(Pathway = factor(Pathway, levels = Pathway[order(GeneRatio, decreasing = TRUE)]))

# ==================== 7. Plot final combiné ====================
x1_clean <- ggplot(COMBI_plot3_join, aes(x = cross, y = Pathway)) +
  geom_col(width = 1, fill = "black", col = "white") +
  coord_flip() +
  theme_minimal() +
  labs(y = "n Pair") +
  theme(
    axis.title   = element_blank(),
    panel.grid   = element_blank(),
    axis.text.x  = element_blank(),
    axis.text.y  = element_text(angle = 0, size = 12)
  )

bottom_row2 <- x1_clean / x2_clean

ggsave("CADD_mutations_SCC_total_mutations_barplot.pdf",
       plot = bottom_row2,
       width = 14, height = 6,
       limitsize = FALSE, device = 'pdf', dpi = 300)



#######################################################################################################

###########################
# ==================== 1. Conversion SYMBOL → ENTREZ ====================
message("Conversion des gènes SYMBOL en ENTREZID...")
gene_df <- bitr(COMBI_ADC_analyze$Gene,
                fromType = "SYMBOL",
                toType   = "ENTREZID",
                OrgDb    = org.Hs.eg.db)

# ==================== 2. Enrichissement KEGG ====================
message("Lancement de l'enrichissement KEGG...")
kegg_enrich <- enrichKEGG(
  gene         = gene_df$ENTREZID,
  organism     = 'hsa',
  pvalueCutoff = 0.05
)

# Rendre les résultats lisibles (ENTREZ → SYMBOL)
kegg_annotated <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_df <- as.data.frame(kegg_annotated)

# ==================== 3. Visualisation KEGG (barplot + dotplot) ====================
pdf("CADD_mutations_mutations_ADC_barplot_kegg_enrich.pdf", width = 10, height = 28)
x1_clean <- barplot(kegg_enrich, showCategory = 331)
dev.off()

pdf("CADD_mutations_mutations_ADC_dotplot_kegg_enrich.pdf", width = 14, height = 6)
x2_clean <- dotplot(kegg_enrich, showCategory = 331, orderBy = "GeneRatio", decreasing = FALSE) +
  coord_flip() +
  scale_x_reverse() +   
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



# ==================== 5. Jointure gènes enrichis/pathways ====================
message("Jointure avec les gènes CADD_COMBI...")

df_merged_ADC <- inner_join(COMBI_ADC, pathway_gene_df, by = c("Gene" = "GeneSymbol"))

# Comptage par pathway
COMBI_ADC_analyze <- df_merged_ADC %>%
  group_by(Pathway) %>%
  mutate(cross = n_distinct(New)) %>%
  ungroup() %>%
  arrange(desc(cross))

# Données pour le plot
COMBI_plot3 <- distinct(COMBI_ADC_analyze[, c("Pathway", "cross")])

# ==================== 6. Récupérer infos du KEGG enrichissement ====================
infos_generatios <- kegg_df[, c("Description", "GeneRatio", "p.adjust")]

COMBI_plot3_join <- inner_join(COMBI_plot3, infos_generatios, by = c("Pathway" = "Description")) %>%
  filter(p.adjust < 0.05) %>%
  mutate(Pathway = factor(Pathway, levels = Pathway[order(GeneRatio, decreasing = TRUE)]))

# ==================== 7. Plot final combiné ====================
x1_clean <- ggplot(COMBI_plot3_join, aes(x = cross, y = Pathway)) +
  geom_col(width = 1, fill = "black", col = "white") +
  coord_flip() +
  theme_minimal() +
  labs(y = "n Pair") +
  theme(
    axis.title   = element_blank(),
    panel.grid   = element_blank(),
    axis.text.x  = element_blank(),
    axis.text.y  = element_text(angle = 0, size = 12)
  )

bottom_row2 <- x1_clean / x2_clean

ggsave("CADD_mutations_ADC_total_genes_barplot.pdf",
       plot = bottom_row2,
       width = 14, height = 6,
       limitsize = FALSE, device = 'pdf', dpi = 300)

message("✅ Script terminé.")

bottom_row2
COMBI_plot3_join$genes<-gsub("(\\d+)\\/\\d+")
test_for_fulladc
#############################################################










######################################################  ######################################################  ######################################################  
######################################################  SELECTED Genes  ##################################################################
######################################################  ######################################################  ######################################################

####################################

COMBI_ADC<-COMBI%>%
  group_by(New,Gene,sig,Chromosome)%>%
  mutate(test=n_distinct(Cancer))%>%
  filter(test<2)%>%
  filter(Cancer=="ADC")

COMBI_SCC<-COMBI%>%
  group_by(New,Gene,sig,Chromosome)%>%
  mutate(test=n_distinct(Cancer))%>%
  filter(test<2)%>%
  filter(Cancer=="SCC")

COMBI_COM<-COMBI%>%
  group_by(New,Gene,sig,Chromosome)%>%
  mutate(test=n_distinct(Cancer))%>%
  filter(test==2)%>%
  filter(Cancer=="ADC")

####################################

COMBI_analyze<-COMBI%>%
  group_by(New,Gene,Cancer)%>%
  mutate(multi=n(),cross2=n_distinct(Mut))%>%
  ungroup()%>%
  group_by(Gene)%>%
  mutate(cross=n_distinct(New))%>%
  ungroup()%>%
  group_by(Gene,Mut)%>%
  mutate(cross3=n_distinct(New,Cancer))%>%
  ungroup()%>%
  mutate(FillColor = ifelse(cross2 > 1, "Multi", Mut))%>%
  arrange(cross)

COMBI_COM_analyze<-COMBI_COM%>%
  group_by(New,Gene,Cancer)%>%
  mutate(multi=n(),cross2=n_distinct(Mut))%>%
  ungroup()%>%
  group_by(Gene)%>%
  mutate(cross=n_distinct(New))%>%
  ungroup()%>%
  group_by(Gene,Mut)%>%
  mutate(cross3=n_distinct(New))%>%
  ungroup()%>%
  mutate(FillColor = ifelse(cross2 > 1, "Multi", Mut))%>%
  arrange(cross)

COMBI_ADC_analyze<-COMBI_ADC%>%
  group_by(New,Gene,Cancer)%>%
  mutate(multi=n(),cross2=n_distinct(Mut))%>%
  ungroup()%>%
  group_by(Gene)%>%
  mutate(cross=n_distinct(New))%>%
  ungroup()%>%
  group_by(Gene,Mut)%>%
  mutate(cross3=n_distinct(New))%>%
  ungroup()%>%
  mutate(FillColor = ifelse(cross2 > 1, "Multi", Mut))%>%
  arrange(cross)

COMBI_SCC_analyze<-COMBI_SCC%>%
  group_by(New,Gene,Cancer)%>%
  mutate(multi=n(),cross2=n_distinct(Mut))%>%
  ungroup()%>%
  group_by(Gene)%>%
  mutate(cross=n_distinct(New))%>%
  ungroup()%>%
  group_by(Gene,Mut)%>%
  mutate(cross3=n_distinct(New))%>%
  ungroup()%>%
  mutate(FillColor = ifelse(cross2 > 1, "Multi", Mut))%>%
  arrange(cross)

####################


COMBI_analyze$PROP<-(COMBI_analyze$cross/12)*100
COMBI_COM_analyze$PROP<-(COMBI_COM_analyze$cross/12)*100
COMBI_ADC_analyze$PROP<-(COMBI_ADC_analyze$cross/12)*100
COMBI_SCC_analyze$PROP<-(COMBI_SCC_analyze$cross/12)*100

#####

############################################################################################################
#####################################      ADC    ##########################################################
############################################################################################################
############################################################################################################
#####

COMBI_plot <- COMBI_ADC_analyze
COMBI_plot2<-unique(COMBI_plot[,c("New","Cancer","Gene","FillColor","cross")])
COMBI_plot3<-unique(COMBI_plot[,c("Gene","cross")])
COMBI_plot_selection<-unique(COMBI_plot[,c("Gene","cross")])

COMBI_plot4<-data.frame(New=c(12:1),VC=c(12:1))


COMBI_plot4$New <- factor(COMBI_plot4$New, levels = 12:1)
COMBI_ADC_analyze$New <- factor(COMBI_ADC_analyze$New, levels = 12:1)
COMBI_plot$New

x0 <- ggplot(COMBI_plot, aes(x = reorder(New, -New), fill = Mut)) +
  geom_bar(position = "stack", col = "white") +
  coord_flip() +  # Barres horizontales
  theme_minimal() +
  scale_x_discrete(position = "top") +  # 👈 Affiche les étiquettes de New à droite après flip
  scale_y_reverse() +  # Inverser l'axe des valeurs : 0 à droite
  scale_fill_manual(
    values = c("CNVs" = "#A6C8F4", "SNPs" = "#90EE90", "Multi" = "black", "INDELs" = "#F4A6A6"),
    na.value = "white"
  ) +
  #labs(x = "Number of mutation") +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    
    axis.text.y = element_text(angle = 0, size = 12,hjust = 0.5 ),     # Valeurs numériques (horizontal)
    axis.ticks.y = element_blank(),
    
    # Supprimer les étiquettes de New à gauche (en bas de x après flip)
    axis.text.x.bottom = element_text(angle = 0, size = 12),
    axis.ticks.x.bottom = element_line(),
    
    # Afficher les étiquettes de New à droite (en haut de x après flip)
    axis.text.x.top = element_text(angle = 0, size = 12),
    axis.ticks.x.top = element_line(),
    
    legend.position = "none"
  )

x0


# # 1) Build a single canonical gene ordering from COMBI_plot3 (adjust asc/desc below if you prefer)
# gene_levels <- COMBI_plot3 %>%
#   distinct(Gene, .keep_all = TRUE) %>%
#   arrange(cross) %>%    # ascending ORDER: smallest PROP -> bottom ; largest -> top
#   pull(Gene)%>%head(50)
COMBI_plot3$cross

library(dplyr)

#COMBI_plot_selection2ADC<-COMBI_plot_selection[!COMBI_plot_selection$Gene%in%COMBI_SCC_analyze$Gene,]
gene_levels <- COMBI_plot_selection %>%
  filter(cross > 2) %>%
  dplyr::select(Gene)

gene_levels_adc <- gene_levels$Gene

COMBI_plot2<-COMBI_plot2[COMBI_plot2$Gene%in%gene_levels$Gene,]
COMBI_plot3<-COMBI_plot3[COMBI_plot3$Gene%in%gene_levels$Gene,]
# COMBI_plot3
# COMBI_plot2 <- COMBI_plot2 %>%
#   mutate(Gene = factor(Gene, levels = gene_levels))
# 
# COMBI_plot3 <- COMBI_plot3 %>%
#   mutate(Gene = factor(Gene, levels = gene_levels))



# 2) Rebuild x1 and x2 so they use y = Gene (the factor) instead of reorder()
x1_clean <- ggplot(COMBI_plot3, aes(x = cross, y = reorder(COMBI_plot3$Gene,(-as.numeric(COMBI_plot3$cross))))) +
  geom_col(width = 1, fill = "black",col="white") +          # change fill if you want a different color
  #  facet_grid(. ~ New, switch = "x") +             # same faceting as x2 so columns align
  coord_flip() +
  scale_fill_manual(
    values = c("CNVs" = "#A6C8F4", "SNPs" = "#90EE90", "Multi" = "black", "INDELs" = "#F4A6A6"), 
    na.value = "white"
  ) +
  theme_minimal() +
  labs(y = "n Pair") +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(angle = 0, size = 12)# hide duplicate gene labels (keep them in x2 if desired)
  )

x1_clean
# Force factor levels 1–12 in both data frames
COMBI_plot4$New <- factor(COMBI_plot4$New, levels = 1:12)
COMBI_plot2$New <- factor(COMBI_plot2$New, levels = 1:12)


# Création de la grille complète
complete_data <- expand_grid(
  Cancer = unique(COMBI_plot2$Cancer),
  Gene = unique(COMBI_plot2$Gene),
  New = c("1","2","3","4","5","6","7","8","9","10","11","12")
)

#View(COMBI_plot2_filled)
# Fusion avec les vraies données
COMBI_plot2_filled <- complete_data %>%
  left_join(COMBI_plot2, by = c("Cancer", "Gene", "New")) %>%
  mutate(x_label = paste(New, Cancer, sep = "_"))

as.numeric(as.factor(COMBI_plot2_filled$Gene))

# Heatmap avec quadrillage continu (sans facette)

COMBI_plot2_filled <- COMBI_plot2_filled %>%
  mutate(
    Gene = factor(Gene, levels = unique(Gene[order(-cross,Gene)]))
  )

COMBI_plot2_filled

x2_clean <- ggplot(COMBI_plot2_filled, aes(
  x = factor(reorder(x_label,-as.numeric(New))),
  y = reorder(Gene, -cross),
  fill = FillColor
)) +
  coord_flip() +
  geom_tile(color = "grey", size = 0.3, na.rm = FALSE) +
  scale_fill_manual(
    values = c("CNVs" = "#A6C8F4", "SNPs" = "#90EE90", "Multi" = "black", "INDELs" = "#F4A6A6"),
    na.value = "white"
  ) +
  theme_minimal() +
  labs(x = "Pair", y = "Gene", fill = "Mutation Type") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,vjust = 0.95, size = 10),
    axis.text.y = element_blank(),
    panel.grid = element_blank()
  )


x2_clean

dev.off()






# 3) Prepare x0: remove its legend (we'll collect a single legend from x2_clean)
x0_clean <- x0 + theme(legend.position = "none")

# 4) Layout:
# top row = x0 | spacer
# bottom row = x2 | x1
top_row    <- plot_spacer()|x1_clean   # spacer prevents x1 from spanning up
top_row2<-top_row+plot_layout( widths = c(4, 17), guides = "collect") 
top_row2
bottom_row <- x0_clean | x2_clean
bottom_row2
bottom_row2<-bottom_row+plot_layout( widths = c(3, 15), guides = "collect") 

final_plot <- top_row2 / bottom_row2 +
  plot_layout(heights = c(1, 4), guides = "collect") &
  theme(legend.position = "none")

# show
final_plot

ggsave("CADD_GENES_ADC_total_genes.pdf", plot = final_plot, width = 14, height = 6, limitsize = FALSE, device = 'pdf', dpi = 300)


#############################################
#############################################
#############################################



############################################################################################################
#####################################      SCC    ##########################################################
############################################################################################################
############################################################################################################
#####

COMBI_plot <- COMBI_SCC_analyze
COMBI_plot2<-unique(COMBI_plot[,c("New","Cancer","Gene","FillColor","cross")])
COMBI_plot3<-unique(COMBI_plot[,c("Gene","cross")])
COMBI_plot_selection<-unique(COMBI_plot[,c("Gene","cross")])

COMBI_plot4<-data.frame(New=c(1:12),VC=c(1:12))


COMBI_plot4$New <- factor(COMBI_plot4$New, levels = 1:12)
COMBI_SCC_analyze$New <- factor(COMBI_SCC_analyze$New, levels = 1:12)
COMBI_plot$Mut
library(dplyr)
library(tidyr)
library(ggplot2)

# On force New comme numérique de 1 à 12 et on complète
COMBI_plot <- COMBI_plot %>%
  mutate(New = as.numeric(New)) %>%
  group_by(New, Mut) %>%
  summarise(Mut_num = n(), .groups = "drop") %>%
  complete(
    New = 1:12,
    fill = list(Mut_num = 0)
  )

# Graphique
x0 <- ggplot(COMBI_plot, aes(x = reorder(factor(New),-New), y = Mut_num, fill = Mut)) +
  geom_col(col = "white") +
  coord_flip() +  # barres horizontales
  theme_minimal() +
  scale_x_discrete(position = "top", breaks = 1:12, labels = 1:12) +  
  scale_y_reverse() +  
  scale_fill_manual(
    values = c("CNVs" = "#A6C8F4", 
               "SNPs" = "#90EE90", 
               "Multi" = "black", 
               "INDELs" = "#F4A6A6"),
    na.value = "white"
  ) +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12, hjust = 0.5),
    axis.ticks.y = element_blank(),
    axis.text.x.bottom = element_text(size = 12),
    axis.ticks.x.bottom = element_line(),
    axis.text.x.top = element_text(size = 12),
    axis.ticks.x.top = element_line(),
    legend.position = "none"
  )

x0



# # 1) Build a single canonical gene ordering from COMBI_plot3 (adjust asc/desc below if you prefer)
# gene_levels <- COMBI_plot3 %>%
#   distinct(Gene, .keep_all = TRUE) %>%
#   arrange(cross) %>%    # ascending ORDER: smallest PROP -> bottom ; largest -> top
#   pull(Gene)%>%head(50)
COMBI_plot3$cross

library(dplyr)
COMBI_SCC_analyze

#COMBI_plot_selection2SCC<-COMBI_plot_selection[!COMBI_plot_selection$Gene%in%COMBI_SCC_analyze$Gene,]



gene_levels <- COMBI_SCC_analyze %>%
  filter(cross > 2) %>%
  dplyr::select(Gene)

gene_levels_SCC <- gene_levels$Gene

COMBI_plot2<-COMBI_plot2[COMBI_plot2$Gene%in%gene_levels$Gene,]
COMBI_plot3<-COMBI_plot3[COMBI_plot3$Gene%in%gene_levels$Gene,]
# COMBI_plot3
# COMBI_plot2 <- COMBI_plot2 %>%
#   mutate(Gene = factor(Gene, levels = gene_levels))
# 
# COMBI_plot3 <- COMBI_plot3 %>%
#   mutate(Gene = factor(Gene, levels = gene_levels))



# 2) Rebuild x1 and x2 so they use y = Gene (the factor) instead of reorder()
x1_clean <- ggplot(COMBI_plot3, aes(x = cross, y = reorder(COMBI_plot3$Gene,(-as.numeric(COMBI_plot3$cross))))) +
  geom_col(width = 1, fill = "black",col="white") +          # change fill if you want a different color
  #  facet_grid(. ~ New, switch = "x") +             # same faceting as x2 so columns align
  coord_flip() +
  scale_fill_manual(
    values = c("CNVs" = "#A6C8F4", "SNPs" = "#90EE90", "Multi" = "black", "INDELs" = "#F4A6A6"), 
    na.value = "white"
  ) +
  theme_minimal() +
  labs(y = "n Pair") +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(angle = 0, size = 12)# hide duplicate gene labels (keep them in x2 if desired)
  )

x1_clean
# Force factor levels 1–12 in both data frames
COMBI_plot4$New <- factor(COMBI_plot4$New, levels = 1:12)
COMBI_plot2$New <- factor(COMBI_plot2$New, levels = 1:12)


# Création de la grille complète
complete_data <- expand_grid(
  Cancer = unique(COMBI_plot2$Cancer),
  Gene = unique(COMBI_plot2$Gene),
  New = c("1","2","3","4","5","6","7","8","9","10","11","12")
)

#View(COMBI_plot2_filled)
# Fusion avec les vraies données
COMBI_plot2_filled <- complete_data %>%
  left_join(COMBI_plot2, by = c("Cancer", "Gene", "New")) %>%
  mutate(x_label = paste(New, Cancer, sep = "_"))

as.numeric(as.factor(COMBI_plot2_filled$Gene))

# Heatmap avec quadrillage continu (sans facette)

COMBI_plot2_filled <- COMBI_plot2_filled %>%
  mutate(
    Gene = factor(Gene, levels = unique(Gene[order(-cross,Gene)]))
  )

COMBI_plot2_filled

x2_clean <- ggplot(COMBI_plot2_filled, aes(
  x = factor(reorder(x_label,-as.numeric(New))),
  y = reorder(Gene, -cross),
  fill = FillColor
)) +
  coord_flip() +
  geom_tile(color = "grey", size = 0.3, na.rm = FALSE) +
  scale_fill_manual(
    values = c("CNVs" = "#A6C8F4", "SNPs" = "#90EE90", "Multi" = "black", "INDELs" = "#F4A6A6"),
    na.value = "white"
  ) +
  theme_minimal() +
  labs(x = "Pair", y = "Gene", fill = "Mutation Type") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,vjust = 0.95, size = 12),
    axis.text.y = element_blank(),
    panel.grid = element_blank()
  )


x2_clean

dev.off()






# 3) Prepare x0: remove its legend (we'll collect a single legend from x2_clean)
x0_clean <- x0 + theme(legend.position = "none")

# 4) Layout:
# top row = x0 | spacer
# bottom row = x2 | x1
top_row    <- plot_spacer()|x1_clean   # spacer prevents x1 from spanning up
top_row2<-top_row+plot_layout( widths = c(4, 17), guides = "collect") 
top_row2
bottom_row <- x0_clean | x2_clean
bottom_row2
bottom_row2<-bottom_row+plot_layout( widths = c(3, 15), guides = "collect") 

final_plot <- top_row2 / bottom_row2 +
  plot_layout(heights = c(1, 4), guides = "collect") &
  theme(legend.position = "none")

# show
final_plot

ggsave("CADD_GENES_SCC_total_genes.pdf", plot = final_plot, width = 10, height = 6, limitsize = FALSE, device = 'pdf', dpi = 300)
ggsave("CADD_GENES_SCC_total_genes2.pdf", plot = final_plot, width = 14, height = 6, limitsize = FALSE, device = 'pdf', dpi = 300)

COMBI_plot3_SCC<-COMBI_plot3







###########################
# ==================== 1. Conversion SYMBOL → ENTREZ ====================
message("Conversion des gènes SYMBOL en ENTREZID...")
gene_df <- bitr(COMBI_SCC_analyze$Gene,
                fromType = "SYMBOL",
                toType   = "ENTREZID",
                OrgDb    = org.Hs.eg.db)

# ==================== 2. Enrichissement KEGG ====================
message("Lancement de l'enrichissement KEGG...")
kegg_enrich <- enrichKEGG(
  gene         = gene_df$ENTREZID,
  organism     = 'hsa',
  pvalueCutoff = 0.05
)

# Rendre les résultats lisibles (ENTREZ → SYMBOL)
kegg_annotated <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_df <- as.data.frame(kegg_annotated)

# ==================== 3. Visualisation KEGG (barplot + dotplot) ====================
pdf("CADD_GENES_mutations_SCC_barplot_kegg_enrich.pdf", width = 10, height = 28)
x1_clean <- barplot(kegg_enrich, showCategory = 331)
dev.off()

pdf("CADD_GENES_mutations_SCC_dotplot_kegg_enrich.pdf", width = 14, height = 6)
x2_clean <- dotplot(kegg_enrich, showCategory = 331, orderBy = "GeneRatio", decreasing = FALSE) +
  coord_flip() +
  scale_x_reverse() +   
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


# ==================== 5. Jointure gènes enrichis/pathways ====================
message("Jointure avec les gènes CADD_COMBI...")

df_merged_SCC <- inner_join(COMBI_SCC, pathway_gene_df, by = c("Gene" = "GeneSymbol"))

# Comptage par pathway
COMBI_SCC_analyze <- df_merged_SCC %>%
  group_by(Pathway) %>%
  mutate(cross = n_distinct(New)) %>%
  ungroup() %>%
  arrange(desc(cross))

# Données pour le plot
COMBI_plot3 <- distinct(COMBI_SCC_analyze[, c("Pathway", "cross")])

# ==================== 6. Récupérer infos du KEGG enrichissement ====================
infos_generatios <- kegg_df[, c("Description", "GeneRatio", "p.adjust")]

COMBI_plot3_join <- inner_join(COMBI_plot3, infos_generatios, by = c("Pathway" = "Description")) %>%
  filter(p.adjust < 0.05) %>%
  mutate(Pathway = factor(Pathway, levels = Pathway[order(GeneRatio, decreasing = TRUE)]))

# ==================== 7. Plot final combiné ====================
x1_clean <- ggplot(COMBI_plot3_join, aes(x = cross, y = Pathway)) +
  geom_col(width = 1, fill = "black", col = "white") +
  coord_flip() +
  theme_minimal() +
  labs(y = "n Pair") +
  theme(
    axis.title   = element_blank(),
    panel.grid   = element_blank(),
    axis.text.x  = element_blank(),
    axis.text.y  = element_text(angle = 0, size = 12)
  )

bottom_row2 <- x1_clean / x2_clean

ggsave("CADD_GENES_SCC_total_genes_barplot.pdf",
       plot = bottom_row2,
       width = 14, height = 6,
       limitsize = FALSE, device = 'pdf', dpi = 300)



#######################################################################################################

###########################
# ==================== 1. Conversion SYMBOL → ENTREZ ====================
message("Conversion des gènes SYMBOL en ENTREZID...")
gene_df <- bitr(COMBI_ADC_analyze$Gene,
                fromType = "SYMBOL",
                toType   = "ENTREZID",
                OrgDb    = org.Hs.eg.db)

# ==================== 2. Enrichissement KEGG ====================
message("Lancement de l'enrichissement KEGG...")
kegg_enrich <- enrichKEGG(
  gene         = gene_df$ENTREZID,
  organism     = 'hsa',
  pvalueCutoff = 0.05
)

# Rendre les résultats lisibles (ENTREZ → SYMBOL)
kegg_annotated <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_df <- as.data.frame(kegg_annotated)

# ==================== 3. Visualisation KEGG (barplot + dotplot) ====================
pdf("CADD_GENES_mutations_ADC_barplot_kegg_enrich.pdf", width = 10, height = 28)
x1_clean <- barplot(kegg_enrich, showCategory = 331)
dev.off()

pdf("CADD_GENES_mutations_ADC_dotplot_kegg_enrich.pdf", width = 14, height = 6)
x2_clean <- dotplot(kegg_enrich, showCategory = 331, orderBy = "GeneRatio", decreasing = FALSE) +
  coord_flip() +
  scale_x_reverse() +   
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



# ==================== 5. Jointure gènes enrichis/pathways ====================
message("Jointure avec les gènes CADD_COMBI...")

df_merged_ADC <- inner_join(COMBI_ADC, pathway_gene_df, by = c("Gene" = "GeneSymbol"))

# Comptage par pathway
COMBI_ADC_analyze <- df_merged_ADC %>%
  group_by(Pathway) %>%
  mutate(cross = n_distinct(New)) %>%
  ungroup() %>%
  arrange(desc(cross))

# Données pour le plot
COMBI_plot3 <- distinct(COMBI_ADC_analyze[, c("Pathway", "cross")])

# ==================== 6. Récupérer infos du KEGG enrichissement ====================
infos_generatios <- kegg_df[, c("Description", "GeneRatio", "p.adjust")]

COMBI_plot3_join <- inner_join(COMBI_plot3, infos_generatios, by = c("Pathway" = "Description")) %>%
  filter(p.adjust < 0.05) %>%
  mutate(Pathway = factor(Pathway, levels = Pathway[order(GeneRatio, decreasing = TRUE)]))

# ==================== 7. Plot final combiné ====================
x1_clean <- ggplot(COMBI_plot3_join, aes(x = cross, y = Pathway)) +
  geom_col(width = 1, fill = "black", col = "white") +
  coord_flip() +
  theme_minimal() +
  labs(y = "n Pair") +
  theme(
    axis.title   = element_blank(),
    panel.grid   = element_blank(),
    axis.text.x  = element_blank(),
    axis.text.y  = element_text(angle = 0, size = 12)
  )

bottom_row2 <- x1_clean / x2_clean

ggsave("CADD_GENES_ADC_total_genes_barplot.pdf",
       plot = bottom_row2,
       width = 14, height = 6,
       limitsize = FALSE, device = 'pdf', dpi = 300)

message("✅ Script terminé.")



































Loas_files3<-read.delim(file = "Mutect2_VCF0.8.maf", header = T, sep = "\t")

# REVEL_score > 0.5 &
#  MutationTaster_pred == "D" &
# Variant_Classification %in% c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site"))

Loas_files_INDELs_mutect2<-Loas_files3[Loas_files3$Variant_Type=="INS"|Loas_files3$Variant_Type=="DEL",]
Loas_files_SNPs<-Loas_files3[Loas_files3$Variant_Type=="SNP"|Loas_files3$Variant_Type=="TNP"|Loas_files3$Variant_Type=="DNP",]



datas_cliniques<-read.delim(file = "./Datas_cliniques.txt", header = T, sep = "\t")
CODING_snps<-Loas_files_SNPs[Loas_files_SNPs$Variant_Classification=="Missense_Mutation"|Loas_files_SNPs$Variant_Classification=="Nonsense_Mutation"|Loas_files_SNPs$Variant_Classification=="Nonstop_Mutation"|Loas_files_SNPs$Variant_Classification=="Translation_Start_Site",]
CODING_INDELs<-Loas_files_INDELs_mutect2[!Loas_files_INDELs_mutect2$Variant_Classification=="Non exonic",]

combined_muts<-rbind(CODING_snps,CODING_INDELs)

###############################################################################
##########################################  InDels ######################################################

Result_SNPs<- data.frame(Paire = integer(),
                         ADC = integer(),
                         SCC = integer(),
                         class = numeric(),
                         New = numeric(),
                         Subtype = integer(),
                         Integration = integer(),
                         Cancer = integer(),
                         Variant_Classification = integer(),
                         prop_commun = numeric(),
                         count_variant = numeric(),
                         statut = integer(),
                         compare = integer(),
                         stringsAsFactors = FALSE)




for (i in 1:length(unique(combined_muts$Paire))){
  my_sample1<-combined_muts[which(combined_muts$Paire==unique(combined_muts$Paire)[i]&combined_muts$Cancer=="SCC"),]
  my_sample3<-combined_muts[which(combined_muts$Paire==unique(combined_muts$Paire)[i]&combined_muts$Cancer=="ADC"),]
  ################
  for (y in 1:length(unique(combined_muts$Paire))){
    combined_muts$Reference_Allele
    my_sample2<-combined_muts[which(combined_muts$Paire==unique(combined_muts$Paire)[y]&combined_muts$Cancer=="ADC"),]
    my_sample<-rbind(my_sample1,my_sample2)

    # combined_muts$Variant_Classification
    SUM_my_sample<-my_sample%>%
      mutate(total_mut=n())%>%
      group_by(Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2,End_Position) %>%
      mutate(iter= n_distinct(Cancer))%>%
      ungroup()%>%
      mutate(n_com=sum(iter == 2),prop_commun=((sum(iter == 2)/2)/(total_mut-(sum(iter == 2)/2)))*100)%>%
      group_by(Paire,Cancer, Variant_Classification,Variant_Type,prop_commun,total_mut,n_com) %>%
      summarise(count_variant = n())
    
    SUM_my_sample_n<- merge(x=datas_cliniques, y= SUM_my_sample, by = "Paire" )
    
    if(nrow(SUM_my_sample_n)>0){
      
      SUM_my_sample_n$statut<-"CTR"
      SUM_my_sample_n$compare<-paste(i,y,sep = "_")
      
      if(i==y){
        
        SUM_my_sample_n$statut<-as.character(i)
        
      }  
      
      Result_SNPs<-rbind(Result_SNPs,SUM_my_sample_n)
      
    }
  }
}


####
Result3 <- Result_SNPs[Result_SNPs$statut != "CTR",]

  


#########################################################################################################
##########################################  fill tmb ######################################################


Result4<-Result3%>%
  group_by(New,Cancer)%>%
  summarise(counting=sum(count_variant))


Result4$counting<-Result4$counting/37






Result4 <- Result4 %>%
  ungroup() %>%                           # drop grouping
  complete(
    New    = 1:12,
    Cancer = c("ADC", "SCC"),
    fill   = list(counting = 0)          # set missing counting values to 0
  )



Result4<- merge(x=datas_cliniques, y= Result4, by = "New" )

Result4$Integration[Result4$Integration=="INTE_VS"]<-"INTE"
Result4$Integration[Result4$Integration=="VS"]<-"EPI"
TMB_plot1 <- ggplot(Result4, aes(x = Cancer, y = counting)) +
  geom_violin(alpha = 1, fill = "white",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  geom_violin(alpha = 0.4, aes(fill = Cancer),col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +  # Violin plot for distributione
  #geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(shape = 16, size = 3, width = 0.15, alpha = 0.7) +
  scale_y_log10() +
  #scale_fill_brewer(palette = "Set2") +
  
  scale_fill_manual(
    values = c(
      "ADC" = "#FF6666",
      "SCC" = "#3399ff",
      "Commun" = "gray"
    )
  )+
  
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  labs(y = "Number of mutation", title = "Integration Statut") +
  
  # 🟡 Add statistical test brackets
  stat_compare_means(method = "t.test", label.y = 1,label="p.format",label.x.npc = "right") 

TMB_plot1

Result4$type[Result4$type=="HPV58"]<-"Other HPVs"
Result4$type[Result4$type=="HPV45"]<-"Other HPVs"


TMB_plot2 <- ggplot(Result4, aes(x = type, y = counting)) +
  geom_violin(alpha = 1, fill = "white",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  geom_violin(alpha = 0.4, fill = "gray",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  #geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(aes(col=Cancer),shape = 16, size = 3, width = 0.15, alpha = 0.7) +  # Add points for individual values
  #scale_y_log10() +  # Use log scale for better visualization
  scale_color_manual(
    values = c(
      "ADC" = "#FF6666",
      "SCC" = "#3399ff",
      "Commun" = "grey"
    )
  )+
  scale_y_log10() +  # Use log scale for better visualization
  scale_fill_brewer(palette = "Set2") +  # Nice color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text( size = 14,angle = 45,hjust = 1,vjust = 1),
    axis.text.y =  element_text( size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  
  labs(y = "Number of mutation", title = "HPV subtypes")+
  stat_compare_means(method = "anova", label.y = 1,label="p.format",label.x.npc = "right") 


TMB_plot3 <- ggplot(Result4, aes(x = FIGO, y = counting)) +
  geom_violin(alpha = 1, fill = "white",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  geom_violin(alpha = 0.4, fill = "gray",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  #geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(aes(col=Cancer),shape = 16, size = 3, width = 0.15, alpha = 0.7) +  # Add points for individual values
  #scale_y_log10() +  # Use log scale for better visualization
  scale_color_manual(
    values = c(
      "ADC" = "#FF6666",
      "SCC" = "#3399ff",
      "Commun" = "gray"
    )
  )+
  scale_y_log10() +  # Use log scale for better visualization
  scale_fill_brewer(palette = "Set2") +  # Nice color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text( size = 14,angle = 45,hjust = 1,vjust = 1),
    axis.text.y =  element_text( size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  labs(y = "Number of mutation", title = "FIGO cancer classification")+
  stat_compare_means(method = "anova", label.y = 1,label="p.format",label.x.npc = "right") 
?stat_compare_means


TMB_plot4 <- ggplot(Result4, aes(x = Integration, y = counting)) +
  geom_violin(alpha = 1, fill = "white",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  geom_violin(alpha = 0.4, fill = "gray",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  #geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(aes(col=Cancer),shape = 16, size = 3, width = 0.15, alpha = 0.7) +  # Add points for individual values
  #scale_y_log10() +  # Use log scale for better visualization
  scale_color_manual(
    values = c(
      "ADC" = "#FF6666",
      "SCC" = "#3399ff",
      "Commun" = "gray"
    )
  )+
  scale_y_log10() +  # Use log scale for better visualization
  scale_fill_brewer(palette = "Set2") +  # Nice color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text( size = 14,angle = 45,hjust = 1,vjust = 1),
    axis.text.y =  element_text( size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  labs(y = "Number of mutation", title = "Integration Statut")+
  stat_compare_means(method = "t.test", label.y = 1,label="p.format",label.x.npc = "right") 



combi_row <- TMB_plot1 + TMB_plot2 + TMB_plot3 + TMB_plot4+
  plot_layout(ncol = 4)

combi_row

# Save the plot
ggsave("17_TMB_plot_by_Strain2_FULLs_select.pdf", width = 12, height = 6, dpi = 300)














