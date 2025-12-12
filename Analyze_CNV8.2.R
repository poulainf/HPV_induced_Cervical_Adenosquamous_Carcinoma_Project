#!/usr/bin/Rscript

library(car)
library(dplyr)
library(Formula)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(ggrepel)
library(ggh4x)
library(ggsignif)
library(grid)
library(Hmisc)
library(lattice)
library(lsa)
library(magrittr)
library(patchwork)
library(plotly)
library(randomcoloR)
library(RColorBrewer)
library(Rcmdr)
library(reshape2)
library(scales)
library(stringr)
library(survival)
library(tidyr)
library(wesanderson)
library(patchwork)
library(maftools)
library(gridExtra)
library(CINmetrics)


Loas_files_POINT<-read.delim(file = "CNVs/Point_CR_COMBINED.txt")

Loas_files_POINT$combi<-paste0(Loas_files_POINT$Paire,"_",Loas_files_POINT$Type)
Loas_files_POINT$Type[Loas_files_POINT$Type=="CIN3"]<-"SCC"
Loas_files_POINT$Type[Loas_files_POINT$Type=="ADC2"]<-"ADC"  
Loas_files_POINT$Type[Loas_files_POINT$Type=="SCC2"]<-"SCC"
Loas_files_POINT$Type[Loas_files_POINT$Type=="SCC3"]<-"SCC" 

Loas_files_POINT$Type_order<-1
Loas_files_POINT$Type_order[which(Loas_files_POINT$Type=="ADC")]<-2
Loas_files_POINT$Type_order[which(Loas_files_POINT$Type=="SCC")]<-3


Loas_files_POINT<-unique(Loas_files_POINT)
My_type=unique(Loas_files_POINT$Type)

Loas_files_POINT$ordre<-gsub('chr(.+)','\\1',Loas_files_POINT$CONTIG) 
Loas_files_POINT$ordre[which(Loas_files_POINT$CONTIG=="X")]<-22
Loas_files_POINT$LOC<-Loas_files_POINT$START+((Loas_files_POINT$END-Loas_files_POINT$START)/2)

Loas_files_POINT$Paire_num<-gsub("Paire_(\\d+)","\\1",Loas_files_POINT$Paire)
Loas_files_POINT<-unique(Loas_files_POINT)

#Loas_files_POINT <- Loas_files_POINT[!(Loas_files_POINT$Paire %in% c(2,5, 7, 9,10,15)), ]


Loas_files<-read.delim(file = "/media/florian2/T7/To_send/CNVs/Test_COMBINED.txt")
Loas_files$Type_order<-1
Loas_files$Type[Loas_files$Type=="CIN3"]<-"SCC"
Loas_files$Type[Loas_files$Type=="ADC2"]<-"ADC"  
Loas_files$Type[Loas_files$Type=="SCC2"]<-"SCC"
Loas_files$Type[Loas_files$Type=="SCC3"]<-"SCC" 

Loas_files$Type_order[which(Loas_files$Type=="ADC")]<-2
Loas_files$Type_order[which(Loas_files$Type=="SCC")]<-3
Loas_files<-unique(Loas_files)



Loas_files$LOC<-Loas_files$START+((Loas_files$END-Loas_files$START)/2)


Loas_files$ordre<-gsub('chr(.+)','\\1',Loas_files$CONTIG)
Loas_files$ordre[which(Loas_files$CONTIG=="X")]<-22






split_segments <- function(data, breakpoints) {
  result <- data.frame()
  for (i in 1:(length(breakpoints) - 1)) {
    start <- breakpoints[i]
    end <- breakpoints[i + 1]
    for (row in 1:nrow(data)) {
      if (data$PosSTART[row] <= start & data$PosSTOP[row] >= end) {
        new_row <- data[row, ]
        new_row$PosSTART <- start
        new_row$PosSTOP <- end
        result <- rbind(result, new_row)
      }
    }
  }
  return(result)
}



my_full<-matrix(ncol = 6)



Loas_files$Type

Loas_files$Paire_num<-gsub("Paire_(\\d+)","\\1",Loas_files$Paire)
# Loas_files <- merge(Loas_files, Loas_samples, by.x = "Paire_num", by.y = "Old")
Loas_files$Paire<-as.numeric(Loas_files$Paire_num)


Loas_files <- Loas_files[!(Loas_files$Paire %in% c(2,5, 7, 9,10,15)), ]
Loas_files$Paire.y<-NULL
Loas_files<-unique(Loas_files)


Loas_files$size<-Loas_files$END-Loas_files$START
Loas_files<-Loas_files[!Loas_files$NUM_POINTS_COPY_RATIO<5,]
#for (hhg in c(0.3,0.4,0.5)) {
hhg<-0.3


 
 Resultes_ctr<- data.frame(Type = integer(),Total_CNV = numeric(),sign = integer(),share = numeric(),commun = numeric(),Paire = integer(), stringsAsFactors = FALSE)
 Resultes<- data.frame(Type = integer(),Total_CNV = numeric(),sign = integer(),share = numeric(),commun = numeric(),Paire = integer(), stringsAsFactors = FALSE)
 My_paire=unique(Loas_files$Paire)
 full_rect<-data.frame()
 Loas_files_POINT$Type
 Loas_files_POINT$Paire_num<-as.numeric(Loas_files_POINT$Paire_num)
 
 Loas_files$Type
 My_paire
 for ( i in 1:length(My_paire)){
   
   for ( y in 1:length(My_paire)){
     
     #i<-8
     #y<-8
     My_paire[i]
     My_paire[y]

     
     Loas_files_select_i<-Loas_files[which(Loas_files$Paire==My_paire[i] & Loas_files$Type=="ADC"),]
     my_ext1_i<-Loas_files_POINT[which(Loas_files_POINT$Paire_num==My_paire[i]& Loas_files_POINT$Type=="ADC"),]
     
     Loas_files_select_y<-Loas_files[which(Loas_files$Paire==My_paire[y]& Loas_files$Type=="SCC"),]
     my_ext1_y<-Loas_files_POINT[which(Loas_files_POINT$Paire_num==My_paire[y]& Loas_files_POINT$Type=="SCC"),]
     
     Loas_files_select<-rbind(Loas_files_select_i,Loas_files_select_y)
     my_ext1<-rbind(my_ext1_i,my_ext1_y)  
     
     my_ext1$chr_pos<-gsub("chr(.+)","\\1",my_ext1$CONTIG)
     my_ext1$chr_pos[which(my_ext1$chr_pos=="X")]<-23
     Loas_files_select$chr_pos<-gsub("chr(.+)","\\1",Loas_files_select$CONTIG)
     Loas_files_select$chr_pos[which(Loas_files_select$chr_pos=="X")]<-23
     
     
     my_ext1$Type[which(my_ext1$Type=="CIN3")]<-"SCC"
     my_ext1$Type[which(my_ext1$Type=="SCC2")]<-"SCC"
     my_ext1$Type[which(my_ext1$Type=="ADC2")]<-"ADC"
     my_ext1$Type[which(my_ext1$Type=="SCC3")]<-"SCC"
     
     Loas_files_select$Type[which(Loas_files_select$Type=="SCC3")]<-"SCC"
     Loas_files_select$Type[which(Loas_files_select$Type=="CIN3")]<-"SCC"
     Loas_files_select$Type[which(Loas_files_select$Type=="SCC2")]<-"SCC"
     Loas_files_select$Type[which(Loas_files_select$Type=="ADC2")]<-"ADC"
     log2(0.8)
     log2(1.2)
     

     ###################
     ##################
     
     LEN<-0
     
     
     # Step 1: Calculate the length of each contig
     TEST3 <- Loas_files_select %>%
       group_by(CONTIG) %>%
       mutate(cont_len = max(END)) %>%
       ungroup()
     
     # Step 2: Create a dataframe with unique contigs and their cumulative lengths
     LEN <- TEST3 %>%
       select(CONTIG, cont_len) %>%
       distinct() %>%
       arrange(CONTIG) %>%
       mutate(CUM = cumsum(as.numeric(lag(cont_len, default = 0))))
     
     # Step 3: Merge the LEN dataframe with the original data to get cumulative lengths
     Loas_files_select <- merge(LEN, TEST3, by = "CONTIG")
     
     # Step 4: Calculate the linear positions
     Loas_files_select <- Loas_files_select %>%
       mutate(PosSTART = START + CUM,
              PosSTOP = END + CUM)
     
     
     ########################
     ########################
     
     Loas_files_select$roundSTART<-round(Loas_files_select$PosSTART/1000000,digits = 1)
     Loas_files_select$roundEND<-round(Loas_files_select$PosSTOP/1000000,digits = 1)
     Loas_files_select$roundSTART
     Loas_files_select$roundEND  
     
     Loas_files_select2<-unique(Loas_files_select)
     my_ext1_select23<-data.frame()
     if (nrow(Loas_files_select2)>0){
       breakpoints <- sort(unique(c(Loas_files_select2$PosSTART, Loas_files_select2$PosSTOP)))
       split_data <- split_segments(Loas_files_select2, breakpoints)
       split_data$sign<-"Amp"
       
       split_data$sign[which(split_data$MEAN_LOG2_COPY_RATIO<0)]<-"Del"
       
       unique(split_data$CALL)
       # abs(abs(mean(Loas_files_select_i$LOG2_COPY_RATIO_POSTERIOR_10))-abs(mean(Loas_files_select_i$LOG2_COPY_RATIO_POSTERIOR_90)))*2
       # abs(abs(mean(Loas_files_select_y$LOG2_COPY_RATIO_POSTERIOR_10))-abs(mean(Loas_files_select_y$LOG2_COPY_RATIO_POSTERIOR_90)))*2
       
       
       split_data_select<-unique(split_data[which(split_data$MEAN_LOG2_COPY_RATIO< -hhg| split_data$MEAN_LOG2_COPY_RATIO>hhg ),])
       
       #split_data_select<-unique(split_data[which(split_data$CALL!="0 "),])
       split_data_select
       
       split_data_select$PosSIZE<-split_data_select$PosSTOP-split_data_select$PosSTART
       split_data_select$PROP_SIZE<-split_data_select$PosSIZE/split_data_select$size
       
       #split_data_select<-split_data_select[split_data_select$PROP_SIZE>0.5,]
  
 
       my_ext1_select23 <- split_data_select %>%
         # per-Type: unique (START, END) count
         group_by(Type) %>%
         mutate(Total_CNV = n_distinct(START, END)) %>%
         ungroup() %>%
         
         # global: sum one Total_CNV per Type (no undercount when counts are equal)
         mutate(Full_CNV = sum(Total_CNV[!duplicated(Type)])) %>%
         
         # per-region: how many Types have this (PosSTART, PosSTOP)
         group_by(PosSTART, PosSTOP) %>%
         mutate(ncount = n_distinct(Type)) %>%
         ungroup() %>%
         
         # flag shared regions (present in >= 2 Types)
         mutate(shared_flag = as.integer(ncount >= 2)) %>%
         
         # per-Type: number of shared regions
         #group_by(Type) %>%
         mutate(share = sum(shared_flag)) %>%
         ungroup() %>%
         
         # % shared relative to non-shared count
         mutate(commun = ((share/2) / (Full_CNV - (share/2)) * 100)) %>%
         
         # keep one row per grouping key (matches your original intent)
         group_by(Type, Total_CNV, sign, share, commun, Paire) %>%
         summarise(.groups = "drop")
       
       
       my_ext1_select23$commun
       my_ext1_select23$share
       
     }
#     
     if (i == y){
#       
       full_rect<-rbind(full_rect,split_data_select)
       split_data_select_up<-unique(split_data[which(split_data$MEAN_LOG2_COPY_RATIO< -hhg| split_data$MEAN_LOG2_COPY_RATIO>hhg ),])
       split_data_select_down<-unique(split_data[which(split_data$MEAN_LOG2_COPY_RATIO>= -hhg& split_data$MEAN_LOG2_COPY_RATIO<=hhg ),])
#       
       split_data_select_down$sign<-"Neutral"
       split_data_select_fig<-rbind(split_data_select_down,split_data_select_up)
       split_data_select_fig$NUM_POINTS_COPY_RATIO<-NULL
       split_data_select_fig$MEAN_LOG2_COPY_RATIO<-NULL
       split_data_select_fig<-unique(split_data_select_fig)
#       
       split_data_select_fig$Type
       # 2. Separate into two data frames
       merged_df <- split_data_select_fig %>%
         filter(sign != "Neutral") %>%
         group_by(CONTIG, PosSTART, PosSTOP) %>%
         mutate(compt = n_distinct(Type)) %>%
         ungroup() %>%
         filter(compt == 2) %>%
         mutate(Type = "Commune") %>%
         distinct()  # to remove duplicates
       
       merged_df$compt<-NULL
       split_data_select_fig<-rbind(split_data_select_fig,merged_df)
       
       
       
       # Build background rectangle data by CONTIG and Type
     background_rects <- split_data %>%
         group_by(Type, CONTIG,chr_pos) %>%
         summarise(
           xmin = min(PosSTART, na.rm = TRUE),
           xmax = max(PosSTOP, na.rm = TRUE),
           .groups = "drop"
         ) %>%
         mutate(
          ymin = 0,
          ymax = 1,
          fill = "background"
        )

      background_rects$Type<-"ADC"
      background_rects<-unique(background_rects)
      background_rects2<-background_rects
      background_rects2$Type<-"SCC"
      background_rects3<-background_rects
      background_rects3$Type<-"Commune"
      background_rects<-rbind(background_rects,background_rects2,background_rects3)

      split_data_select_fig$ordre_split<-1
      split_data_select_fig$ordre_split[split_data_select_fig$Type=="SCC"]<-2
      split_data_select_fig$ordre_split[split_data_select_fig$Type=="Commune"]<-3

      background_rects$ordre_split<-1
      background_rects$ordre_split[background_rects$Type=="SCC"]<-2
      background_rects$ordre_split[background_rects$Type=="Commune"]<-3





      split_data_select_fig<-split_data_select_fig[!split_data_select_fig$sign=="Neutral",]

       x2<-ggplot() +
      
         # 1. Grey background rectangles per CONTIG
         geom_rect(
           data = background_rects,
           aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
           fill = "gray95",  # or another grey tone
           color = NA
         ) +
      
         # 2. CNV-specific rectangles (Amp/Del/Neutral)
         geom_rect(
           data = split_data_select_fig,
           aes(xmin = PosSTART, xmax = PosSTOP, ymin = 0, ymax = 1, fill = sign),
           color = NA
         ) +
      
         # 3. Manual color palette
         scale_fill_manual(
           values = c(
             "Amp" = "#9ACBB3",
             "Del" = "#C89CCB",
             "Neutral" = "gray90"
           )
         ) +
      
         # 4. Labels and layout
         xlab("Position") +
         ggtitle(paste("Copy Number", My_paire[i])) +
      
         facet_grid(reorder(Type, ordre_split) ~ reorder(CONTIG, as.numeric(chr_pos)),
                    scales = "free_x", space = "free_x", switch = "x") +
      
         # 5. Theme
         theme(
           axis.text.x = element_text(size = 3, face = "bold", color = "gray",angle = 45, hjust = 0.5,vjust = 0),
           axis.title.x = element_text(size = 18, face = "bold", color = "gray"),
           plot.title = element_text(color = "black", size = 12, face = "bold"),
           axis.text.y = element_blank(),
           axis.title.y.left = element_text(size = 18, face = "bold", color = "black",angle = 90),
           panel.background = element_blank(),
           axis.title = element_blank(),
           axis.ticks.y = element_blank(),
           plot.margin = unit(c(1, 1, 0, 1), "cm"),
           strip.text.x = element_text(angle = 90, hjust = 0.5, size = 10, face = "bold"),
           strip.text.y.left =  element_text(angle = 90, hjust = 0.5, size = 12, face = "bold"),
           strip.background = element_blank()
         )
      
      
        ggsave(paste0(hhg,"MAP_SEG_paire",My_paire[i],".pdf"), plot = x2, width = 16, height = 3, limitsize = FALSE, device = 'pdf')
      
      
     }
     
      my_ext1_select23 <- my_ext1_select23 %>%
        complete(Type = c("ADC","SCC"),
                 sign = c("Amp","Del"),
                 fill = list(Total_CNV = 0, share = 0, commun = 0, Paire = 0))


     
     
     if(nrow(my_ext1_select23)==0){
       
       nSCCdel<-nrow(split_data_select[which(split_data_select$Type=="SCC"& split_data_select$sign=="Del"),])
       nADCdel<-nrow(split_data_select[which(split_data_select$Type=="ADC"& split_data_select$sign=="Del"),])
       nSCCamp<-nrow(split_data_select[which(split_data_select$Type=="SCC"& split_data_select$sign=="Amp"),])
       nADCamp<-nrow(split_data_select[which(split_data_select$Type=="ADC"& split_data_select$sign=="Amp"),])
       my_ext1_select23<-data.frame(Type=c("ADC","SCC","ADC","SCC"),Total_CNV=c(nADCamp,nSCCamp,nADCdel,nSCCdel),sign=c("Amp","Amp","Del","Del"),
                                    share=c(0,0,0,0),commun=c(0,0,0,0),Paire=c(0,0,0,0))
     }
     
     
     
     if( nrow(my_ext1_select23[which(my_ext1_select23$Type=="ADC"& my_ext1_select23$sign=="Amp"),])<1){
       nADCamp<-sum(split_data_select[which(split_data_select$Type=="ADC"& split_data_select$sign=="Amp"),])
       liness<-data.frame(Type=c("ADC"),Total_CNV=c(nADCamp),sign=c("Amp"),share=c(0),commun=c(0),Paire=c(0))
       my_ext1_select23<-rbind(my_ext1_select23,liness)
     }
     
     if( nrow(my_ext1_select23[which(my_ext1_select23$Type=="SCC"& my_ext1_select23$sign=="Amp"),])<1){
       nSCCamp<-nrow(split_data_select[which(split_data_select$Type=="SCC"& split_data_select$sign=="Amp"),])
       liness<-data.frame(Type=c("SCC"),Total_CNV=c(nADCamp),sign=c("Amp"),share=c(0),commun=c(0),Paire=c(0))
       my_ext1_select23<-rbind(my_ext1_select23,liness)
     }
     
     if( nrow(my_ext1_select23[which(my_ext1_select23$Type=="ADC"& my_ext1_select23$sign=="Del"),])<1){
       nADCdel<-nrow(split_data_select[which(split_data_select$Type=="ADC"& split_data_select$sign=="Del"),])
       liness<-data.frame(Type=c("ADC"),Total_CNV=c(nADCamp),sign=c("Del"),share=c(0),commun=c(0),Paire=c(0))
       my_ext1_select23<-rbind(my_ext1_select23,liness)
     }
     
     if( nrow(my_ext1_select23[which(my_ext1_select23$Type=="SCC"& my_ext1_select23$sign=="Del"),])<1){
       nSCCdel<-nrow(split_data_select[which(split_data_select$Type=="SCC"& split_data_select$sign=="Del"),])
       liness<-data.frame(Type=c("SCC"),Total_CNV=c(nADCamp),sign=c("Del"),share=c(0),commun=c(0),Paire=c(0))
       my_ext1_select23<-rbind(my_ext1_select23,liness)
     }
     
     
     
     if (i == y){
       
       my_ext1_select23$Paire<-My_paire[i]
       Resultes<-rbind(Resultes,my_ext1_select23)
       
     }else{
       my_ext1_select23$Paire<-"CTR"
       Resultes_ctr<-rbind(Resultes_ctr,my_ext1_select23)
       
     }
     
     
   }
 }
#############################################################
#############################################################
#############################################################
#############################################################
# 
# 
 Loas_files_POINT$CONTIG[!grepl("chr.+",Loas_files_POINT$CONTIG)]<-paste0("chr",Loas_files_POINT$CONTIG[!grepl("chr.+",Loas_files_POINT$CONTIG)])

 for ( i in 1:length(My_paire)){
   My_paire[i]
     Loas_files_select<-Loas_files[which(Loas_files$Paire==My_paire[i] ),]
     Loas_files_select$Type

     my_ext1<-Loas_files_POINT[which(Loas_files_POINT$Paire_num==My_paire[i]),]

     my_ext1$chr_pos<-gsub("chr(.+)","\\1",my_ext1$CONTIG)
     my_ext1$chr_pos[which(my_ext1$chr_pos=="X")]<-23
     Loas_files_select$chr_pos<-gsub("chr(.+)","\\1",Loas_files_select$CONTIG)
     Loas_files_select$chr_pos[which(Loas_files_select$chr_pos=="X")]<-23


     my_ext1$Type[which(my_ext1$Type=="CIN3")]<-"SCC"
     my_ext1$Type[which(my_ext1$Type=="SCC2")]<-"SCC"
     my_ext1$Type[which(my_ext1$Type=="ADC2")]<-"ADC"
     my_ext1$Type[which(my_ext1$Type=="SCC3")]<-"SCC"

     Loas_files_select$Type[which(Loas_files_select$Type=="SCC3")]<-"SCC"
     Loas_files_select$Type[which(Loas_files_select$Type=="CIN3")]<-"SCC"
     Loas_files_select$Type[which(Loas_files_select$Type=="SCC2")]<-"SCC"
     Loas_files_select$Type[which(Loas_files_select$Type=="ADC2")]<-"ADC"


     Loas_files_select$MEAN_LOG2_COPY_RATIO[Loas_files_select$MEAN_LOG2_COPY_RATIO>1]<-1
     Loas_files_select$MEAN_LOG2_COPY_RATIO[Loas_files_select$MEAN_LOG2_COPY_RATIO< -1]<- -1
     Loas_files_select<-unique(Loas_files_select)

     x1 <- ggplot(data = my_ext1, aes(x = LOC, y = LOG2_COPY_RATIO)) +
       ylim(c(-1.1, 1.1)) +
       theme_classic() +

       # Scatter and segments
       geom_point(col = "grey", size = 0.5, alpha = 0.05) +
       geom_segment(data = Loas_files_select, aes(x = START, xend = END,
                                                  y = MEAN_LOG2_COPY_RATIO, yend = MEAN_LOG2_COPY_RATIO, col = Type),
                    linewidth = 2) +

       # Colors and axis labels
       scale_color_manual(values = c("indianred1", "#EE9A49", "steelblue3", "#f781bf")) +
       ylab("Log2 Copy Number") +
       geom_hline(yintercept = hhg) +
       geom_hline(yintercept = -hhg) +
       xlab("Position") +
       ggtitle(paste("Copy Number", My_paire[i])) +

       # Faceting only by CONTIG
       facet_grid(reorder(Type,Type_order) ~ reorder(CONTIG, as.numeric(chr_pos)),
                  scales = "free_x", space = "free_x", switch = "x") +

       # Theming
       theme(
         axis.text.x = element_text(size = 3, face = "bold", color = "gray",angle = 45, hjust = 0.5,vjust = 0),
         axis.title.x = element_text(size = 18, face = "bold", color = "gray"),
         plot.title = element_text(color = "black", size = 20, face = "bold"),
         axis.text.y = element_text(size = 18, face = "bold", color = "black"),
         axis.title.y.left = element_text(size = 18, face = "bold", color = "black"),
         panel.background = element_blank(),

         axis.title = element_blank(),
         plot.margin = unit(c(1, 1, 0, 1), "cm"),

         # Rotate CONTIG strip text
         strip.text.x = element_text(angle = 90, hjust = 0.5, size = 10, face = "bold"),
         strip.background = element_blank()
       )


   # dev.off()
   myOutFile_plot <- paste0(hhg,"new_",My_paire[i],"_plot_.png")
   ggsave(myOutFile_plot, width=16, height=8, limitsize = FALSE, device='png')



 }


############################################################################
#############################################################
 
 
Loas_files <- read.delim(file = "Full_CNVs.txt", header = TRUE, sep = "\t")
Loas_files$sig<-"AMP"
Loas_files$sig[which(Loas_files$CopyNumber<0)]<-"DEL"

Loas_files$Paire_num<-gsub("Paire_(\\d+)","\\1",Loas_files$Paire)
Loas_files$Paire_num<-as.numeric(Loas_files$Paire_num)
datas_cliniques <- read.delim(file = "./Datas_cliniques.txt", header = TRUE, sep = "\t")
Loas_files_annotated <- merge(datas_cliniques, Loas_files, by.x = "Paire", by.y = "Paire_num")
Loas_files_annotated<-Loas_files_annotated[Loas_files_annotated$Threshold==1,]


Loas_files<-read.delim(file = "Test_COMBINED.txt")
Loas_files$Type_order<-1
Loas_files$Type[Loas_files$Type=="CIN3"]<-"SCC"
Loas_files$Type[Loas_files$Type=="ADC2"]<-"ADC"  
Loas_files$Type[Loas_files$Type=="SCC2"]<-"SCC"
Loas_files$Type[Loas_files$Type=="SCC3"]<-"SCC" 
Loas_files$Type_order[which(Loas_files$Type=="ADC")]<-2
Loas_files$Type_order[which(Loas_files$Type=="SCC")]<-3
Loas_files<-unique(Loas_files)
Loas_files$LOC<-Loas_files$START+((Loas_files$END-Loas_files$START)/2)
Loas_files$ordre<-gsub('chr(.+)','\\1',Loas_files$CONTIG)
Loas_files$ordre[which(Loas_files$CONTIG=="X")]<-22
Loas_files$Paire_num<-gsub("Paire_(\\d+)","\\1",Loas_files$Paire)
Loas_files$Paire<-as.numeric(Loas_files$Paire_num)
Loas_files <- Loas_files[!(Loas_files$Paire %in% c(2,5, 7, 9,10,15)), ]
Loas_files$Paire.y<-NULL
Loas_files<-unique(Loas_files)
Loas_files$size<-Loas_files$END-Loas_files$START
Loas_files_seg<-Loas_files[!Loas_files$NUM_POINTS_COPY_RATIO<5,]
Loas_files_seg<-Loas_files_seg[which(Loas_files_seg$MEAN_LOG2_COPY_RATIO>1|Loas_files_seg$MEAN_LOG2_COPY_RATIO< -1),]

Loas_files_annotated2<-merge(Loas_files_annotated,Loas_files_seg, by.x=c("Paire","Size","Cancer"),by.y = c("Paire","size","Type"))

#
Ful_datas<-read.delim(file = "combined_files.txt", header = T, sep = "\t")
Ful_datas$Paire_num<-gsub("Paire_(\\d+)","\\1",Ful_datas$PAIR)
Ful_datas$Paire<-as.numeric(Ful_datas$Paire_num)
Ful_datas<-Ful_datas[!Ful_datas$Num_Probes<5,]
Ful_datas<-Ful_datas[which(Ful_datas$Segment_Mean>0.3|Ful_datas$Segment_Mean< -0.3),]
Ful_datas$Size<-Ful_datas$end-Ful_datas$start
Ful_datas$TYPE[Ful_datas$TYPE=="CIN"]<-"SCC"
Ful_datas$TYPE[Ful_datas$TYPE=="ADC2"]<-"ADC"  
Ful_datas$TYPE[Ful_datas$TYPE=="SCC2"]<-"SCC"
Ful_datas$TYPE[Ful_datas$TYPE=="SCC3"]<-"SCC" 
Ful_datas$sig<-"AMP"
Ful_datas$sig[which(Ful_datas$Segment_Mean<0)]<-"DEL"

Ful_datas <- merge(datas_cliniques, Ful_datas, by.x = "Paire", by.y = "Paire")

Loas_files_annotated2 <- Ful_datas %>%
  mutate(n_UD = "1") %>%  # Valeur par défaut
  group_by(Paire,chr, start,sig) %>%
  mutate(
    n_comm2 = n_distinct(TYPE),
    n_UD = if_else(n_comm2 > 1, paste0(Paire, start), n_UD)
  ) %>%
  ungroup() %>%
  group_by(Paire,chr, end,sig) %>%
  mutate(
    n_comm3 = n_distinct(TYPE),
    n_UD = if_else(n_comm3 > 1, paste0(Paire, end), n_UD)
  ) %>%
  ungroup() %>%
  group_by(Paire, start_gene, end_gene,sig) %>%
  mutate(
    n_comm = n_distinct(TYPE),
    n_UD = if_else(n_comm > 1, paste0(Paire, start_gene, end_gene), n_UD)
  ) %>%
  ungroup() %>%
  mutate(
    annotation_genes = case_when(
      genes >= 500 & genes < 5000 ~ "500-2500",
      genes >= 50 & genes < 500 ~ "50-500",
      genes >= 5 & genes < 50 ~ "5-50",
      TRUE ~ "1-10"
    )
  )


Loas_files_annotated2_commun <- Loas_files_annotated2 %>%
  filter(n_comm == 2 | n_comm2 == 2 | n_comm3 == 2)

Loas_files_annotated2_commun<-Loas_files_annotated2_commun%>%
  group_by(n_UD)%>%
  mutate(genes=mean(genes))%>%
  ungroup()%>%
  mutate(
    annotation_genes = case_when(
      genes >= 500 & genes < 5000 ~ "500-2500",
      genes >= 50 & genes < 500 ~ "50-500",
      genes >= 5 & genes < 50 ~ "5-50",
      TRUE ~ "1-10"
    )
  )

Loas_files_annotated2_commun_test <- Loas_files_annotated2_commun %>%
  filter( n_comm == 1)

# Filtrage des autres cas (aucun == 2)
Loas_files_annotated2_nocommun <- Loas_files_annotated2 %>%
  filter(!(n_comm == 2 | n_comm2 == 2 | n_comm3 == 2))


Loas_files_annotated2_commun$TYPE<-"Commun"
Loas_files_annotated2_commun<-Loas_files_annotated2_commun%>%
  dplyr::select(TYPE,sig,annotation_genes,n_UD)%>%
  distinct()%>%
  dplyr::select(TYPE,sig,annotation_genes)

Loas_files_annotated2_nocommun<-Loas_files_annotated2_nocommun%>%
  dplyr::select(TYPE,sig,annotation_genes)


Loas_files_annotated3<-rbind(Loas_files_annotated2_commun,Loas_files_annotated2_nocommun)

Loas_files_annotated3.5<-Loas_files_annotated3%>%
  group_by(TYPE,sig,annotation_genes)%>%
  summarise(count=n())
head(Loas_files_annotated3.5)

Loas_files_annotated3.5_filled <- Loas_files_annotated3.5 %>%
  ungroup() %>%
  complete(
    #genes = 1:12,
    TYPE  = c("ADC", "SCC","Commun"),
    sig     = c("AMP", "DEL"),
    annotation_genes = c("1-10", "5-50","50-500","500-2500"),
    fill    = list(count = 0)
  )

max_count<-max(Loas_files_annotated3.5_filled$count)

Loas_files_annotated3.5_filled$type_num<-1
Loas_files_annotated3.5_filled$type_num[Loas_files_annotated3.5_filled$TYPE=="SCC"]<-2
Loas_files_annotated3.5_filled$type_num[Loas_files_annotated3.5_filled$TYPE=="Commun"]<-3


sum(Loas_files_annotated3.5_filled$count)

sum(Loas_files_annotated3.5_filled$count[Loas_files_annotated3.5_filled$TYPE=="Commun"])

x2a <- ggplot(data = Loas_files_annotated3.5_filled,
              aes(x = annotation_genes, 
                  y = count , 
                  fill = sig )) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  geom_hline(aes(yintercept = 0)) +
  ylim(0,as.numeric(max_count))+
  scale_fill_manual(
    values = c(
      "Amp" = "#9ACBB3",
      "Del" = "#C89CCB",
      "AMP" = "#9ACBB3",
      "DEL" = "#C89CCB"
    )
  ) + 
  theme(
    axis.text.x = element_text(size = 16, face = "bold", color = "gray", vjust = 0,angle = 45,hjust = 1),
    axis.title.x = element_blank(),
    
    panel.background = element_blank(),
    
    axis.ticks.x = element_blank(),
    axis.title = element_blank())+
  facet_wrap(. ~ reorder(TYPE,type_num), nrow = 1,scale="free_x", strip.position = "bottom") +
  #scale_y_continuous(trans = "reverse") +
  theme(
    
    plot.margin = unit(c(1, 1, 1.5, 1.5), "cm"),
    strip.text = element_text(size = 18, face = "bold", color = "gray", vjust = 1,angle = 0,hjust = 0.5)
  )

ggsave(paste0("Genes_commun_select.pdf"), plot = x2a, width = 12, height = 8, limitsize = FALSE, device = 'pdf', dpi = 300)













################################################################################################################################################################"

Loas_files<-read.delim(file = "Test_COMBINED.txt")
Loas_files$Type_order<-1
Loas_files$Type[Loas_files$Type=="CIN3"]<-"SCC"
Loas_files$Type[Loas_files$Type=="ADC2"]<-"ADC"  
Loas_files$Type[Loas_files$Type=="SCC2"]<-"SCC"
Loas_files$Type[Loas_files$Type=="SCC3"]<-"SCC" 

Loas_files$Type_order[which(Loas_files$Type=="ADC")]<-2
Loas_files$Type_order[which(Loas_files$Type=="SCC")]<-3
Loas_files<-unique(Loas_files)



Loas_files$LOC<-Loas_files$START+((Loas_files$END-Loas_files$START)/2)


Loas_files$chr_pos<-gsub('chr(.+)','\\1',Loas_files$CONTIG)
Loas_files$chr_pos[which(Loas_files$CONTIG=="X")]<-23


Loas_files$Paire_num<-gsub("Paire_(\\d+)","\\1",Loas_files$Paire)
# Loas_files <- merge(Loas_files, Loas_samples, by.x = "Paire_num", by.y = "Old")
Loas_files$Paire<-as.numeric(Loas_files$Paire_num)


Loas_files <- Loas_files[!(Loas_files$Paire %in% c(2,5, 7, 9,10,15)), ]
Loas_files$Paire.y<-NULL
Loas_files<-unique(Loas_files)


Loas_files$size<-Loas_files$END-Loas_files$START
Loas_files<-Loas_files[!Loas_files$NUM_POINTS_COPY_RATIO<5,]
# 
Ful_datas$end
# 
Ful_datas<-read.delim(file = "/media/florian2/T7/To_send/CNVs/combined_files.txt", header = T, sep = "\t")
Ful_datas$Paire_num<-gsub("Paire_(\\d+)","\\1",Ful_datas$PAIR)
Ful_datas$Paire<-as.numeric(Ful_datas$Paire_num)
Ful_datas<-Ful_datas[!Ful_datas$Num_Probes<5,]
Ful_datas$Size<-Ful_datas$end-Ful_datas$start
Ful_datas$TYPE[Ful_datas$TYPE=="CIN"]<-"SCC"
Ful_datas$TYPE[Ful_datas$TYPE=="ADC2"]<-"ADC"  
Ful_datas$TYPE[Ful_datas$TYPE=="SCC2"]<-"SCC"
Ful_datas$TYPE[Ful_datas$TYPE=="SCC3"]<-"SCC" 
Ful_datas$sig<-"AMP"
Ful_datas$sig[which(Ful_datas$Segment_Mean<0)]<-"DEL"

Ful_datas$chr_pos<-gsub("chr(.+)","\\1",Ful_datas$chr)
Ful_datas$chr_pos[which(Ful_datas$chr_pos=="X")]<-23

datas_cliniques <- read.delim(file = "Datas_cliniques.txt", header = TRUE, sep = "\t")
Loas_files_annotated <- merge(datas_cliniques, Ful_datas, by.x = "Paire", by.y = "Paire")


Resultes_ctr<- data.frame(Type = integer(),Total_CNV = numeric(),sig = integer(),share = numeric(),commun = numeric(),Paire = integer(),ids= integer(), stringsAsFactors = FALSE)
Resultes<- data.frame(Type = integer(),Total_CNV = numeric(),sig = integer(),share = numeric(),commun = numeric(),Paire = integer(), ids=integer(), stringsAsFactors = FALSE)
My_paire=unique(Loas_files_annotated$Paire)

My_paire

for ( i in 1:length(My_paire)){
  
  for ( y in 1:length(My_paire)){
    # 
      # i<-6
      # y<-6
    My_paire[i]
    My_paire[y]
    Loas_files_annotated$TYPE
    
    
    oLoas_files_select_i<-Loas_files[which(Loas_files$Paire==My_paire[i] & Loas_files$Type=="ADC"),]
    oLoas_files_select_y<-Loas_files[which(Loas_files$Paire==My_paire[y]& Loas_files$Type=="SCC"),]
    oLoas_files_select<-rbind(oLoas_files_select_i,oLoas_files_select_y)
    
    
    Loas_files_select_i<-Loas_files_annotated[which(Loas_files_annotated$Paire==My_paire[i] & Loas_files_annotated$TYPE=="ADC"),]
    Loas_files_select_y<-Loas_files_annotated[which(Loas_files_annotated$Paire==My_paire[y]& Loas_files_annotated$TYPE=="SCC"),]
    Loas_files_select<-rbind(Loas_files_select_i,Loas_files_select_y)
    
    Loas_files_select_select<-Loas_files_select[which(Loas_files_select$Segment_Mean>0.3|Loas_files_select$Segment_Mean< -0.3),]
    
    library(dplyr)
    
    my_ext1_select23 <- Loas_files_select_select %>% 
      # per-Type: unique (START, END) count 
      group_by(TYPE) %>% 
      mutate(Total_CNV = n_distinct(start, end)) %>% 
      ungroup() %>% 
    # global: sum one Total_CNV per Type (no undercount when counts are equal) 
      mutate(Full_CNV = sum(Total_CNV[!duplicated(TYPE)]))%>% 
      # per-region: how many Types have this (PosSTART, PosSTOP) 
      group_by(start_gene, end_gene,sig) %>% 
      mutate(ncount = n_distinct(TYPE)) %>% 
      ungroup() %>% 
      group_by(chr, start,sig) %>% 
      mutate(ncount2 = n_distinct(TYPE)) %>% 
      ungroup() %>% 
      group_by(chr, end,sig) %>% 
      mutate(ncount3 = n_distinct(TYPE)) %>% 
      ungroup() %>% 
      # flag shared regions (present in >= 2 Types) 
      mutate(shared_flag = as.integer(ncount >= 2 |ncount2 >= 2 | ncount3 >= 2)) %>% 
      # per-Type: number of shared regions #group_by(TYPE) %>% 
      mutate(share = sum(shared_flag)) %>% ungroup() %>% 
      # % shared relative to non-shared count 
      mutate(commun = ((share/2) / (Full_CNV - (share/2)) * 100)) %>% 
      # keep one row per grouping key (matches your original intent) 
      group_by(TYPE, Total_CNV, sig, share, commun, Paire) %>% 
      summarise(.groups = "drop")
      
    
    my_ext1_select23$Paire<-"CTR"
    
    my_ext1_select23$ids<-paste(sep = "_",i,y)

    
    if (i == y){
      
      my_ext1_select23$Paire<-My_paire[i]
      Resultes<-rbind(Resultes,my_ext1_select23)
      
    }else{
      
      Resultes_ctr<-rbind(Resultes_ctr,my_ext1_select23)
      
    }
    
  }
  
}

# colnames(Resultes)

Resultes_filled <- Resultes %>%
  ungroup() %>%
  complete(
    Paire = 1:18,
    TYPE  = c("ADC", "SCC","Commun"),
    sig     = c("AMP", "DEL"),
    fill    = list(Total_CNV = 0,share = 0,commun = 0)
  )







datas_cliniques <- read.delim(file = "Datas_cliniques.txt", header = TRUE, sep = "\t")
Resultes2 <- merge(datas_cliniques, Resultes_filled, by.x = "Paire", by.y = "Paire")
corrected_Loas_files2<-unique(Resultes2[,c(5,11,17)])
df <- data.frame(x = factor(), y = numeric())

corrected_Loas_files2$ctr<-max(Resultes_ctr$commun)

Result2a<-corrected_Loas_files2[which(corrected_Loas_files2$FIGO=="1a1"),]
Result2b<-corrected_Loas_files2[which(corrected_Loas_files2$FIGO=="1b1"),]
Result2c<-corrected_Loas_files2[which(corrected_Loas_files2$FIGO=="1b2"),]

max_prop<-as.numeric(max(corrected_Loas_files2$commun))



df <- data.frame(x = factor(), y = numeric())

# Start the plot
x10<-ggplot(df, aes(x, y)) +
  geom_bar(stat = "identity") +  # Won't actually draw anything because df is empty
  theme_void() +
  labs(y = "Proportion of common CNVs") +
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
                                  y = commun)) +
  geom_bar(stat = "identity", fill = "gray") +
  ylim(0,max_prop)+
  geom_hline(aes(yintercept = ctr), linetype = "dashed") +
  geom_hline(aes(yintercept = 0)) +
  geom_hline(yintercept =max(Resultes_ctr$commun), linetype = "dashed")+
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
                                   y = commun)) +
  geom_bar(stat = "identity", fill = "gray") +
  ylim(0,max_prop)+
  geom_hline(aes(yintercept = 0)) +
  geom_hline(yintercept =max(Resultes_ctr$commun), linetype = "dashed") +
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
Result2c$commun
x1c <- ggplot(data = Result2c, aes(x = as.character(Result2c$New),
                                   y = commun)) +
  geom_bar(stat = "identity", fill = "gray") +
  ylim(0,max_prop)+
  geom_hline(aes(yintercept = max(Resultes_ctr$commun)), linetype = "dashed") +
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

Loas_files_epure <- Loas_files %>%
  filter(NUM_POINTS_COPY_RATIO >= 5) %>%
  filter(MEAN_LOG2_COPY_RATIO > 0.3 | MEAN_LOG2_COPY_RATIO < -0.3) %>%
  mutate(sign = ifelse(MEAN_LOG2_COPY_RATIO > 0, "Amp", "Del")) %>%
  group_by(Type, Paire, sign) %>%
  summarise(Count = n(), .groups = "drop")

Loas_files_epure_link <- merge(datas_cliniques, Loas_files_epure, by.x = "Paire", by.y = "Paire")
colnames(Loas_files_epure_link)
for (i in 1:12) {
  for (y in c("ADC", "SCC")) {
    
    # Vérifie si une ligne "Amp" existe
    if (nrow(subset(Loas_files_epure_link, New == i & Type == y & sign == "Amp")) < 1) {
      liny <- datas_cliniques[datas_cliniques$New == i, ] %>%
        dplyr::mutate(Type = y, sign = "Amp", Count = 0)
      
      Loas_files_epure_link <- rbind(Loas_files_epure_link, liny)
    }
    
    # Vérifie si une ligne "Del" existe
    if (nrow(subset(Loas_files_epure_link, New == i & Type == y & sign == "Del")) < 1) {
      liny <- datas_cliniques[datas_cliniques$New == i, ] %>%
        dplyr::mutate(Type = y, sign = "Del", Count = 0)
      
      Loas_files_epure_link <- rbind(Loas_files_epure_link, liny)
    }
    
  }
}

sum(Loas_files_epure_link$Count[Loas_files_epure_link$Type=="SCC"])

Result3a<-Loas_files_epure_link[which(Loas_files_epure_link$FIGO=="1a1"),]
Result3b<-Loas_files_epure_link[which(Loas_files_epure_link$FIGO=="1b1"),]
Result3c<-Loas_files_epure_link[which(Loas_files_epure_link$FIGO=="1b2"),]

Loas_files_epure_link<-Loas_files_epure_link%>%
  group_by(New)%>%
  mutate(max=sum(Count))

max_count<-max(Loas_files_epure_link$max)
# Start the plot
x20<-ggplot(df, aes(x, y)) +
  geom_bar(stat = "identity") +  # Won't actually draw anything because df is empty
  theme_void() +
  labs(y = "Number of CNVs") +
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
Result3a$Count

x2a <- ggplot(data = Result3a,
              aes(x = Type, 
                  y = Count , 
                  fill = sign )) +
  geom_bar(stat = "identity", color = "black") +
  geom_hline(aes(yintercept = 0)) +
  ylim(as.numeric(max_count),0)+
  scale_fill_manual(
    values = c(
      "Amp" = "#9ACBB3",
      "Del" = "#C89CCB",
      "AMP" = "#9ACBB3",
      "DEL" = "#C89CCB"
    )
  ) + 
  facet_wrap(~ reorder((New),(New)), nrow = 1,scale="free_x") +
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
              aes(x = Type, 
                  y = Count, 
                  fill = sign)) +
  geom_bar(stat = "identity", color = "black") +
  geom_hline(aes(yintercept = 0)) +
  ylim(as.numeric(max_count),0)+
  scale_fill_manual(
    values = c(
      "Amp" = "#9ACBB3",
      "Del" = "#C89CCB",
      "AMP" = "#9ACBB3",
      "DEL" = "#C89CCB"
    )
  ) + 
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
              aes(x = Type, 
                  y = Count, 
                  fill = sign)) +
  geom_bar(stat = "identity", color = "black") +
  geom_hline(aes(yintercept = 0)) +
  ylim(as.numeric(max_count),0)+
  scale_fill_manual(
    values = c(
      "Amp" = "#9ACBB3",
      "Del" = "#C89CCB",
      "AMP" = "#9ACBB3",
      "DEL" = "#C89CCB"
    )
  ) + 
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

ggsave(paste0("_newFig2_CNV_commun_select22.pdf"), plot = final_plot, width = 12, height = 8, limitsize = FALSE, device = 'pdf', dpi = 300)



########################
############################
################################

y_limits <- c(1e2, max(Loas_files$size))  # adjust min/max according to your data
Ful_datas<-Ful_datas[which(Ful_datas$Segment_Mean>0.3|Ful_datas$Segment_Mean< -0.3),]
Loas_files<-Loas_files[which(Loas_files$MEAN_LOG2_COPY_RATIO>0.3|Loas_files$MEAN_LOG2_COPY_RATIO< -0.3),]

Loas_files_annotated <- merge(datas_cliniques, Ful_datas, by.x = "Paire", by.y = "Paire")




Loas_files_annotated<-Loas_files_annotated%>%
  mutate(n_UD="1")%>%
  group_by(Paire,chr, start,sig) %>%
  mutate(
    n_comm2 = n_distinct(TYPE),
    n_UD = if_else(n_comm2 > 1, paste0(Paire, start), n_UD)
  ) %>%
  ungroup() %>%
  group_by(Paire,chr, end,sig) %>%
  mutate(
    n_comm3 = n_distinct(TYPE),
    n_UD = if_else(n_comm3 > 1, paste0(Paire, end), n_UD)
  ) %>%
  ungroup() %>%
  group_by(Paire, start_gene, end_gene,sig) %>%
  mutate(
    n_comm = n_distinct(TYPE),
    n_UD = if_else(n_comm > 1, paste0(Paire, start_gene, end_gene), n_UD)
  ) %>%
  mutate(n_comm=n_distinct(TYPE))%>%
  
  filter(n_comm==2|n_comm2==2|n_comm3==2)%>%
  ungroup()%>%
  group_by(n_UD)%>%
  mutate(size=mean(Size),ratio=mean(Segment_Mean))%>%
  dplyr::select(start,end,Num_Probes,Size,size,Segment_Mean,ratio)%>%
  distinct()


Loas_files_epure1 <- Loas_files 
Loas_files_epure1$selection<-paste0(Loas_files_epure1$START,
                                    Loas_files_epure1$END,
                                    Loas_files_epure1$NUM_POINTS_COPY_RATIO,
                                    Loas_files_epure1$MEAN_LOG2_COPY_RATIO)

Loas_files_annotated$selection<-paste0(Loas_files_annotated$start,
                                       Loas_files_annotated$end,
                                       Loas_files_annotated$Num_Probes,
                                       Loas_files_annotated$Segment_Mean)

Loas_files_epure12<-Loas_files_epure1[!Loas_files_epure1$selection%in%Loas_files_annotated$selection,]

Loas_files_epure2 <- Loas_files_epure12 %>%
  mutate(sign = ifelse(MEAN_LOG2_COPY_RATIO > 0, "Amp", "Del")) %>%
  dplyr::select(Type, Paire, sign, size)

Loas_files_epure_link2 <- merge(datas_cliniques, Loas_files_epure2, by.x = "Paire", by.y = "Paire")


TMB_plot1b <- ggplot(Loas_files_epure2[Loas_files_epure2$Type=="ADC",], aes(x = sign, y = size,fill = sign)) +
  geom_violin(alpha = 0.6) +  # Violin plot for distributione
  #geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(shape = 16, size = 3, width = 0.15, alpha = 0.7) +  # Add points for individual values
  geom_hline(yintercept = 1000000,linetype = "dashed")+
  scale_y_log10(limits = y_limits) +  # same y-axis for all plots
  scale_fill_brewer(palette = "Set2") +  # Nice color palette
  theme_void() +
  scale_fill_manual(
    values = c(
      "Amp" = "#9ACBB3",
      "Del" = "#C89CCB",
      "AMP" = "#9ACBB3",
      "DEL" = "#C89CCB"
    )
  ) + 
  theme(
    axis.text.x = element_text( size = 14,angle = 0,hjust = 1,vjust = 1),
    axis.text.y =  element_text( size = 14),
    axis.title.x = element_text(size = 16, face = "bold",angle = 0),
    axis.title.y = element_text(size = 16, face = "bold",angle = 90,vjust = 1),
    legend.position = "top"
  ) +
  labs(y = "Size of CNVs (bp)", x = "ADC")+
  stat_compare_means(method = "anova",label="p.format",label.x.npc = "right") 

TMB_plot2b <- ggplot(Loas_files_epure2[Loas_files_epure2$Type=="SCC",], aes(x = sign, y = size,fill = sign)) +
  geom_violin(alpha = 0.6) +  # Violin plot for distributione
  #geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(shape = 16, size = 3, width = 0.15, alpha = 0.7) +  # Add points for individual values
  geom_hline(yintercept = 1000000,linetype = "dashed")+
  scale_y_log10(limits = y_limits) +  # same y-axis for all plots
  scale_fill_brewer(palette = "Set2") +  # Nice color palette
  theme_void() +
  scale_fill_manual(
    values = c(
      "Amp" = "#9ACBB3",
      "Del" = "#C89CCB",
      "AMP" = "#9ACBB3",
      "DEL" = "#C89CCB"
    )
  ) + 
  theme(
    axis.title.y = element_blank(),  # remove y-axis title
    axis.text.y = element_blank(),   # remove y-axis text
    axis.ticks.y = element_blank(),  # remove y-axis ticks
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    legend.position = "top"
  ) +
  labs(y = "Size of CNVs (bp)", x = "SCC")+
  stat_compare_means(method = "anova",label="p.format",label.x.npc = "right") 

Loas_files_annotated <- merge(datas_cliniques, Ful_datas, by.x = "Paire", by.y = "Paire")

Loas_files_annotated3<-Loas_files_annotated%>%
  mutate(n_UD="1")%>%
  group_by(Paire,chr, start,sig) %>%
  mutate(
    n_comm2 = n_distinct(TYPE),
    n_UD = if_else(n_comm2 > 1, paste0(Paire, start), n_UD)
  ) %>%
  ungroup() %>%
  group_by(Paire,chr, end,sig) %>%
  mutate(
    n_comm3 = n_distinct(TYPE),
    n_UD = if_else(n_comm3 > 1, paste0(Paire, end), n_UD)
  ) %>%
  ungroup() %>%
  group_by(Paire, start_gene, end_gene,sig) %>%
  mutate(
    n_comm = n_distinct(TYPE),
    n_UD = if_else(n_comm > 1, paste0(Paire, start_gene, end_gene), n_UD)
  ) %>%
  mutate(n_comm=n_distinct(TYPE))%>%
  filter(n_comm==2|n_comm2==2|n_comm3==2)%>%
  ungroup()%>%
  group_by(n_UD)%>%
  mutate(size=mean(Size),ratio=mean(Segment_Mean))%>%
  ungroup()%>%
  dplyr::select(sig,size)%>%
  distinct()

Loas_files_annotated3<-Loas_files_annotated%>%
  mutate(sign = ifelse(Segment_Mean > 0, "Amp", "Del")) %>%
  group_by(Paire,start_gene,end_gene)%>%
  mutate(n_comm=n_distinct(TYPE))%>%
  filter(n_comm==2)%>%
  ungroup()%>%
  group_by(Paire,start_gene,end_gene)%>%
  mutate(size=mean(Size))%>%
  select(sign,size)

TMB_plot3b <- ggplot(Loas_files_annotated3, aes(x = sig, y = size, fill = sig)) +
  geom_violin(alpha = 0.6) +
  geom_jitter(shape = 16, size = 3, width = 0.15, alpha = 0.7) +
  geom_hline(yintercept = 1000000, linetype = "dashed") +
  scale_y_log10(limits = y_limits) +  # same y-axis for all plots
  scale_fill_manual(
    values = c(
      "Amp" = "#9ACBB3",
      "Del" = "#C89CCB",
      "AMP" = "#9ACBB3",
      "DEL" = "#C89CCB"
    )
  ) +
  theme_void() +  # base theme with no axes
  theme(
    axis.title.y = element_blank(),  # remove y-axis title
    axis.text.y = element_blank(),   # remove y-axis text
    axis.ticks.y = element_blank(),  # remove y-axis ticks
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    legend.position = "top"
  ) +
  labs(x = "Commun") +
  stat_compare_means(method = "anova", label = "p.format", label.x.npc = "right")

combi_row <- TMB_plot1b + TMB_plot2b + TMB_plot3b + 
  plot_layout(ncol = 3, guides = "collect") &    # & applies theme to all subplots
  theme(legend.position = "top")                 # place the collected legend on top

# Save the combined plot
ggsave(paste0("lol_plot2_plot_by_CNV_select.pdf"), combi_row, width = 12, height = 10, dpi = 30






################################" VOlcano ####################

Loas_files_epure2 <- Loas_files_epure12 %>%
  mutate(sign = ifelse(MEAN_LOG2_COPY_RATIO > 0, "Amp", "Del")) %>%
  dplyr::select(Type, Paire, sign, size,MEAN_LOG2_COPY_RATIO)


Loas_files_annotated <- merge(datas_cliniques, Ful_datas, by.x = "Paire", by.y = "Paire")

Loas_files_annotated1<-Loas_files_annotated%>%
  group_by(Paire,start_gene,end_gene)%>%
  mutate(n_comm=n_distinct(TYPE))%>%
  filter(n_comm==2)%>%
  
  ungroup()%>%
  group_by(Paire,start_gene,end_gene)%>%
  mutate(size=mean(Size),ratio=mean(Segment_Mean))%>%
  dplyr::select(TYPE,start,end,Num_Probes,Size,size,Segment_Mean,ratio)

Loas_files_annotated1<-Loas_files_annotated%>%
  mutate(n_UD="1")%>%
  group_by(Paire,chr, start,sig) %>%
  mutate(
    n_comm2 = n_distinct(TYPE),
    n_UD = if_else(n_comm2 > 1, paste0(Paire, start), n_UD)
  ) %>%
  ungroup() %>%
  group_by(Paire,chr, end,sig) %>%
  mutate(
    n_comm3 = n_distinct(TYPE),
    n_UD = if_else(n_comm3 > 1, paste0(Paire, end), n_UD)
  ) %>%
  ungroup() %>%
  group_by(Paire, start_gene, end_gene,sig) %>%
  mutate(
    n_comm = n_distinct(TYPE),
    n_UD = if_else(n_comm > 1, paste0(Paire, start_gene, end_gene), n_UD)
  ) %>%
  mutate(n_comm=n_distinct(TYPE))%>%
  
  filter(n_comm==2|n_comm2==2|n_comm3==2)%>%
  ungroup()%>%
  group_by(n_UD)%>%
  mutate(size=mean(Size),ratio=mean(Segment_Mean))%>%
  ungroup()%>%
  dplyr::select(TYPE,Paire,Num_Probes,size,ratio)


full_rect3<-Loas_files_annotated1%>%
  ungroup()%>%
  mutate(sig = ifelse(ratio > 0, "Amp", "Del")) %>%
  dplyr::select(TYPE, Paire, sig, size,ratio)
full_rect3$TYPE<-"Commun"
colnames(full_rect3)<-colnames(Loas_files_epure2)
full_rect3<-unique(full_rect3)

full_rect4<-rbind(full_rect3,Loas_files_epure2)
full_rect5 <- merge(datas_cliniques, full_rect4, by.x = "Paire", by.y = "Paire")

volcano<-ggplot(data = full_rect5,aes(y = full_rect5$size, x = full_rect5$MEAN_LOG2_COPY_RATIO,colour = Type))+
  geom_point(size=3,alpha=1,fill="black")+
  geom_vline(xintercept =log2(3/2))+
  geom_vline(xintercept =log2(1/2))+
  theme_minimal()+
  scale_y_log10() + 
  scale_color_manual(
    values = c(
      "ADC" = "#FF6666",
      "SCC" = "#3399ff",
      "Commun" = "gray"
    )
  )+
  theme(
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 14),   # remove y-axis text
    axis.ticks.y = element_blank(),  # remove y-axis ticks
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    legend.position = "top"
  ) +
  labs(x = "Log2 copy ratio",y = "Size pb") 


ggsave("volcano_cnvs.pdf", volcano, width = 6, height = 6, dpi = 300)


Loas_files_annotated3.5_filled5<-Loas_files_annotated3.5_filled%>%
  group_by(New,TYPE)%>%
  mutate(Total_number=n())%>%
  group_by(New,TYPE,sig)%>%
  filter(n_comm=="AMP")%>%
  mutate(Amp_number=n())%>%
  ungroup()


table_dats<-Ful_datas%>%
  select(New,TYPE,sig,chr,start,end,start_gene,end_gene,Segment_Mean,Num_Probes,genes)

################################

Ful_datas<-read.delim(file = "combined_files.txt", header = T, sep = "\t")
Ful_datas$Paire_num<-gsub("Paire_(\\d+)","\\1",Ful_datas$PAIR)
Ful_datas$Paire<-as.numeric(Ful_datas$Paire_num)
Ful_datas<-Ful_datas[!Ful_datas$Num_Probes<5,]
Ful_datas$Size<-Ful_datas$end-Ful_datas$start
Ful_datas$TYPE[Ful_datas$TYPE=="CIN"]<-"SCC"
Ful_datas$TYPE[Ful_datas$TYPE=="ADC2"]<-"ADC"  
Ful_datas$TYPE[Ful_datas$TYPE=="SCC2"]<-"SCC"
Ful_datas$TYPE[Ful_datas$TYPE=="SCC3"]<-"SCC" 
Ful_datas$sig<-"AMP"
Ful_datas$sig[which(Ful_datas$Segment_Mean<0)]<-"DEL"
Ful_datas$chr_pos<-gsub("chr(.+)","\\1",Ful_datas$chr)
Ful_datas$chr_pos[which(Ful_datas$chr_pos=="X")]<-23
datas_cliniques <- read.delim(file = "./Datas_cliniques.txt", header = TRUE, sep = "\t")
Ful_datas<-Ful_datas[which(Ful_datas$Segment_Mean>0.3|Ful_datas$Segment_Mean< -0.3),]
Loas_files_annotated <- merge(datas_cliniques, Ful_datas, by.x = "Paire", by.y = "Paire")

Loas_files_annotated_full<-Loas_files_annotated%>%
  mutate(n_UD="1")%>%
  group_by(Paire,chr, start,sig) %>%
  mutate(
    n_comm2 = n_distinct(TYPE),
    n_UD = if_else(n_comm2 > 1, paste0(Paire, start), n_UD)
  ) %>%
  ungroup() %>%
  group_by(Paire,chr, end,sig) %>%
  mutate(
    n_comm3 = n_distinct(TYPE),
    n_UD = if_else(n_comm3 > 1, paste0(Paire, end), n_UD)
  ) %>%
  ungroup() %>%
  group_by(Paire, start_gene, end_gene,sig) %>%
  mutate(
    n_comm = n_distinct(TYPE),
    n_UD = if_else(n_comm > 1, paste0(Paire, start_gene, end_gene), n_UD)
  ) %>%
  #mutate(TYPE = ifelse(n_comm > 1, "Commun", TYPE)) %>%
  ungroup()%>%
  mutate(Clonal = ifelse(n_comm==2|n_comm2==2|n_comm3==2, "Yes", "No")) %>%
  ungroup()%>%
  dplyr::select(New,TYPE,sig,chr,start,end,start_gene,end_gene,Segment_Mean,Num_Probes,genes,Clonal)


Loas_files_annotated_full$New<-as.numeric(Loas_files_annotated_full$New)
table_dats_part1 <- Loas_files_annotated_full %>%
  filter(Clonal == "Yes") %>%
  mutate(TYPE = "Shared") %>%
  # supprimer les colonnes indésirables
  dplyr::select(-`Segment_Mean`, -`Num_Probes`, -`start`, -`end`) %>%
  # enlever les lignes en double
  distinct() %>%
  group_by(New, TYPE, sig) %>%
  summarise(CNVs = n(), .groups = "drop") %>%
  complete(New = 1:12,
           TYPE = "Shared",
           sig = c("AMP", "DEL"),
           fill = list(CNVs = 0))


table_dats_part2<-Loas_files_annotated_full%>%
  #filter(Clonal=="No")%>%
  group_by(New,TYPE,sig,)%>%
  summarise(CNVs=n())%>%
  ungroup()%>%
  complete(New = c(1:12),
           TYPE = c("ADC","SCC"),
           sig = c("AMP","DEL"),
           fill = list(CNVs = 0)
  )


table_dats<-rbind(table_dats_part1,table_dats_part2)


colnames(table_dats)<-c("Paire","Cancer","Type of variation","CNVs")


table_dats <- table_dats %>%
  arrange(Paire, Cancer)


write.table(x = table_dats,file = "CNVs_summary.csv",quote = F,sep = ",",row.names = F,col.names = T)

colnames(Loas_files_annotated_full)<-c("Paire","Cancer","Type of variation","Chromosome","Start","End","Start gene","End gene","Log2 fold change","Number of probe","number of Genes","Clonal")
write.table(x = Loas_files_annotated_full,file = "Full_CNVs_summary.csv",quote = F,sep = ",",row.names = F,col.names = T)


##################################################

Loas_files<-read.delim(file = "Test_COMBINED.txt")
Loas_files$Type_order<-1
Loas_files$Type[Loas_files$Type=="CIN3"]<-"SCC"
Loas_files$Type[Loas_files$Type=="ADC2"]<-"ADC"  
Loas_files$Type[Loas_files$Type=="SCC2"]<-"SCC"
Loas_files$Type[Loas_files$Type=="SCC3"]<-"SCC" 
Loas_files$Type_order[which(Loas_files$Type=="ADC")]<-2
Loas_files$Type_order[which(Loas_files$Type=="SCC")]<-3
Loas_files<-unique(Loas_files)
Loas_files$ordre<-gsub('chr(.+)','\\1',Loas_files$CONTIG)
Loas_files$Paire_num<-gsub("Paire_(\\d+)","\\1",Loas_files$Paire)
Loas_files$Paire<-as.numeric(Loas_files$Paire_num)
Loas_files <- Loas_files[!(Loas_files$Paire %in% c(2,5, 7, 9,10,15)), ]
Loas_files$Paire.y<-NULL
Loas_files<-unique(Loas_files)
Loas_files<-Loas_files[!Loas_files$NUM_POINTS_COPY_RATIO<5,]

datas_cliniques <- read.delim(file = "datas_cliniques.txt", header = TRUE, sep = "\t")
Loas_files_annotated <- merge(datas_cliniques, Loas_files, by.x = "Paire", by.y = "Paire")
Loas_files_annotated$Sample<-paste(sep = "_",Loas_files_annotated$Paire,Loas_files_annotated$Type)
Loas_files_annotated<-Loas_files_annotated[which(Loas_files_annotated$MEAN_LOG2_COPY_RATIO< -0.03| Loas_files_annotated$MEAN_LOG2_COPY_RATIO>0.03 ),]
Loas_files_formated<-Loas_files_annotated[,c("Sample","START","END","NUM_POINTS_COPY_RATIO","MEAN_LOG2_COPY_RATIO")]
colnames(Loas_files_formated)<-c("Sample", "Start", "End","Num_Probes", "Segment_Mean")


modified.tai.cancer <- taiModified(Loas_files_formated,segmentMean = 0.2)
cinmetrics.cancer <- CINmetrics(Loas_files_formated,genomeSize_fga=37000000,segmentMean_tai=0,segmentMean_cna = 0,
                                                        segmentMean_base_segments = 0,segmentMean_break_points =0, segmentMean_fga = 0 )

cinmetrics.cancer <- inner_join(cinmetrics.cancer, modified.tai.cancer, by = "sample_id")
cinmetrics.cancer$Paire<-gsub("(\\d+)_(\\w+)","\\1",cinmetrics.cancer$sample_id)
cinmetrics.cancer$Cancer<-gsub("(\\d+)_(\\w+)","\\2",cinmetrics.cancer$sample_id)
cinmetrics.cancer_annotated <- merge(datas_cliniques, cinmetrics.cancer, by.x = "Paire", by.y = "Paire")
cinmetrics.cancer_annotated<-cinmetrics.cancer_annotated[!cinmetrics.cancer_annotated$Cancer=="CTR",]
Full_cinmetrics.cancer_annotated<-cinmetrics.cancer_annotated
cinmetrics.cancer_annotated<-Full_cinmetrics.cancer_annotated
cinmetrics.cancer_annotated$counting<-cinmetrics.cancer_annotated$tai
cinmetrics.cancer_annotated_reduced<-cinmetrics.cancer_annotated[,c("New","Cancer","counting")]

cinmetrics.cancer_annotated_reduced <- cinmetrics.cancer_annotated_reduced %>%
  ungroup() %>%                           # drop grouping
  complete(
    New    = 1:12,
    Cancer = c("ADC", "SCC"),
    fill   = list(counting = 0)          # set missing counting values to 0
  )

cinmetrics.cancer_annotated_reduced<- merge(datas_cliniques, cinmetrics.cancer_annotated_reduced, by.x = "New", by.y = "New")
cinmetrics.cancer_annotated<-cinmetrics.cancer_annotated_reduced
cinmetrics.cancer_annotated$Type_order<-1
cinmetrics.cancer_annotated$Type_order[which(cinmetrics.cancer_annotated$Type=="ADC")]<-2
cinmetrics.cancer_annotated$Type_order[which(cinmetrics.cancer_annotated$Type=="SCC")]<-3

TMB_plot1 <- ggplot(cinmetrics.cancer_annotated, aes(x = reorder(Cancer,Type_order), y = counting)) +
  geom_violin(alpha = 1, fill = "white",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  geom_violin(alpha = 0.4, aes(fill = Cancer),col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") + 
  

#geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(shape = 16, size = 3, width = 0.15, alpha = 0.7) +
  #scale_y_log10() +
  scale_fill_brewer(palette = "Set2") +
  
  scale_fill_manual(
    values = c(
      "ADC" = "#FF6666",
      "SCC" = "#3399ff",
      "Commun" = "gray"
    )
  ) + 
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  labs(y = "TAI score", title = "Cancer") +
  
  stat_compare_means(method = "t.test",label="p.format",label.x.npc = "right") 

cinmetrics.cancer_annotated$Integration[cinmetrics.cancer_annotated$Integration=="INTE_VS"]<-"INTE"
cinmetrics.cancer_annotated$Integration[cinmetrics.cancer_annotated$Integration=="VS"]<-"EPI"
unique(cinmetrics.cancer_annotated$New[cinmetrics.cancer_annotated$Cancer=="ADC"])
cinmetrics.cancer_annotated$type[cinmetrics.cancer_annotated$type=="HPV58"]<-"Other HPVs"
cinmetrics.cancer_annotated$type[cinmetrics.cancer_annotated$type=="HPV45"]<-"Other HPVs"


TMB_plot2 <- ggplot(cinmetrics.cancer_annotated, aes(x = type, y = counting)) +
  geom_violin(alpha = 1, fill = "white",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  geom_violin(alpha = 0.4, fill = "gray",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  #geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(aes( col = Cancer),shape = 16, size = 3, width = 0.15, alpha = 0.7) +  # Add points for individual values
 # scale_y_log10() +  # Use log scale for better visualization
  scale_color_manual(
    values = c(
      "ADC" = "#FF6666",
      "SCC" = "#3399ff",
      "Commun" = "gray"
    )
  )+
  theme_minimal() +
  theme(
    axis.text.x = element_text( size = 14,angle = 45,hjust = 1,vjust = 1),
    axis.text.y =  element_text( size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  
  labs(y = "TAI score", title = "HPV types")+
  stat_compare_means(method = "anova",label="p.format",label.x.npc = "right") 

TMB_plot3 <- ggplot(cinmetrics.cancer_annotated, aes(x = FIGO, y = counting)) +
  geom_violin(alpha = 1, fill = "white",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  geom_violin(alpha = 0.4, fill = "gray",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  #geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(aes(col = Cancer),shape = 16, size = 3, width = 0.15, alpha = 0.7) +  # Add points for individual values
 # scale_y_log10() +  # Use log scale for better visualization
  scale_color_manual(
    values = c(
      "ADC" = "#FF6666",
      "SCC" = "#3399ff",
      "Commun" = "gray"
    )
  )+
  theme_minimal() +
  theme(
    axis.text.x = element_text( size = 14,angle = 45,hjust = 1,vjust = 1),
    axis.text.y =  element_text( size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  labs(y = "TAI score", title = "FIGO cancer classification")+
  stat_compare_means(method = "anova" ,label="p.format",label.x.npc = "right") 

TMB_plot4 <- ggplot(cinmetrics.cancer_annotated, aes(x = Integration,y = counting)) +
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
  theme_minimal() +
  theme(
    axis.text.x = element_text( size = 14,angle = 45,hjust = 1,vjust = 1),
    axis.text.y =  element_text( size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  labs(y = "TAI score", title = "Integration Statut")+
  stat_compare_means(method = "t.test",label="p.format",label.x.npc = "right") 

combi_row <- TMB_plot1 + TMB_plot2 + TMB_plot3 + TMB_plot4+
  plot_layout(ncol = 4)

# Save the plot
ggsave("TAI_plot_by_Strain2_FULLs00000_select.pdf", width = 12, height = 6, dpi = 300)


dev.off()



###########################break_points #####################################"

cinmetrics.cancer_annotated<-Full_cinmetrics.cancer_annotated
cinmetrics.cancer_annotated$counting<-cinmetrics.cancer_annotated$break_points
cinmetrics.cancer_annotated_reduced<-cinmetrics.cancer_annotated[,c("New","Cancer","counting")]
cinmetrics.cancer_annotated_reduced <- cinmetrics.cancer_annotated_reduced %>%
  ungroup() %>%                           # drop grouping
  complete(
    New    = 1:12,
    Cancer = c("ADC", "SCC"),
    fill   = list(counting = 0)          # set missing counting values to 0
  )

cinmetrics.cancer_annotated_reduced<- merge(datas_cliniques, cinmetrics.cancer_annotated_reduced, by.x = "New", by.y = "New")
cinmetrics.cancer_annotated<-cinmetrics.cancer_annotated_reduced

cinmetrics.cancer_annotated$Type_order<-1
cinmetrics.cancer_annotated$Type_order[which(cinmetrics.cancer_annotated$Type=="ADC")]<-2
cinmetrics.cancer_annotated$Type_order[which(cinmetrics.cancer_annotated$Type=="SCC")]<-3


TMB_plot1 <- ggplot(cinmetrics.cancer_annotated, aes(x = reorder(Cancer,Type_order), y = counting)) +
  geom_violin(alpha = 1, fill = "white",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  geom_violin(alpha = 0.4, aes(fill = Cancer),col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") + 
  
  
  #geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(shape = 16, size = 3, width = 0.15, alpha = 0.7) +
  #scale_y_log10() +
  scale_fill_brewer(palette = "Set2") +
  
  scale_fill_manual(
    values = c(
      "ADC" = "#FF6666",
      "SCC" = "#3399ff",
      "Commun" = "gray"
    )
  ) + 
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  labs(y = "Break points", title = "Cancer") +
  
  stat_compare_means(method = "t.test",label="p.format",label.x.npc = "right") 


cinmetrics.cancer_annotated$Integration[cinmetrics.cancer_annotated$Integration=="INTE_VS"]<-"INTE"
cinmetrics.cancer_annotated$Integration[cinmetrics.cancer_annotated$Integration=="VS"]<-"EPI"

cinmetrics.cancer_annotated$type[cinmetrics.cancer_annotated$type=="HPV58"]<-"Other HPVs"
cinmetrics.cancer_annotated$type[cinmetrics.cancer_annotated$type=="HPV45"]<-"Other HPVs"


TMB_plot2 <- ggplot(cinmetrics.cancer_annotated, aes(x = type, y = counting)) +
  geom_violin(alpha = 1, fill = "white",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  geom_violin(alpha = 0.4, fill = "gray",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  #geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(aes( col = Cancer),shape = 16, size = 3, width = 0.15, alpha = 0.7) +  # Add points for individual values
  # scale_y_log10() +  # Use log scale for better visualization
  scale_color_manual(
    values = c(
      "ADC" = "#FF6666",
      "SCC" = "#3399ff",
      "Commun" = "gray"
    )
  )+
  theme_minimal() +
  theme(
    axis.text.x = element_text( size = 14,angle = 45,hjust = 1,vjust = 1),
    axis.text.y =  element_text( size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  
  labs(y = "Break points", title = "HPV types")+
  stat_compare_means(method = "anova",label="p.format",label.x.npc = "right") 
TMB_plot2

TMB_plot3 <- ggplot(cinmetrics.cancer_annotated, aes(x = FIGO, y = counting)) +
  geom_violin(alpha = 1, fill = "white",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  geom_violin(alpha = 0.4, fill = "gray",col="white",draw_quantiles = 0.5,quantile.linewidth = 1.5,quantile.color = "black") +
  #geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
  geom_jitter(aes(col = Cancer),shape = 16, size = 3, width = 0.15, alpha = 0.7) +  # Add points for individual values
  # scale_y_log10() +  # Use log scale for better visualization
  scale_color_manual(
    values = c(
      "ADC" = "#FF6666",
      "SCC" = "#3399ff",
      "Commun" = "gray"
    )
  )+
  theme_minimal() +
  theme(
    axis.text.x = element_text( size = 14,angle = 45,hjust = 1,vjust = 1),
    axis.text.y =  element_text( size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  labs(y = "Break points", title = "FIGO cancer classification")+
  stat_compare_means(method = "anova" ,label="p.format",label.x.npc = "right") 



TMB_plot4 <- ggplot(cinmetrics.cancer_annotated, aes(x = Integration,y = counting)) +
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
  theme_minimal() +
  theme(
    axis.text.x = element_text( size = 14,angle = 45,hjust = 1,vjust = 1),
    axis.text.y =  element_text( size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  labs(y = "Break points", title = "Integration Statut")+
  stat_compare_means(method = "t.test",label="p.format",label.x.npc = "right") 

combi_row <- TMB_plot1 + TMB_plot2 + TMB_plot3 + TMB_plot4+
  plot_layout(ncol = 4)

# Save the plot
ggsave("Break_points_plot_by_Strain2_FULLs00000_select.pdf", width = 12, height = 6, dpi = 300)

dev.off()




