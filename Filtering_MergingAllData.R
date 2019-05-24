### This script filters PAR-CLIP data based on RNA-seq expression > 1 FPKM and then catenates RIP, PAR and RNA-seq data for further analysis 

library(tidyverse)
library(dplyr)
library(reshape2)
library(ggplot2)
library(data.table)
library(mltools)

###Import PARCLIP cluster data per gene. File contains number of binding sites for each transcript

naive_parclipdata<-read.csv("~/RNASeq_Analysis/PARCLIP/ELAVL1_Strict_2018/minuscGAMPTHP1/untrimmed1.gene_cl.csv") %>% 
  select(GeneName, X5.utr,Intron,Exon, X3.utr)
IRF3_parclipdata<-read.csv("~/RNASeq_Analysis/PARCLIP/ELAVL1_Strict_2018/pluscGAMPTHP1/untrimmed2.gene_cl.csv") %>% 
  select(GeneName, X5.utr,Intron,Exon, X3.utr)

##renaming gene id columns
colnames(IRF3_parclipdata)[1] <- "gene_short_name"
colnames(naive_parclipdata)[1] <- "gene_short_name"

#### Importing Expression data
expressiondata<-read_tsv("/Volumes/RAID1/RNASeq_Analysis/KR_THP1/EncapsulatedProject/CuffDiffOutput/deNOVO_run/cds.fpkm_tracking") %>%
  select(gene_short_name,THP1_0H_1_1_FPKM, THP1_16H_1_1_FPKM)

### Remove duplicate gene names
#expressiondata2<-expressiondata[!duplicated(expressiondata[c('gene_short_name')]),]

expressiondata_max <- expressiondata %>%
  group_by(gene_short_name) %>%
  summarize(THP1_0hr = max(THP1_0H_1_1_FPKM),
            THP1_16hr = max(THP1_16H_1_1_FPKM))


### Remove lowly expressing genes. Genes must have a FPKM > 2 in one of the conditions
expressiondata3<-expressiondata_max %>% filter(expressiondata_max$THP1_0hr > 2 | expressiondata_max$THP1_16hr > 2)

#### Alternative options: Remove FPKM < 2 for each condition seperately 
#expressiondata_naive<-expressiondata_max %>% filter(expressiondata_max$THP1_0hr > 2) 
#expressiondata_IRF3<-expressiondata_max %>% filter(expressiondata_max$THP1_16hr > 2)

### log2 transform and then change all created infinities to zero 
expressiondata3[, 2:3] <- log2(expressiondata3[,2:3])
#expressiondata_naive[, 2:3] <- log2(expressiondata_naive[,2:3])
#expressiondata_IRF3[, 2:3] <- log2(expressiondata_IRF3[,2:3])

### chaning all infinites to zero
expressiondata4<-expressiondata3 %>% mutate_at(c(2:3), funs(replace(., is.infinite(.), 1))) 
#expressiondata_naive_2<-expressiondata_naive %>% mutate_at(c(2:3), funs(replace(., is.infinite(.), 1))) 
#expressiondata_IRF3_2<-expressiondata_IRF3 %>% mutate_at(c(2:3), funs(replace(., is.infinite(.), 1))) 

#### joining PARCLIP and Expression Data
naive_dat<-full_join(naive_parclipdata, expressiondata4[,1:2], by="gene_short_name")
IRF3_dat<-full_join(IRF3_parclipdata, expressiondata4[,c(1,3)], by="gene_short_name")

##Testing
#only_expression <- expressiondata4 %>% filter(!(gene_short_name %in% naive_parclipdata$gene_short_name))

#### Replacing only in the PARCLIP data NAs with Zeros
naive_data<- naive_dat %>% mutate_at(c(2:5), funs(replace(., is.na(.), 0))) 
naive_data2<- na.omit(naive_data)

IRF3_data<- IRF3_dat %>% mutate_at(c(2:5), funs(replace(., is.na(.), 0)))
IRF3_data2<- na.omit(IRF3_data)

###Importing RIP-seq data
RIP_data<-read_tsv("/Volumes/RAID1/RNASeq_Analysis/RIPseq_Analysis/CuffDiffOutput/cds.fpkm_tracking") %>%
  select("gene_short_name", ends_with( "_FPKM"))

RIP_data2<-RIP_data[!duplicated(RIP_data[c('gene_short_name')]),]

### Subtracting IgG control from FLAG in RIPseq data. Also log2 transforming RIPseq data

RIP_data2[,"q1corrected"]<-log2(RIP_data2$q1_FPKM-RIP_data2$q5_FPKM)
RIP_data2[,"q2corrected"]<-log2(RIP_data2$q2_FPKM-RIP_data2$q6_FPKM)
RIP_data2[,"q3corrected"]<-log2(RIP_data2$q3_FPKM-RIP_data2$q8_FPKM)
RIP_data2[,"q4corrected"]<-log2(RIP_data2$q4_FPKM-RIP_data2$q9_FPKM)

RIP_data2<-RIP_data2%>% mutate_at(c(12:15), funs(replace(., is.na(.), 1)))
RIP_data2<-RIP_data2%>% mutate_at(c(12:15), funs(replace(., is.infinite(.), 1)))

###Joining RIP to PAR_Expression data
full_naive<-full_join(naive_data2,RIP_data2[,c(1,12)], by ="gene_short_name")
full_naive2<- na.omit(full_naive)
full_naive_mRNA<- full_naive2 %>% filter (X5.utr > 0 | Intron > 0 | X3.utr >0 | Exon > 0)
##write_csv(full_naive2, "naiveTHP1_PAR_RIP_Exp_data")

full_IRF3<-full_join(IRF3_data2, RIP_data2[,c(1,13)], by="gene_short_name")
full_IRF3_2<-na.omit(full_IRF3)
full_IRF3_mRNA<- full_IRF3_2 %>% filter (X5.utr > 0 | Intron > 0 | X3.utr >0 | Exon > 0) 
##write_csv(full_IRF3_2, "IRF3THP1_PAR_RIP_Exp_data")

x<-(setdiff(full_naive_mRNA$gene_short_name, full_IRF3_mRNA$gene_short_name))
str(x)
### Importing specific Clustering Information 
clusterdataNaive<-read.csv("~/RNASeq_Analysis/PARCLIP/ELAVL1_Strict_2018/minuscGAMPTHP1/untrimmed1.clusters.csv") %>%
  select('GeneName', 'Chr' , 'Start','End', 'ClusterSequence', 'Aligned.to', 'T2Cfraction','Strand')

clusterdataIRF3<-read.csv("~/RNASeq_Analysis/PARCLIP/ELAVL1_Strict_2018/pluscGAMPTHP1/untrimmed2.clusters.csv") %>%
  select('GeneName', 'Chr' , 'Start','End', 'ClusterSequence', 'Aligned.to', 'T2Cfraction','Strand')

colnames(clusterdataNaive)[1] <- "gene_short_name"
colnames(clusterdataIRF3)[1] <- "gene_short_name"


full_cluster_naive<-left_join(full_naive2,clusterdataNaive, by = "gene_short_name")
write.csv(full_cluster_naive, "full_cluster_naive")
full_cluster_IRF3<-left_join(full_IRF3_2,clusterdataIRF3, by = "gene_short_name")
write.csv(full_cluster_IRF3, "full_cluster_IRF3")

#### Make a bedfile of clusters that have been filtered: FPKM > 2 in one of the two conditions 
bed_naive_clusters<-full_cluster_naive %>% select('Chr', 'Start','End','Strand', 'ClusterSequence', 'Aligned.to')
write.csv(bed_naive_clusters, "bed_naive_clusters")
bed_IRF3_clusters<-full_cluster_IRF3 %>% select('Chr', 'Start','End','Strand', 'ClusterSequence', 'Aligned.to')
write.csv(bed_IRF3_clusters, "bed_IRF3_clusters")

length(which(full_IRF3_mRNA$q2corrected > 0))


minus_THP1_specificClusters<-read_tsv("~/RNASeq_Analysis/PARCLIP/ELAVL1_Strict_2018/BEDFILES/naive_specific", col_names = FALSE)
plus_THP1_specificClusters<-read_tsv("~/RNASeq_Analysis/PARCLIP/ELAVL1_Strict_2018/BEDFILES/IRF3_specific", col_names = FALSE)

colnames(minus_THP1_specificClusters)<-c("Chr", "Start", "End", "Strand")
colnames(plus_THP1_specificClusters)<-c("Chr", "Start", "End", "Strand")


specifc_naive_clusters<-left_join(minus_THP1_specificClusters, full_cluster_naive, by= "Start")
write.csv(specifc_IRF3_clusters, "specifc_IRF3_clusters")
specifc_IRF3_clusters<-left_join(plus_THP1_specificClusters, full_cluster_IRF3, by= "Start")


summary(specifc_IRF3_clusters)
summary(specifc_naive_clusters)

lines(density(specifc_IRF3_clusters$THP1_16hr))
plot(density(specifc_naive_clusters$THP1_0hr), col="red")
