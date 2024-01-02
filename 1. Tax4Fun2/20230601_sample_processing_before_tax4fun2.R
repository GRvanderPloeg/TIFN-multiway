wd <- "."
setwd(wd)

library(vegan)
library(Biostrings)
library(tidyverse)

counts <- read.delim("../Data/20221005_wp2/count-table.tsv")
taxonomy <- read.delim("../Data/20221005_wp2/taxonomic-classification.tsv")
rare_threshold = 40000 # also consider 15k,  25k, 40k

counts <- counts[,!colnames(counts) %in% tax$asv[which(taxonomy$Order=="Chloroplast" | taxonomy$Family == "Mitochondria")]]
taxonomy <- taxonomy[-which(taxonomy$Order=="Chloroplast" | taxonomy$Family == "Mitochondria"),]
otus <- DNAStringSet(taxonomy$representative_sequence)
names(otus) <- taxonomy$asv

# Remove test individuals
counts_control = counts[counts$group == "control",]

# Diagnostics
gm <- apply(counts_control[,6:ncol(counts_control)],1,function(x) exp(mean(log(x[x>0]))))
gm2 <- apply(counts_control[,6:ncol(counts_control)],1,function(x) exp(mean(log(c(x[x>0],rep(0.5,length(which(x==0))))))))

hist(rowSums(counts_control[,6:ncol(counts_control)]))
summary(rowSums(counts_control[,6:ncol(counts_control)]))

# Boxplot of all sequencing depths
cbind(counts_control$niche, rowSums(counts_control[,6:ncol(counts_control)])) %>%
  as_tibble() %>%
  mutate(V1 = as.factor(V1), V2 = as.numeric(V2)) %>%
  ggplot(aes(x=V1,y=V2)) + 
  geom_boxplot(width=0.5) + 
  xlab("Sampling location") + 
  ylab("Sample depth") + 
  geom_hline(yintercept=rare_threshold, col="red") +
  ggtitle("Sequencing depth per sampling site")

# Rarification curves
#rarecurve(counts_control[,6:ncol(counts_control)], step=1000, label=FALSE)

# IF we rarefy to n reads, what do we remove?
tab10 <- counts_control[which(rowSums(counts_control[,6:ncol(counts_control)])>rare_threshold),]
sn10 <- table(tab10$subject,tab10$niche)
rm10 <- cbind(rep(levels(as.factor(tab10$subject)),times=length(levels(as.factor(tab10$niche)))),
              rep(levels(as.factor(tab10$niche)),each=length(levels(as.factor(tab10$subject)))))[which(!sn10 %in% c(0,7)),]

counts_filtered <- tab10[apply(tab10[,c("subject","niche")],1,function(x) length(which(rm10[,1]==x[1]&rm10[,2]==x[2]))==0),]
counts_filtered <- counts_filtered[,c(1:5,5+which(colSums(counts_filtered[,6:ncol(counts_filtered)])>0))]
taxonomy_filtered <- taxonomy[which(taxonomy$asv %in% colnames(tabu10)),]

# Rarefication
set.seed(607)
counts_rarefied <- data.frame(counts_filtered[,1:5],rrarefy(counts_filtered[,6:ncol(counts_filtered)],rare_threshold))
counts_rarefied <- counts_rarefied[,c(1:5,5+which(colSums(counts_rarefied[,6:ncol(counts_rarefied)])>0))]
taxonomy_rarefied <- taxonomy_filtered[which(taxonomy_filtered$asv %in% colnames(counts_rarefied)),]
otu_rarefied <- DNAStringSet(taxonomy_rarefied$representative_sequence)
names(otu_rarefied) <- taxonomy_rarefied$asv
counts_rarefied_numeric <- data.frame(t(counts_rarefied[,6:ncol(counts_rarefied)]))
colnames(counts_rarefied_numeric) <- counts_rarefied_numeric$sample
counts_rarefied_numeric$OTU <- rownames(counts_rarefied_numeric)

# Diagnostics
gm.10 <- apply(tabr10[,6:ncol(tabr10)],1,function(x) exp(mean(log(c(x[x>0])))))
gm2.10 <- apply(tabr10[,6:ncol(tabr10)],1,function(x) exp(mean(log(c(x[x>0],rep(0.5,length(which(x==0))))))))

saveRDS(counts_rarefied_numeric,"TIFN_Tax4Fun2_rarefied_40k.RDS")
writeXStringSet(otu_rarefied,"TIFN_Tax4Fun2_rarefied_40k_OTUs.fa")
