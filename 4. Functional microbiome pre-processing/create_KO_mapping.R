library(tidyverse)
library(stringr)

# Load input data
metabolomics_metadata = read.csv("../0. Raw data input/20221005_wp2/metabolomics-features.tsv", sep="\t")
compound2pathway = read.csv("./kegg_compound2pathway.txt", header=FALSE, sep="\t")
colnames(compound2pathway) = c("pathway", "KEGG")
ko2pathway = read.csv("./kegg_ko2pathway.txt", header=FALSE, sep="\t")
colnames(ko2pathway) = c("pathway", "KO")

# Map metabolomics compounds to KEGG pathways
mappedPathways = unique(compound2pathway[paste0("cpd:",metabolomics_metadata$KEGG) %in% compound2pathway$KEGG,]$pathway)

# Select KOs mapping to only these pathways
KOselection = ko2pathway[ko2pathway$pathway %in% mappedPathways,]$KO
KOselection = str_split_fixed(KOselection, "ko:", 2)[,2]

# Create output file
write.table(KOselection, "./KOselection.csv", row.names=FALSE, col.names=FALSE)
