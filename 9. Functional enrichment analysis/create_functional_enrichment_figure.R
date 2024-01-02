library(tidyverse)
library(ggpubr)

source("../R scripts/summary.functions.R")
source("../R scripts/PARAFAC_functions.R")
set.seed(0)

#tongueCombResult = readRDS("./tongue_pathway_enrichment_result.RDS")
#lowlingCombResult = readRDS("./lowling_pathway_enrichment_result.RDS")
#lowinterCombResult = readRDS("./lowinter_pathway_enrichment_result.RDS")
#uplingCombResult = readRDS("./upling_pathway_enrichment_result.RDS")
#upinterCombResult = readRDS("./upinter_pathway_enrichment_result.RDS")
#salivaCombResult = readRDS("./saliva_pathway_enrichment_result.RDS")

#allResults = list(tongueCombResult, lowlingCombResult, lowinterCombResult, uplingCombResult, upinterCombResult, salivaCombResult)
nicheNames = c("tongue", "lower jaw, lingual", "lower jaw, interproximal", "upper jaw, lingual", "upper jaw, interproximal", "saliva")

#pathwayResults = ""
#maxPvalue = 0.1
#for(i in 1:length(allResults)){
#  result = allResults[[i]]
#  pathwayResults = c(pathwayResults, result %>% filter(rowMean < maxPvalue | rowMedian < maxPvalue) %>% select(pathway) %>% pull)
#}
#pathwayResults = pathwayResults[2:length(pathwayResults)]

#pathwaysOfInterest = pathwayResults %>% table() %>% as_tibble() %>% filter(n >= 1) %>% select(".") %>% pull

#plottableData = tongueCombResult %>% filter(pathway %in% pathwaysOfInterest) %>% select(pathway, rowMedian, rowStd) %>% mutate(niche = "tongue", logValue = -log10(rowMedian))

#for(i in 2:length(allResults)){
#  result = allResults[[i]]
#  name = nicheNames[i]
#  addition = result %>% filter(pathway %in% pathwaysOfInterest) %>% select(pathway, rowMedian, rowStd) %>% mutate(niche = name, logValue = -log10(rowMedian))
#  plottableData = rbind(plottableData, addition) %>% as_tibble()
#}

#plottableData %>% ggplot(aes(x=as.factor(pathway), y=logValue, fill=as.factor(niche))) +
#  geom_bar(stat="identity", position=position_dodge()) +
#  geom_hline(yintercept=(-log10(0.05)), col="red") +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#  xlab("KEGG Pathway") +
#  ylab("-log10(p-value)")  


# Alternative based on setRank

tongueSetRankResult = read.csv("./20230704_run_onlyVarExp/Tongue_SetRank_result_pathways.txt", sep="\t")
lowlingSetRankResult = read.csv("./20230704_run_onlyVarExp/Lowling_SetRank_result_pathways.txt", sep="\t")
lowinterSetRankResult = read.csv("./20230704_run_onlyVarExp/Lowinter_SetRank_result_pathways.txt", sep="\t")
uplingSetRankResult = read.csv("./20230704_run_onlyVarExp/Upling_SetRank_result_pathways.txt", sep="\t")
upinterSetRankResult = read.csv("./20230704_run_onlyVarExp/Upinter_SetRank_result_pathways.txt", sep="\t")
salivaSetRankResult = read.csv("./20230704_run_onlyVarExp/Saliva_SetRank_result_pathways.txt", sep="\t")

allResults2 = list(tongueSetRankResult, lowlingSetRankResult, lowinterSetRankResult, uplingSetRankResult, upinterSetRankResult, salivaSetRankResult)
minPvalue = 0.005
pathwayResults2 = ""

for(i in 1:length(allResults2)){
  result = allResults2[[i]]
  pathwayResults2 = c(pathwayResults2, result %>% filter(correctedPValue < minPvalue) %>% select(description) %>% pull)
}

minNiches = 1
selectedPathways = pathwayResults2 %>% table() %>% as_tibble() %>% filter(n >= minNiches)
selectedPathways = selectedPathways[,1] %>% pull()

plottableData2 = tongueSetRankResult %>% filter(description %in% selectedPathways) %>% select(description, correctedPValue) %>% mutate(niche = "tongue", logValue = -log10(correctedPValue))

for(i in 2:length(allResults2)){
  result = allResults2[[i]]
  name = nicheNames[i]
  addition = result %>% filter(description %in% selectedPathways) %>% select(description, correctedPValue) %>% mutate(niche = name, logValue = -log10(correctedPValue))
  plottableData2 = rbind(plottableData2, addition) %>% as_tibble()
}

plottableData2 %>% ggplot(aes(x=as.factor(description), y=logValue, fill=as.factor(niche))) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") + 
  geom_hline(yintercept=(-log10(0.05)), col="red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("KEGG Pathway") +
  ylab("-log10(p-value)") +
  scale_fill_discrete(name = "Sampling location") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# test figures - load models
microb_featureNames = c("KO", "description")
metab_featureNames = c("SUPER_PATHWAY","SUB_PATHWAY","BIOCHEMICAL")
microb_days = c(-14,0,2,5,9,14,21)
metab_days = c(0,2,5,9,14)

path = "../7. Functional microbiome modeling/20230607_run/PARAFAC models/"

tongueFuncModel = importPARAFAC(paste0(path,"Tongue_func"), microb_featureNames, microb_days, "KO")
lowlingFuncModel = importPARAFAC(paste0(path,"Low_ling_func"), microb_featureNames, microb_days, "KO")
lowinterFuncModel = importPARAFAC(paste0(path,"Low_inter_func"), microb_featureNames, microb_days, "KO")
uplingFuncModel = importPARAFAC(paste0(path,"Up_ling_func"), microb_featureNames, microb_days, "KO")
upinterFuncModel = importPARAFAC(paste0(path,"Up_inter_func"), microb_featureNames, microb_days, "KO")
salivaFuncModel = importPARAFAC(paste0(path, "Saliva_func"), microb_featureNames, microb_days, "KO")

metabolomicsModel = importPARAFAC("../6. Metabolome modeling/20230605_run/PARAFAC models/Metabolomics", metab_featureNames, metab_days, "BIOCHEMICAL")

rf_data = read.csv("../0. Raw data input/RFdata.csv")
colnames(rf_data) = c("subject", "id", "fotonr", "day", "group", "RFgroup", "MQH", "SPS(tm)", "Area_delta_R30", "Area_delta_Rmax", "Area_delta_R30_x_Rmax", "gingiva_mean_R_over_G", "gingiva_mean_R_over_G_upper_jaw", "gingiva_mean_R_over_G_lower_jaw")
rf_data = rf_data %>% as_tibble()

rf_data[rf_data$subject == "VSTPHZ", 1] = "VSTPH2"
rf_data[rf_data$subject == "D2VZH0", 1] = "DZVZH0"
rf_data[rf_data$subject == "DLODNN", 1] = "DLODDN"
rf_data[rf_data$subject == "O3VQFX", 1] = "O3VQFQ"
rf_data[rf_data$subject == "F80LGT", 1] = "F80LGF"
rf_data[rf_data$subject == "26QQR0", 1] = "26QQrO"

rf_data2 = read.csv("../0. Raw data input/red_fluorescence_data.csv") %>% as_tibble()
rf_data2 = rf_data2[,c(2,4,181:192)]
rf_data = rf_data %>% left_join(rf_data2)

rf = rf_data %>% select(subject, RFgroup) %>% unique()

# Mapping of KO K-number to pathways
KO2pathway = read.csv("../4. Functional microbiome pre-processing/kegg_ko2pathway.txt", sep="\t", header=FALSE) %>% as_tibble()
colnames(KO2pathway) = c("pathway", "KO")
dummy = unlist(strsplit(KO2pathway["KO"] %>% pull, "ko:"))
KO2pathway["KO"] = dummy[dummy != ""]

# Mapping of pathway ID to pathway name
pathway2name = read.csv("../4. Functional microbiome pre-processing/kegg_pathwayNames.txt", sep="\t", header=FALSE) %>% as_tibble()
colnames(pathway2name) = c("pathway", "name")
#pathway2name2 = pathway2name
#pathway2name2$pathway = str_replace(pathway2name$pathway, "map", "ko")

# Mapping of compound to pathway
compound2pathway = read.csv("../4. Functional microbiome pre-processing/kegg_compound2pathway.txt", sep="\t", header=FALSE) %>% as_tibble()
colnames(compound2pathway) = c("pathway", "compound")
dummy = unlist(strsplit(compound2pathway["compound"] %>% pull, "cpd:"))
compound2pathway["compound"] = dummy[dummy != ""]

# test figures, mark important pathways

