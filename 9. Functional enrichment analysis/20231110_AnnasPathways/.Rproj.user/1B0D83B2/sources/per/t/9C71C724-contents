setwd("~/Documents/UvA/Emission/Roel1/20231102_setRankInput")

pwn <- read.delim("pwn.txt",
                 header=F)
colnames(pwn) <- c("pw","name")
pwn$pw <- gsub("map","ko",pwn$pw)

kon <- read.delim("kon.txt",
                  header=F,quote="")
colnames(kon) <- c("koID","description")

pw2ko <- read.delim(file.path("ko2pw.txt"),header=F) 
#pw2ko <- pw2ko[c(grep("path:map00",pw2ko[,1]),grep("path:map01",pw2ko[,1]),
#                 grep("path:ko00",pw2ko[,1]),grep("path:ko01",pw2ko[,1])),]
pw2ko[,1] <- gsub("path:","",pw2ko[,1])
pw2ko[,2] <- gsub("ko:","",pw2ko[,2])
colnames(pw2ko) <- c("pw","ko")
pw2ko$pw <- gsub("map","ko",pw2ko$pw)


pw2cpd <- read.delim(file.path("cmp2pw.txt"),header=F) 
#pw2cpd <- pw2cpd[c(grep("path:map00",pw2cpd[,1]),grep("path:map01",pw2cpd[,1])),]
pw2cpd[,1] <- gsub("path:","",pw2cpd[,1])
pw2cpd[,2] <- gsub("cpd:","",pw2cpd[,2])
colnames(pw2cpd) <- c("pw","cpd")
pw2cpd$pw <- gsub("map","ko",pw2cpd$pw)

con <- read.delim("cmpn.txt",
                  header=F)
colnames(con) <- c("cpd","name")

ll <- read.delim("../20231102_setRankInput/lowlingSetRankInput.csv",sep=",")
ul <- read.delim("../20231102_setRankInput/uplingSetRankInput.csv",sep=",")
ui <- read.delim("../20231102_setRankInput/upinterSetRankInput.csv",sep=",")
sa <- read.delim("../20231102_setRankInput/salivaSetRankInput.csv",sep=",")
to <- read.delim("../20231102_setRankInput/tongueSetRankInput.csv",sep=",")

gs <- list(ll,ul,ui,sa,to)
names(gs) <- c("ll","ul","ui","sa","to")

pws <- c("ko01051","ko01240","ko01250","ko01200","ko00562","ko00540",
         "ko00670","ko00190","ko00550","ko00230","ko00240","ko02024",
         "ko03018","ko00900")

# KOs per pathway
         
sapply(1:length(pws), function(x){
  sapply(1:length(gs), function(y){
    length(unique(pw2ko$ko[pw2ko$pw==pws[x] &
                             pw2ko$ko %in% gs[[y]]$geneID]))
  })
})

# compounds per pathway

sapply(1:length(pws), function(x){
  sapply(1:length(gs), function(y){
    length(unique(pw2cpd$cpd[pw2cpd$pw==pws[x] &
                             pw2cpd$cpd %in% gs[[y]]$geneID]))
  })
})

# both per pathway

sapply(1:length(pws), function(x){
  sapply(1:length(gs), function(y){
    length(unique(c(pw2cpd$cpd[pw2cpd$pw==pws[x] &
                               pw2cpd$cpd %in% gs[[y]]$geneID],
                    pw2ko$ko[pw2ko$pw==pws[x] &
                               pw2ko$ko %in% gs[[y]]$geneID])))
  })
})

#       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
# [1,]    1  102   24   60    5    6   11   18   17    52    43    31     4     8
# [2,]    1   84   25   59    5    5   11   14   15    51    47    37     4     4
# [3,]    1   96   25   59    4   11   11   16   16    53    45    30     4     8
# [4,]    2   90   28   58    4   15   11   16   17    51    47    30     4     9
# [5,]    2   88   29   54    4   17   11   16   17    51    46    32     4     9

pathCol_KO <- lapply(gs,function(x){
  out <- x$score[grep("^K",x$geneID)]
  names(out) <- x$geneID[grep("^K",x$geneID)]
  out
})    
pathCol_cmp <- lapply(gs,function(x){
  out <- x$score[grep("^C",x$geneID)]
  names(out) <- x$geneID[grep("^C",x$geneID)]
  out
})  

for(cpw in pws[c(1,5,7:14)]){
  for(i in 1:length(gs)){
    pv <- pathview(pathCol_KO[[i]],
                   cpd.data = pathCol_cmp[[i]],
                   pathway.id=cpw,species="ko",
                   gene.idtype="KEGG",
                   both.dirs = list(gene = F, cpd = F),
                   out.suffix=names(gs)[i],
                   low = list(gene = "#2C7BB6", cpd = "#2C7BB6"), 
                   mid =list(gene = "#FFFFBF", cpd = "#FFFFBF"), 
                   high = list(gene = "#D7191C", cpd ="#D7191C"),
                   limit=list(gene=c(0.7,1), cpd=c(0.7,1)))
}}

# 1: too few genes/no compounds  
# 2/3/4/6: cofactors, nucleotide sugars, Cmetabolism, LPS don't work
# 5: Inositol phosphate, just 4-5 genes/compounds, myo-inositol has a strong signal, but enzymes aren't directly connected
# 7: 1-C by folate: several KOs, no compounds
# 8: ox phos: mainly ATPase
# 9: peptidoglycan: strong, most enzymes (no compounds)
# 10: purine/11: pyrimidine : urea strong in compounds, some of the other compounds and enzymes
# 12: quorum sensing: OPP (oligopeptide binding protein - Streptococcus)
# 13: RNA degradation: enolase, just 4 genes
# 14: terpenoids: non-mevalonate pathway

# saliva: co-factors, C metabolism, ox phos, peptidoglycan, purine (without tongue), 
#         quorum sensing (without tongue)

# plaques: nucleotide sugars - 1C-folate - ox phos - peptidoglycan - QS (all), 
#          co-factors & C metabolism, isoprenoids (ui), 

# LPS: tongue
