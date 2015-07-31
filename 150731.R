#install.packages()
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("GenomicFeatures","AnnotationDbi"))
#prompts to update all some none
biocLite(c("rgl","biomaRt","Rcade"))
library("Rcade","rgl","biomaRt")
dir<-file.path(system.file("extdata", package="Rcade"), "STAT1")
dir
DE<-read.csv(file.path(dir, "DE.csv"))
DE
DElookup<- list(GeneID="ENSG", logFC="logFC", B="B", "Genes.Location","Symbol")
DElookup
dir(dir,pattern=".bam")
targets<-read.csv(file.path(dir,"targets.csv"), as.is=TRUE)
anno<-read.csv(file.path(dir,"anno.csv"))
anno<-anno[order(anno$chromosome_name),]
anno
colnames(anno) <- c("ENSG","chr","start","end","str")
library("biomaRt")
library("GenomicRanges")

ChIPannoZones <- defineBins(anno, zone=c(-1500,1500), geneID="ENSG")
#depends on libraries

#Rcade's prior belief is that each gene's DE and ChIP seq status are independent.
#unlike to be true in real data
DE.prior =0.01
prior.mode="keepChIP"
prior=c("D|C" = 0.05, "D|notC" = 0.005)

#analyse
Rcade<- RcadeAnalysis(DE,ChIPannoZones, annoZoneGeneidName = "ENSG",
                      ChIPtargets=targets, ChIPfileDir=dir, DE.prior=DE.prior,
                      prior.mode=prior.mode, prior=prior, DElookup=DElookup)
Rcade

xDE<-getDE(Rcade)
xChIP<-getChIP(Rcade)
xRcade<-getRcade(Rcade)

plotPCA(Rcade)
plotMM(Rcade)
library(rgl)
plotBBB(Rcade)

exportRcade(Rcade, directory="RcadeOutput", cutoffArg=2000)
exportRcade(Rcade, directory="RcadeOutput", cutoffMode="top", cutoffArg=1000,
            justGeneID=FALSE, removeDuplicates="beforeCutoff")
