# Required libraries
library(plyr)
library(edgeR)
library(goseq)
library(GO.db)
library(GenomicRanges)
library(lattice)
library(RColorBrewer)
library(cluster)
library(xtable)
library(limma)

library(edgeR)
table.summary<-read.table(system.file("extdata","Li_sum.txt", package="goseq"), 
                          sep="\t", header=TRUE, stringsAsFactors=FALSE)
counts <- table.summary[,-1]
counts
rownames(counts)=table.summary[,1]
rownames(counts)
grp=factor(rep(c("Control", "Treated"), times=c(4,3)))
grp
summarized=DGEList(counts,lib.size=colSums(counts),group=grp)

#user edgeR to estimate the biological dispersion and 
#calculate differtial expression using a negative binomial model
disp=estimateCommonDisp(summarized)
disp$common.dispersion
disp
str(disp)
tested=exactTest(disp)
topTags(tested)
#comparison treated - control

#Format into a vector
?p.adjust
genes=as.integer(p.adjust(tested$table$PValue[tested$table$logFC!=0],
                          method="BH")<.05)
# caution! .05 not 0.5
#BH aka fdr - implements Benjamini Hochberg method | BY - B, Yekutieli
head(genes)
names(genes)=row.names(tested$table[tested$table$logFC!=0,])
names(genes)
table(genes)

supportedGenomes()
#error - need library
library(goseq)
#supportedGenomes()
head(supportedGenomes())[,1:5]
head(supportedGenomes(),n=12)[,1:4]
#12 lines, cols 1-4
pwf=nullp(genes,"hg19","ensGene")
?nullp
head(pwf)
#prob weighting func for each gene depending on its length (built in)
#give dataframe for all genes - name, DE or not (0=not, 1=true), bias.value based on length and diff expr, result is pwf, plotted
#x is bias score in 500 gene bins 

#http://rpackages.ianhowson.com/bioc/goseq/man/plotPWF.html
#The bias.data column is a quantification of the quantity for 
#which there is a bias in detecting DE for the associated gene, 
#this is usually gene length or the number of counts associated with a gene

#GO category overrepresentation amongst DE genes
GO.wall=goseq(pwf,"hg19","ensGene")
#fetches GO terms
head(GO.wall)
#use random sampling to generate the null distribution for category membership
GO.samp=goseq(pwf,"hg19","ensGene",method="Sampling",repcnt=1000)
head(GO.samp)
#limit analysis of GO to single category (GO molecular function)
GO.MF=goseq(pwf,"hg19","ensGene",test.cats=c("GO:MF"))
head(GO.MF)
enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue,
                                      method="BH")<.05]
#fdr correction
head(enriched.GO)

GOTERM
#load GO.db
library(GO.db)
GOTERM
#get info about each enriched term from GO.db
for(go in enriched.GO[1:5]){
  print(GOTERM[[go]])
  cat("------------------------------\n")
}

#KEGG pathway analysis
pwf=nullp(genes,"hg19","ensGene")
KEGG=goseq(pwf,"hg19","ensGene",test.cats="KEGG")
head(KEGG)

#cluster analysis
setwd("~/Course_Materials")
prDat <- read.table("Day4/RNAseq/GSE4051_data.tsv",
                    header=TRUE, row.names=1)
#load original normalised photo receptor gene exp data
#exp vals from 
str(prDat, max.level=0)
prDat
?readRDS
prDes<-readRDS("Day4/RNAseq/GSE4051_design.rds")
str(prDes)
sort(unique(prDes$sidNum))
unique(prDes$devStage)
unique(prDes$gType)
sprDat<-(t(scale(t(prDat))))
#str(sprDat, max.level=0, give.attr=FALSE)
round(data.frame(avgBefore=rowMeans(head(prDat)),
                avgAfter=rowMeans(head(sprDat)),
                varBefore=apply(head(prDat),1,var),
                varAfter=apply(head(sprDat),1,var)),2)

#clustering
#distance metric Euclidian
pr.dis<-dist(t(sprDat), method="euclidian")
head(pr.dis)
#create new factor representing integration of gType (genotype) and developmental stage
prDes$grp<- with(prDes, interaction(gType,devStage))
summary(prDes$grp)
#compute hierarchical clustering using different linkage types
pr.hc.s<- hclust(pr.dis, method='single')
pr.hc.c<- hclust(pr.dis, method='complete')
pr.hc.a<- hclust(pr.dis, method='average')
pr.hc.w<- hclust(pr.dis, method='ward.D2')
#plot
op<-par(mar=c(0,4,4,2),mfrow=c(2,2))
#plot multi
plot(pr.hc.s, labels=FALSE, main="Single", xlab="")
plot(pr.hc.c, labels=FALSE, main="Complete", xlab="")
plot(pr.hc.a, labels=FALSE, main="Average", xlab="")
plot(pr.hc.w, labels=FALSE, main="Ward", xlab="")
par(op)

#K-Means clustering
# say how many clusters you want to define
op<-par(mar=c(1,4,4,1))
plot(pr.hc.w, labels=prDes$grp, cex=0.6, main="Ward showing 10 Clusters")
rect.hclust(pr.hc.w, k=10)
par(op)

#heatmap
library(RColorBrewer)
library(xtable)
?colorRampPalette
GreyFun<- colorRampPalette(brewer.pal(n=9, "Greys"))
gTypeCols<- brewer.pal(n=11, "RdGy")[c(4,7)]
?brewer.pal
brewer.pal(n=11,"RdGy")
gTypeCols
heatmap(as.matrix(sprDat),
        Rowv=NA,
        col=GreyFun(256), 
        hclustfun=function(x) hclust(x, method='ward.D'),
        scale="none", 
        labCol=prDes$grp,
        labRow=NA,
        margins=c(8,1),
        ColSideColor=gTypeCols[unclass(prDes$gType)]
        )
legend("topright",
       legend=levels(prDes$gType),
       col=gTypeCols,
       lty=1,
       lwd=5,
       cex=0.5
       )
#legend not working? margins error?

GnBuFun<- colorRampPalette(brewer.pal(n=9, "GnBu"))
gTypeCols<- brewer.pal(n=9, "RdGy")[c(4,7)]
heatmap(as.matrix(sprDat),
        Rowv=NA,
        col=GnBuFun(256), 
        hclustfun=function(x) hclust(x, method='average'),
        scale="none", 
        labCol=prDes$grp,
        labRow=NA,
        margins=c(8,1),
        ColSideColor=gTypeCols[unclass(prDes$gType)]
)
legend("topright",
       legend=levels(prDes$gType),
       col=gTypeCols,
       lty=1,
       lwd=5,
       cex=0.5
)

set.seed(31)
k<-5
pr.km <-kmeans(t(sprDat), centers=k, nstart=50)

#look at sum of square of each cluster
pr.km$withinss

#look at comparison of each cluster
pr.kmTable <- data.frame(devStage=prDes$devStage, cluster=pr.km$cluster)
?xtable #create export table - convert to printable xtable object from R object
prTable<-xtable(with(pr.kmTable, table(devStage, cluster)),
                caption='Number of samples 
                from each developmental stage within each k-means cluster')
#align(prTable)<-"lccccc"
align(prTable)
#print(prTable, type='html', caption.placement='top')

#PAM algorithm