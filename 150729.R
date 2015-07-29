rm(list=ls()) # clear environment
file.choose()
# "/home/participant/Course_Materials/Day3/countMatrix.RData"
setwd("~/Course_Materials")
library(edgeR)
load("Day3/Counts.RData")
Counts<- tmp$counts
Counts
colnames(Counts)
head(Counts)
#matrix of counts

# rename col names to shorter parameters for sample number, and tumours/normals
colnames(Counts) <- c("16N","16T","18N","18T","19N","19T")
colnames(Counts)
#str(tmp)
#head(tmp)
dim(Counts)
head(Counts)
rownames(Counts)

# create DE gene list object, provide the counts and the genes
# only counting reads mapped to regions containing genes (from 150728 TopHat)
dgList <- DGEList(counts=Counts, genes=rownames(Counts))
dgList
dgList$samples
colnames(dgList)
head(dgList$counts) #many rows!
head(dgList$genes) #Likewise!
#approx 26000 genes in list

#Filter genes if low counts, many not expressed or having enough reads to contribute
#e.g. 0 counts in normal, 1 count in tumour!
#retain only >1cpm
#must be blind to DE groups, else bias results between tumour and normal

?cpm
#work with CPM else would be biased by library size (take counts * library size / million)
countsPerMillion <- cpm(dgList)
summary(countsPerMillion)
#View(summary(countsPerMillion))
countCheck <- countsPerMillion >1
countCheck
#gives matrix of T or F with counts >1m
head(countCheck)

keep<- which(rowSums(countCheck)>=2)
#gives matrix where T or F has at least 2 cases (i.e. >2 samples being T or F), keep them
summary(keep)
head(keep)
length(keep) # only 14820 genes kept
#Bioc means can behave like list even though it is a matrix!
dgList<-dgList[keep,]
length(dgList)
View(summary(cpm(dgList)))
#Trimmed mean of M-values (TMM) normalise both within and between samples
?calcNormFactors
dgList<- calcNormFactors(dgList,method="TMM")
dgList$sample # new column added - weighting compared to total amount of RNA we have
# always working with raw counts (cpm), NOT divided by total number of reads
# sample norm.factors gives a new fixed offset to be used for our models and an intercept

plotMDS(dgList)
#like PCA but details different, nonetheless taking into account all differences in expression of all genes kept
#x dim very discrete T vs N, y dim unclear with only 3 samples
#plot inter-sample relationships based on multidimensional scaling

#Paried experiment - so need paired design matrix
#factor = 3individual (number of levels-1 so 2 cols), factor = tumour or normal (2 levels so 1 col)
sampleType<- rep("N",ncol(dgList)) #N=normal #T=tumour
sampleType
sampleType[grep("T",colnames(dgList))] <- "T" # apply T to sampleType 
sampleType
sampleReplicate <- paste("S", rep(1:3, each=2), sep="") # makes vector S1, S1, S2, S2, S3, S3
sampleReplicate 
designMat<-model.matrix(~sampleReplicate + sampleType) # can remove intercept by adding a 0 - see 150728
designMat
#sampleType
#dgList
View(designMat) #beautiful

#Per gene variability not respected by Poisson model
#Need overdispersed Poisson model -> this is the negative binomial
#Estimate dispersion for negative binomial model
#use empirical Bayes method to shrink the genewise dispersion estimates towards the common dispersion (tagwise dispersion)
dgList <- estimateGLMCommonDisp(dgList, design=designMat) # want variability within each group, not between samples?
dgList <- estimateGLMTrendedDisp(dgList, design=designMat) # genes more abundant may have more variability
dgList <- estimateGLMTagwiseDisp(dgList, design=designMat) # variability within each gene handled individually
#Plot biological coefficient of variation
?plotBCV
#The BCV is the sqrt of dispersion parameter in the negative binomial
plotBCV(dgList)
#housekeeping genes will have low variability so p value very small so will appear highly diff expressed
#housekeeping genes have BCV increased towards common/trend mean (in middle of plot middle)
#trend is better than common
#may therefore wish to use trended and tagwise only

#now find diff expressed genes
fit<-glmFit(dgList, designMat)
head(coef(fit)) 
# shows baseline, samples vs baseline, T vs N diff exp, on average, taking into account all replicates +in tumour vs normals. 
#Make sure in design matrix the 1 is in the tumour not normal else -ve values overexpressed in tumour

#Calculate ratios
?glmFit # genewise negative binomial generalised linear models
fit
lrt<-glmLRT(fit, coef=4)
head(coef(lrt))

?glmLRT # as above, conducts likelihood ratio tests for one or more coefficients in the linear model
#After fitting the model we can use topTags() fn to explore the results and set threshold to identify subset of DE genes
topTags(lrt)
#logFC -ve is diff exp in normal
# FDR vs P value: 0.05 rule for P value not correct as error way higher - so do multiple testing doing FDR
# need to see number of genes good enough to live with expected number of false positives to set the FDR
# if select genes for more testing, depends how many you want to test
# FDR - some publications 0.1, some 0.05 
# i.e. expect 1:10 false positives - no good if only have 5 genes as 2 false pos!, 
# vs 1000 genes where you would have 100 false +ve
# FDR is just a measure of expected number of false positives
# it is different to just a significance value
edgeR_result <- topTags(lrt)

?topTags
res<-topTags(lrt)$table
save(res, file='Day3/edgeR_Result.RData')
View(res)
# finally plot the log-fold changes of all the genes and highlight those diff expressed!
?decideTests #multiple testing across genes and contrasts - classity a series of related t-stats as up down or n/s.
deGenes <- decideTestsDGE(lrt,p=0.001)
rownames(lrt) # Ensembl gene IDs
deGenes <- rownames(lrt)[as.logical(deGenes)]
?plotSmear #edgeR method same as MA plot - de.tags = rownames for gene identified as being diff expressed, use exact test or glmLRT
# MA allows examination of the relationship between intensity and difference between two data store (currently selected prob list)
# x = average quantitated value across the currently selected probe list, y is the difference between
plotSmear(lrt, de.tags=deGenes)
abline(h=c(-1,1), col=2)
?cpm
#CPM axis = counts per million or reads per kb per million values

topTags(lrt) # shows top DE with coefficient sampleTypeT
# geme 11420 has significant log fold change in tumours vs normals (underexpressed as negative)
# find location genome to view IGV and reads!
tmp$annotation[1,]
which(tmp$annotation$GeneID == "284254")
#tmp[11420,] # error
tmp$annotation[11420,]
#copy output chr18 52258390 -> add colon chr18:52258390
#load IGV
#load Bam files from tumour and normals
#enter location


## supplementary
rm(list=ls())
file.choose()
rawdata <- read.delim("Day3/rnaSeq/OralCarcinoma/TableS1.txt", check.names=FALSE, stringsAsFactors=FALSE)
head(rawdata)
library(edgeR)
# for easy manipulation put the data into a DGEList object
y<- DGEList(counts=rawdata[,4:9], genes=rawdata[,1:3])
dim(y) #15668
#Old study so not all RefSeq IDs provided match RefSeqIDs currently in use. We retain only those transcripts with IDs in the current NCBI annotation, which is provided by org.Hs.eg.db package.
library(org.Hs.eg.db)
head(mappedRkeys(org.Hs.egREFSEQ))
head(y$genes)
head(y$genes$RefSeqID)
idfound<-y$genes$RefSeqID %in% mappedRkeys(org.Hs.egREFSEQ)
y<- y[idfound,]
dim(y) #15559 genes left

#?toTable - don't ask
#Add entrez gene IDs to the annotation 
egREFSEQ <- toTable(org.Hs.egREFSEQ)
summary(egREFSEQ)
head(egREFSEQ)
m<- match(y$genes$RefSeqID, egREFSEQ$accession)
m
y$genes$EntrezGene # null
y$genes$EntrezGene<-egREFSEQ$gene_id[m]
dim(y)

y$genes$EntrezGene # null
#Use the egIDs to update the gene symbols
egSYMBOL <-toTable(org.Hs.egSYMBOL)
head(egSYMBOL)
m<- match(y$genes$EntrezGene, egSYMBOL$gene_id)
y$genes$Symbol # not null - gene names (Symbols)
y$genes$Symbol<-egSYMBOL$symbol[m]
head(y$genes)

# different RNA transcripts for the same gene symbol count predominantly the same reads. 
# So we keep one transcript for each gene symbol. We choose the transcript with the highest overall count.
head(y$counts)
o<-order(rowSums(y$counts))
head(o)
head(y) # has three vectors - counts, samples and genes
y<-y[o,]
head(y)
nrow(y)
?duplicated # BiocGenerics method
d<-duplicated(y$genes$Symbol)
y<-y[!d,]
nrow(y)

#normally we would also filter lowly expressed genes. All transcripts here already have at least 50 reads for all samples of at least one of the tissue types
# check library sizes
y$samples$lib.size
#use Entrez Gene IDs as row
y$genes$EntrezGene
y$genes
rownames(y$genes)
rownames(y$genes)<-y$genes$EntrezGene
rownames(y$counts) <- rownames(y$genes) <- y$genes$EntrezGene
dim(y$genes$EntrezGene)
y$genes
summary(isNA(y$genes$EntrezGene))
which(isNA(y$genes$EntrezGene))
y$genes$EntrezGene

#TMM normalisation
y<- calcNormFactors(y)
y$samples

plotMDS(y)
?factor
Patient<- factor(c(8,8,33,33,51,51))
Tissue<- factor(c("N","T","N","T","N","T"))
data.frame(Sample=colnames(y),Patient,Tissue)
design<-model.matrix(~Patient+Tissue)
rownames(design)<-colnames(y)
colnames(y)
head(design)
design
y<-estimateGLMCommonDisp(y, design, verbose=TRUE) #adds $common dispersion
head(y) 
length(y$common.dispersion)
y<-estimateGLMTrendedDisp(y,design)
head(y)
length(y$trended.dispersion)
y<-estimateGLMTagwiseDisp(y,design)
head(y)
length(y$tagwise.dispersion)
plotBCV(y)

#Determine differentially expressed genes
#fit genewise generalised linear models
?glmFit # genewise negative binomial GLM
fit<-glmFit(y, design) #fit a negative binomial generalised log-linear model to the readcounts for each gene
?glmLRT #produces also table object (df with same rows as y containing the log2-fold changes, likelihood ratio stats and p values ready to be displayed by topTags.. and comparison object with char string of the coefficient or contrast being tested
lrt<-glmLRT(fit) # has conducted test for last coeffienct in the linear model, tumour v normal 
lrt
topTags(lrt)
colnames(design)

o<-order(lrt$table$PValue)
cpm(y)[o[1:10],]
summary(de <- decideTestsDGE(lrt))
detags<-rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-1,1),col="blue")
#for fun 
heatmap(design, de.tags=detags)

### 150729 PM ###
library(biomaRt)
head(listMarts(),5)
ensembl<-useMart("ensembl")

ensembl
ensembl<-useMart("ensembl", dataset="hsapiens_gene_ensembl")
ensembl
head(listDatasets(ensembl),10)
head(listFilters(ensembl),5)
flt<-listFilters(ensembl)
dim(flt) #use chromosome_name to query
flt
flt[grep("entrez",flt[,1]),]
#myinfo can be slow to get for 100s or 1000s
head(listAttributes(ensembl),25)
entrez<-c("673","837")
myfilter<-"entrezgene"
attr=c("entrezgene","hgnc_symbol","ensembl_gene_id","description")
allAttr<-listAttributes(ensembl)
attr %in% allAttr[,1]
myInfo<-getBM(filters="entrezgene", values=entrez, attributes=attr,mart=ensembl)
myInfo
myfilters <- c("chromosome_name", "start", "end")
#filter by called start, later attribute called start_position
myvalues <- list(16, 1100000, 1250000)
head(allAttr[grep("start", allAttr[,1]),])
attr<-c("ensembl_gene_id","hgnc_symbol","entrezgene",)
#chromosome name is 16 not chr16
#biomart server can block you. or be down. and you need to be online
#Bioconductor packages then available instead, every 6 months
#org.Hs.eg.db
#also loads DBI - deals with databases
#bigger than standard package as it's a database
#uses keytypes vs filters
keytypes(org.Hs.eg.db)
#has fewer, but the common types
#gives chromosome location, number and strand

#may wish to know more about gene structure 
# e.g. genomic features - take data from biomaRt and others and makes gene models for each organism
library(TxDb.Hsapiens.UCSC.hg19.knownGene) #transcript database for given organism and resource version
txdb<-TxDb.Hsapiens.UCSC.hg19.knownGene
txdb
keytypes(txdb)
columns(txdb)
select(txdb,keys=entrez, keytype="GENEID",columns=c("TXID","TXCHROM","TXSTART","TXEND"))
#exons reversed in genomic coordinate order?
# if we have bed file or gff file - packages rtracklayer imports these
# bed file is gene location
# wig file
# gtf file
# result is GRanges object
# for targeted sequencing check this has been efficient
# Roche NimbleGen

## 150729PM Practical ##
library(biomaRt)
listMarts()
ensemblLatest<-useMart("ensembl")
datasets<-listDatasets(ensemblLatest)
datasets[which(datasets[,1]=="hsapiens_gene_ensembl"),]
#unfortunately our RNAseq was aligned and annotated against hg19. Different versions of biomaRt and genome vesion can have diff attributes and underlying annotation
#connect to Feb 2014 hg archive - this result was from Mark from google
ensembl_75 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org",path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
datasets <- listDatasets(ensembl_75)
#verify the version is GRCh37.p13
datasets[which(datasets[,1]=="hsapiens_gene_ensembl"),]
load("Day3/topHits_FWER_0.01.RData")
head(topHits)
dim(topHits)
topHits<-topHits[1:100,] # just taking 
length(topHits)
dim(topHits)
#use biomaRt to annotate the results of our RNAseq analysis.
#the first col in the topHits dataframe are EntrezIDs.
flt<- listFilters(ensembl_75)
head(flt)
flt
#this is the filter for entrez gene ID
flt[grep("entrez",flt[,1]),]
attr<-listAttributes(ensembl_75)
attr[1:50,]
attr[grep("symbol",attr[,1]),] # hgnc_symbol is the name to use for HGNC symbol
extraInfo <- getBM(attributes = c("entrezgene","hgnc_symbol","description"),
                   filters="entrezgene",values=topHits[,1], mart=ensembl_75)
head(extraInfo) #annotate the topHits
# !careful - it has orderd by alphanumeric not by p-value as obtained from topHits!
annotatedHits<-merge(topHits, extraInfo, by.x =1, by.y =1, sort=FALSE)
# specify column in these two dataframes that have the same identifier in x=1 and y=1
# not resorted based on extra info, but as it was in topHits
dim(annotatedHits) # has grown by 1
dim(topHits)
which(duplicated(annotatedHits[,1])) # which is duplicated
annotatedHits[80:85,] # 82 is the first instance, 83 the second
head(topHits) # gene column is 1
head(annotatedHits) # gene column is 1
#the resulting annotatedHits table can be used to ease the interpretation of the RNA-seq analysis and can be shared with collaborators
dim(annotatedHits)
tail(annotatedHits)

#retrieve gene sequences
head(topHits)
head(topHits[,1])
posInfo<-getBM(attributes=c("entrezgene","chromosome_name","start_position",
                            "end_position","strand"),
               filters="entrezgene",
               values=topHits[,1],
               mart=ensembl_75)
posInfo
head(posInfo[posInfo[,2] %in% c(1:22, "X","Y"),])
posInfo<-posInfo[posInfo[,2] %in% c(1:22, "X","Y"),] # get rid of extra unassembled chromosomes, in practice wouldnt use biomaRt
head(posInfo)
#NB chr number is a numeric value not chr-prefixed string and strand is 1 or -1 not + -

#Create GRanges representation of these gene coordinates. Make sure that the sequence names and strand info is acceptable.
library(GenomicRanges)
posInfo[1:10,5]
strand=ifelse(posInfo[,5]==1, "+", "-") #clever
strand[1:10]
#paste0("chr",posInfo[,2]) - is adding chr to every posInfo[,2]
genePos <- GRanges(paste0("chr",posInfo[,2]),
                   IRanges(posInfo[,3],posInfo[,4],names=posInfo[,1]),strand)
#paste0 concatenates with no spaces, default paste adds a space
genePos
#with function allows columns from a data frame to be accessed by name vs column index- result the same
genePos<- with(posInfo, GRanges(paste0("chr",chromosome_name),
                                IRanges(start_position, end_position, names=entrezgene),
                                strand))
#paste0 has no spaces
genePos

#load package providing seq for hg19 and get DNA seq for the top hits. Translate into aa
library(BSgenome.Hsapiens.UCSC.hg19)
hg19<-BSgenome.Hsapiens.UCSC.hg19
myseqs<-getSeq(hg19,genePos)
myseqs
translate(myseqs)

#get seq 100 bp upstream of each gene and write these to a fasta file
promSeqs <- getSeq(hg19, flank(genePos,100))
writeXStringSet(promSeqs , file="topGenesPromoters.fa")
promSeqs
#View(promSeqs)

##annotation using pre built packages
library(org.Hs.eg.db)
hs<-org.Hs.eg.db
columns(hs)
keytypes(hs)
chr22Genes <- select(hs, columns="ENTREZID", keytype="CHR", keys="22")
#complains about CHR being deprecated
head(chr22Genes)
chr22ID<-chr22Genes[,2]
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
tx<-TxDb.Hsapiens.UCSC.hg19.knownGene
exo<-exonsBy(tx,"gene")
#takes a sec
exo
length(exo)
chr22ID
exo22<-exo[names(exo) %in% chr22ID]

#import the aligned reads for sample "16N" using GenomicAlignments
library(GenomicAlignments)
bam<-readGAlignments(file="Day2/bam/16N_aligned.bam",
                     param=ScanBamParam(which=GRanges("chr22", IRanges(1,51304566))),
                     use.name=TRUE)
bam
#now can overlap reads with exons
#findOverlaps will report the amount
#start with overlaps of single gene to make sure we understand output
counts<- countOverlaps(exo22[["4627"]],bam)
counts
so<- summarizeOverlaps(exo22, bam)
# extends findOverlaps by providing options to resolve reads that operlap multiple features.
# Each read is counted a maximum of once. Different modes of counting available.
# ?summarizeOverlaps
so

?assays
gCounts<-assays(so)$counts
dim(gCounts)
head(gCounts)

#visualisation
geneRegions<- unlist(range(exo))
#select names of statistically significant genes from the edgeR output in the usual manner
load("Day3/topHits_FWER_0.01.RData")
sigResults<-topHits
head(sigResults)
sigGeneRegions<- geneRegions[na.omit(match(sigResults[,1],names(geneRegions)))]
# for some reason some of the genes reported in the edgeR output do not have associated genomic features. Hence we have to discard such genes from the output using na.omit
head(sigGeneRegions)

#the bed format is also able to colour each range according to some property of the analysos (e.g. direction or magnitude of change)
# a score can also be displayed when a particular region is clicked on.

#can attach metadata to GRanges
head(sigResults[,1])
mcols(sigGeneRegions)<-sigResults[match(names(sigGeneRegions), sigResults[,1]),]
#stores adjusted p-values and log fold-change in the ranges object
head(sigGeneRegions)
head(mcols(sigGeneRegions))

# use the human organism package to add gene symbols
sigGeneRegions$genes
anno <- select(org.Hs.eg.db,keys=sigGeneRegions$genes,
               keytype="ENTREZID", columns="SYMBOL")
head(anno)
tail(anno)
mcols(sigGeneRegions)$Symbol<-anno[,2]
mcols(sigGeneRegions)$Symbol

#tidy up the genomic ranges object so that only chr 1 to 22 and x y are included
seqlevels(sigGeneRegions)
?seqlevels
sigGeneRegions<-keepSeqlevels(sigGeneRegions, paste0("chr",c(1:22,"X","Y")))

#create a score from the p-values that will be displayed under each region and colour scheme for the regions based on the fold change. For convenience restrict fold changes to be within the region -3 to 3
Score<- -log10(sigGeneRegions$FWER)
rbPal<- colorRampPalette(c("red", "blue"))

logfc<- pmax(sigGeneRegions$logFC, -3)
logfc<- pmin(logfc, 3)

Col<- rbPal(10)[as.numeric(cut(logfc, breaks=10))]

#the colours and score can be saved in the GRanges object and score and itemRgb columns respectively. 
mcols(sigGeneRegions)$score<-Score
mcols(sigGeneRegions)$itemRgb<-Col
library(rtracklayer)
export(sigGeneRegions, con="topHits.bed")
#can also export bed wig gff

#ggplot2
load("Day3/edgeRAnalysis_ALLGENES.RData")
head(y)
y$significance<- -10*log10(y$FDR)
edgeRresults<-y

library(ggplot2)
ggplot(edgeRresults, aes(x=logFC, y=significance)) + geom_point()
edgeRresults$DE <- edgeRresults$FDR < 0.01
ggplot(edgeRresults, aes(x=logFC, y=significance,col=DE)) + geom_point()
#make an MA plot
ggplot(edgeRresults, aes(x=logCPM, y=logFC, col=DE)) + geom_point()
ggplot(edgeRresults, aes(x=logCPM, y=logFC, col=DE)) + geom_point(alpha=0.4) + scale_color_manual(values = c("black","red"))
ggplot(edgeRresults, aes(x=logCPM, y=logFC, col=DE,size=significance)) + geom_point(alpha=0.4) +
  scale_colour_manual(values=c("black","red"))

library(ggbio)
plotGrandLinear(sigGeneRegions, aes(y=score))
mcols(sigGeneRegions)$Up <- logfc> 0
plotGrandLinear(sigGeneRegions, aes(y=logFC, col=Up))
myGene <- which(seqnames(sigGeneRegions)=="chr22" & mcols(sigGeneRegions)$logFC <0)[1]
sigGeneRegions[myGene]
entrez<- sigGeneRegions$genes[myGene]
#plot gene structure
autoplot(tx, which=exo22[[entrez]])
myreg<-flank(reduce(exo22[[entrez]]), 1000, both=T)
bamSubset<-bam[bam %over% myreg]
#plot numer of reads in this region
autoplot(bamSubset, which=myreg)
#repeat plot with a smoothed coverage
autoplot(bamSubset , stat="coverage")
#make ideogram for Chr22
idPlot<-plotIdeogram(genome="hg19",subchr="chr22")
idPlot
#make combined plot of aligned reads and transcripts for selected gene
geneMod<-autoplot(tx,which=myreg)
reads1<-autoplot(bamSubset, stat="coverage")
tracks(idPlot,geneMod,reads1)

#compare coverage of several samples on the same plot.
#loop for a vector of file names
mytracks<-alist()
fls<-paste0("Day2/bam/", dir("Day2/bam/",pattern=".bam"))
fls<-fls[-grep("bai",fls)]
names(fls)<-basename(fls)
for (i in 1:length(fls)){
  bam<-readGAlignments(file=fls[i],
                       param=ScanBamParam(which=GRanges("chr22",IRanges(1,51304566))),
                       use.name=TRUE)
  bamSubset<-bam[bam %over% myreg]
  mytracks[[i]]<-autoplot(bamSubset, stat="coverage") + ylim(0,20)
}
tracks(idPlot,geneMod,mytracks[[1]], mytracks[[2]], mytracks[[3]], mytracks[[4]], mytracks[[5]], mytracks[[6]])
#doesnt display separate legend for each bam file