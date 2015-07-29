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
head(y) - has three vectors - counts, samples and genes
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
