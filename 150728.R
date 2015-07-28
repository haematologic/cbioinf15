### 150728AM ###
library(IRanges)
ir <- IRanges(
)
?shift
library(Biostrings)
myseq <- DNAStringSet(randomStrings)
### Practical
library(BSgenome)
available.genomes()
library(BSgenome.Hsapiens.UCSC.hg19)
hg19 <- BSgenome.Hsapiens.UCSC.hg19
hg19
str(hg19)
dim(hg19)
seqnames(hg19)
names(hg19)
length(seqnames(hg19))
hg19[["chr22"]]
length(hg19[["chr22"]])
mylen <- length(hg19[["chr22"]])
mylen
?alphabetFrequency
alphabetFrequency(hg19[["chr22"]],"A")
alphabetFrequency(hg19[["chr22"]])
af <- alphabetFrequency(hg19[["chr22"]])
af
af
af[["A"]]
af[["G"]]
freqA <- af[["A"]]/mylen
freqG <- af[["G"]]/mylen
freqC <- af[["C"]]/mylen
freqT <- af[["T"]]/mylen
freqA
freqG
freqC
freqT
af[["G"]]/af[["N"]]
af <- alphabetFrequency(hg19[["chr22"]],baseOnly=TRUE)
af
af / mylen
#make sure last start position for 1000 base window is 1000 less then length of chr
startpos <- seq(1,mylen-1000,1000)
startpos
#1 based so last window is xxxx001 and ends 000
endpos <- startpos + 999
ir <- IRanges(startpos,endpos)
ir
#dont need to explicitly specify position, just start and width
#startpos is a vector of starting positions
ir <- IRanges(startpos,width=100)
ir
library(GenomicRanges)
rngs <- GRanges("chr22",ir)
rngs
sillyrng<- GRanges("sausages",ir)
sillyrng
winSeq <- getSeq(hg19,rngs)
winSeq
winAf <- alphabetFrequency(winSeq,baseOnly=TRUE)
winAf
head(winAf)
summary(winAf[,1])
#order for window with most bases first decreasing
aOrd <- order(winAf[,1],decreasing=TRUE)
head(winAf[aOrd,])
rngs[aOrd]
#which interval on the genome does the range with most bases correspond to
rngs[which.max(aOrd)] #? rubbish line
##
pdf("BaseDist.pdf",width=10,height=10)
plot(winAf[,1] ,type="n")
points(winAf[,1], pch="A",col="red")
points(winAf[,2], pch="C",col="green")
points(winAf[,3], pch="G",col="blue")
points(winAf[,4], pch="T",col="purple")
points(winAf[,5], pch="N",col="black")
dev.off()
ordByA <- winAf[order(winAf[,1],decreasing=TRUE),]
ordByA
#which regions have >60% GC content, which windows correspond
gcPerc <- (winAf[,2] + winAf[,3]) / 1000
gcPerc
head(gcPerc)
hist(gcPerc)
winAf
head(winAf)
winAf[gcPerc > 0.6,]
rngs[gcPerc > 0.6]
### Find sequences
s1 <- "aaaatgcagtaacccatgccc"
matchPattern("atg", s1)
# masked regions of the genome are represented by Ns
#find all occ of 10 consec Ns, create GRanges of result
n.matches <- matchPattern("NNNNNNNNNN", hg19[["chr22"]])
n.matches
maskRegions <- GRanges("chr22", IRanges(start=start(n.matches),end=end(n.matches)))
maskRegions
#reduce to non-overlapping regions
maskRegions <- reduce(maskRegions)
#get non masked regions
unmaskRegions <- gaps(maskRegions)
unmaskRegions
wholechr <- GRanges("chr22",IRanges(1,mylen))
setdiff(wholechr, maskRegions)
?setdiff
?gaps
## Part 2 AM practical
library(GenomicAlignments)
file.choose()
myfile <- file.choose()
mybam <- myfile
bam <- readGAlignments(file=mybam)
bam
bam <- readGAlignments(file=mybam,param=ScanBamParam(what=c("seq","mapq","flag")),use.names=TRUE)
bam
mcols(bam)
hist(mcols(bam)$mapq)
#filter mapq to gtt 30
bamFilt <- bam[mcols(bam)$mapq > 30]
bamFilt
hist(mcols(bamFilt)$mapq)
#see https://broadinstitute.github.io/picard/explain-flags.html
t<-table(mcols(bamFilt)$flag)
View(t)
# e.g. 1187     read paired, read mapped in proper pair, mate reverse strand, second in pair, read is PCR or optical duplicate
# e.g. 1024     read is PCR or optical duplicate
dupReads <- bamFilt[mcols(bamFilt)$flag==1024]
bamFilt <- bamFilt[mcols(bamFilt)$flag != 1024]
bamFilt
View(table(mcols(bamFilt)$flag))
bam.nodups <- readGAlignments(file=mybam, param=ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE)),use.names=TRUE)
names(dupReads) %in% names(bam.nodups)
names(dupReads) %in% names(bam)
names(dupReads)
# dody reads also discov by looking at base freq
alignedSeqs <- mcols(bamFilt)
alignedSeqs <- mcols(bamFilt)$seq
alignedSeqs
# coverage 2.2
cov <- coverage(bamFilt)
cov
names(cov)
cov[["22"]]
seqlevels(bamFilt)
#Rle is run lenght encoding system value of 0 repeated 16114211 times
#regions with more than 400 reads - just arbitary
slice(cov[["22"]],400)
HighCov <-ranges(slice(cov[["22"]],400))
#resize regions, some may only be a few bases long, want e.g. 100 bases long
#%over% shortcut for finding overlaps
HighCov
HighCov <- resize(HighCov,100)
HighCov <- reduce(HighCov)
HighCovGR <- GRanges("22",HighCov)
HighCovGR
?paste # concatenates strings
names(HighCovGR) <- paste("Region", 1:length(HighCov),sep="")
bam.sub <- bamFilt[bamFilt %over% HighCovGR]
bam.sub
length(bam.sub)
countOverlaps(HighCovGR,bamFilt)
sum(countOverlaps(HighCovGR,bamFilt))

#not informative
countOverlaps(bamFilt,HighCovGR)

bam.sub2 <- readGAlignments(file=mybam,use.names=TRUE,
                            param=ScanBamParam(which=HighCovGR,
                                               what=c("seq","mapq","flag")))
bam.sub2
# quicker to read partic genomic region with Which rather than creating subseq manually
bam.sub2 <- bam.sub2[mcols(bam.sub2)$mapq >30 &mcols(bam.sub2)$flag != 1024]
length(bam.sub)
length(bam.sub2)
all(names(bam.sub)==names(bam.sub2))
length(hg19[["chr22"]])

###
chrLen <- NULL
cnames <- seqnames(hg19)
# need to define seqnames - see ranges.R
cnames


### 150728 PM ###
N<- 60 # participants (variables)
N
## repeat block start
K<-2 # two possible treatments (levels)
# N/K 30 
G <- rep(c("Treatment A", "Control"), N/K)
#View(G)
#table(G)
G
Y <- rep(NA, N) #values for expression not yet controlled
Y
#when pt received treatment A give value 2
Y[which(G=="Treatment A")] <- 2
#when pt received control give value 1.5
Y[which(G=="Control")] <- 1.5
Y
# this is expected value
# add stochastic part - random bias - of biological differences between patients
# random error terms - normal distribution fits well, not biased
# generate random values n=N, from normal distribution with 0 mean, and 0.5 SD
# it's pseudorandom - can predict next number of sequence from previous
# can set seed so we all get the same
#set.seed(353)
error <- rnorm (N, 0, 0.5)
error
summary(error)
Y <- Y + error
Y
# tilde is symbol of linear model or statistical model in R
boxplot (Y ~ G)
## end block repeat
## repeat to see random varibility
#using baseline
BASELINE <- 1.5
#create second vector of expression that will be baseline for all
Y2 <- rep(BASELINE, N)
# add effect of differential expression (a)
Y2[which(G=="Treatment A")] <- BASELINE + 0.5
# add error as before
Y2 <- Y2 + error
boxplot (Y2 ~ G)
# no change, model same, just defined in terms of differential expression
# make design matrix
model.matrix ( ~ G) # uses intercept
model.matrix ( ~ 0 + G) # without intercept
###
set.seed(353)
N <- 120
K1 <- 2 # factor 1, 2 levels
K2 <- 2 # factor 2, 2 levels
G1 <- rep(c("Treatment A", "No Treatment"), N/K1)
G1
# do 3 into 2
G2 <- rep(c("ER-","ER+"), c(N/K2,N/K2))
G2
rep(c("Treatment A", "No Treatment"), 2)
rep(c("ER-","ER+"), c(2,2))
?rep
table(G1,G2)
BASELINE<-1
Y<- rep(BASELINE, N)
Y[which(G1=="Treatment A")]<-BASELINE +1
Y[which(G2=="ER+")]<-Y[which(G2=="ER+")] - 1.5
error<- rnorm(N,0,0.5)
Y
Y<- Y+error
Y
Y2<- rep(BASELINE, N)
Y2[which(G1=="Treatment A")]<-BASELINE +1 # effect of Rx A
Y2[which(G2=="ER+")]<-Y2[which(G2=="ER+")] -1.5 # effect of ER+
Y2[which(G1=="Treatment A" & G2=="ER+")]<-Y2[which(G1=="Treatment A" & G2=="ER+")] +1 -1.5 +3 # effects of Rx A, ER+, interaction
Y2
#get means of each of the groups to make graphs like factor a and b vs expression
interaction.plot(G1,G2,Y)
#sample size changes can make things different
interaction.plot(G1,G2,Y2) # see this looks more obvious the effect of Rx vs ER status
#try with different sample sizes
### Using lm function to estimate model SE etc.
?lm
m1 <- lm(Y ~ G1 + G2) # without interaction
m1
summary(m1)
m1 <- lm(Y ~ G1 * G2) # with interaction
m1
summary(m1)
#95%CI is estimate +-2*SE
#beta hat
m2 <- lm(Y2 ~ G1 + G2)
m2
summary(m2)
m2 <- lm(Y2 ~ G1 * G2)
m2
summary(m2)
?fitted
fitted(m1)
plot(m1)
par(ask=FALSE)
plot(fitted(m1))
?interaction
plot(fitted(m1)~interaction(G1,G2))
# fitted values yhat = beta hat . x
plot(resid(m1))
plot(resid(m1), pch=16, col="red")
#residuals (observed errors)
#allow for more complex data with non normal distributed (errors)

### Bowtie2 / TopHat2 alignment ###
library(Rsubread)
getwd()
setwd("/home/participant/Course_Materials/Day2/")
filesToCount <- dir("bam", pattern=".bam$", full.names=T)
tmp <- featureCounts(filesToCount, annot.inbuilt = "hg19", ignoreDup = F)
save(tmp, file="../Day3/countMatrix.RData")
