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
names(HighCovGR) <- paste()
# quicker to read partic genomic region with Which rather than creating subseq manually
