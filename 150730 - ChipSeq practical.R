library(GenomicAlignments)
library(ChIPQC)
library(GenomicRanges)
getwd()
dataDir<- "Day4/data_for_practical"
bamFiles<- dir(file.path(getwd(), dataDir), pattern="*.bam$", full.names=T)
?gsub # pattern matching and replacement
names(bamFiles)<- gsub(".bam","",basename(bamFiles))
bamFiles
?BamFileList #Rsamtools manage a list of BamFile instances
bfList<-BamFileList(bamFiles)
bfList

#extract SAM header for a single bam file
path(bfList)
samHeader <- scanBamHeader(path(bfList["TF_1"]))
str(samHeader,max.level=2)
#contains two components indicated by $ - targets and text
head(samHeader) # can miss is with this cf str

#explore the information provided in the sam header, how sorted, aligner used, what species
samHeader[[1]]$targets
samHeader[[1]]$text
samHeader[[1]]$text["@HD"]
samHeader[[1]]$text["@PG"]
samHeader[[1]]$text["@CO"]

#extract info about chr1
samHeader[[1]]
chr1dat <- samHeader[[1]]$targets["chr1"]
chr1dat
chr1range <- GRanges(seqnames=names(chr1dat),ranges=IRanges(1,chr1dat))
head(chr1range)
param<-ScanBamParam(which=chr1range)
param
#uses default scanBamParam flag argument
#could use to remove dups or artifacts

#having selected an area of interest, read data from single BAM file, convert to GR oject
#calc average length
alignDat <- readGAlignments(path(bfList["TF_1"]), param=param)
alignGR<- granges(alignDat)
seqlevels(alignGR)<-"chr1"
median(width(alignGR))

#MACS has options including broader peaks
#load the peaks called from MACS2
#NarrowPeak format is an extension of the BED format. UCSC webpage had infor on FAQs
mypeaks<-read.delim("Day4/data_for_practical/mypeaks_peaks.narrowPeak", header=F)
colnames(mypeaks)<-c("chrom","chromStart","chromEnd", "name", "score",
                     "strand", "fold.enrichment", "log10.pvalue", "log10.qvalue", "peak")
head(mypeaks)

#create genomic ranges object
library(GenomicRanges)
GRanges(mypeaks$chrom, IRanges(mypeaks$chromStart+1, mypeaks$chromEnd))

#setwd("/home/participant/Desktop/Course_Materials/")
#Exmaple usgin ChiPQsample()
library(ChIPQC)
bamFile<-file.path(getwd(), "Day4/data_for_practical/TF_1.bam")
exampleExp = ChIPQCsample(bamFile,peaks=NULL,annotation=NULL,chromosomes="chr1")
#speed up as just comp to chr1
class(exampleExp)
exampleExp
QCmetrics(exampleExp)
#.5M reads in example BAM file. Large number from suffic complex lib increases chance of finding true binding site
#optimal number of reads depends on enrichment, ab quality, fraction of the genome, visualize the variety of quality metrics
peaksFile<-file.path(getwd(), "Day4/data_for_practical/TF_1_peaks.bed")
data(blacklist_hg19)
exampleExp = ChIPQCsample(bamFile,peaks = peaksFile, blacklist = blacklist.hg19,
                          annotation="hg19",chromosomes="chr1")
QCmetrics(exampleExp)
#RiP=FRiP, RelCC rel x-cvg
plotCC(exampleExp)
#peak CC score is at 155bp corresponding to fragment size. 
fragmentlength(exampleExp)
#NB no camel case!
readlength(exampleExp)
FragmentLengthCrossCoverage(exampleExp)/ReadLengthCrossCoverage(exampleExp)
RelativeCrossCoverage(exampleExp)
#rel CC above 1 indicating successful ChIP
#RiBL is reads in blacklist

#FRIP
plotFrip(exampleExp)
frip(exampleExp)
#larger FRiP values indicate higher SNR
#>5% generally acceptable (unless very few expected binding sites when lower acceptable)

#distribution of signal in annotated regions
#regi=prop of reads in interval / prop of genome in interval ?genome i.e. features
regi(exampleExp)
plotRegi(exampleExp)

#hallmark of good ChIP quality is the inequality of coverage that occurs when there is enrichment of reads to certain portions of the genome
coveragehistogram(exampleExp)[1:10]
plotCoverageHist(exampleExp)
#not cutoff at 100bp for visualisation purposes
#if want to replot the whole histo you can extract the data using coveragehistogram func or extend x-axis per Advanced Topics
#shows significant stretch of high signal pile-up

plotSSD(exampleExp)
#slightly affected by blacklisted regions, taken together with coverage histogram indicates the TF chip has clear chip-signal above background
#signal pile up above that within artifact regions and this is pos assoc with TSS regions

#example with several bam files.
#paths relative to directory of samples.csv sheet, else need complete path
samples=read.csv("Day4/samples.csv")
samples
#exampleExp2 = ChIPQC(samples,annotation="hg19")
#computed with all BAM files and for all chromosomes - takes half an hour.
#save(exampleExp2, file="ChIPQCexample1.RData")
#precomputed result from ChIPQC
Rdata<-file.path(getwd(), "Day4/data_for_practical/ChIPQCexample1.RData")
load(Rdata)
QCmetrics(exampleExp2)
plotRap(exampleExp2,facetBy=c("Factor","Tissue"))
#vis reads in peaks. rel no reads that overlap peaks vs background reads
#vis is low %RiP for GATA3 in A549

#Assesing sample similarity with Diffbind
#for assessing ChiP quality - relate to correlation btw binding events across samples within an expt
#ChIPQC will perform clustinrg based on cooccurrence of peaks using methods available with DB
plotCorHeatmap(exampleExp2)
#as expected the GATA3 and CTCF samples form two distinct clusters

###
#loadl pre-computed ChIPQC objects.
Rdata<-file.path(getwd(), "Day4/data_for_practical/BCell_Examples.RData")
load(Rdata)
resExperiment
QCmetrics(resExperiment)
FragmentLengthCrossCoverage(resExperiment)
RelativeCrossCoverage(resExperiment)
plotCC(resExperiment,colourBy="Tissue",facetBy="Factor",lineBy="Replicate")
#some samples have much hight scores hence efficency, and fragment length appears v different from each other.
#DNAse sample for example has a frag length less than half of the PU1 sample
#Now established the diff in total efficiency, we can look at overall shape of cross coverage scores and relationship
#between signal peaks and artifact peaks in the cross coverage scores.
#To help better visualise we will apply a further facet wrap to the ggplot2 object in order to compare within factors
ccplot<-plotCC(resExperiment,colourBy="Tissue",facetBy="Factor",lineBy="Replicate")
ccplot +facet_wrap(~Factor,scales="free_y")
#see appropriate differential binding of Ikaros in DPT vs others - see text

plotFrip(resExperiment)
plotFribl(resExperiment)
plotRegi(resExperiment)
regi(resExperiment)["All3utrs",]
plotRegi(resExperiment)+scale_fill_gradient2(low="white",high="red",
                                             mid="white",midpoint=regi(resExperiment)
                                             ["All3utrs","Input_Ch12"])
#scale replace warn
plotCoverageHist(resExperiment,facetBy="Factor",colourBy="Tissue",lineBy="Replicate")
plotSSD(resExperiment)
plotCorHeatmap(resExperiment)
#visualise clustering in correlation heatmap based on correlation values for all the peaks in each sample.
