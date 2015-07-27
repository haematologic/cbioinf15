sampleMat <-pData(sample.ExpressionSet)
dim(sampleMat)
head(sampleMat)
hist(evals[,1])
boxplot(evals[,1:5])
plot(evals[,1],evals[,2])
plot(evals[1,],evals[2,])
M <- log2(evals[,1]) - log2(evals[,2])
A <- 0.5*(log2(evals[,1]) + log2(evals[,2]))
plot(A,M)
M <- log2(evals[,1]) - log2(evals[,2])
A <- 0.5*(log2(evals[,1]) + log2(evals[,2]))
plot(A,M)
x<1:10
x<-1:10
y<-2*x
plt(x,y)
plot(x,y)
x <- runif(100)
hist(x)
dd <- data.frame(A=rnorm(100, mean = 2),
                 B = rnorm(100, mean=5))
boxplot(dd)
runif(100)
man(runif)
plot(x,y,xlab="My X label",ylab="My Y label"
     ,main="My Title",col="steelblue",pch=16,xlim=c(0,20))
points(x,x,col="red",pch=17)
grid()
abline(0,2,lty=2)
abline(0,1,lty=2)
text(12,10, label="y = x")
text(12,20, label="y = 2x")
values <- rnorm(10)
?runif
values <0
values >0
values
data <- data.frame(Counts = values, Name = rep(c("A","B")))
datga
data
data[values>0,]
data[data$Name == "A"]
data[data$Name == "A",]
data[data$Name != "A",]
data[data$Name == "A" & values > 0,]
data[data$Name == "A" | values > 0,]
data$Name == "A" | values > 0
.libPaths
.libPaths()
.Library
for (i in 1:3){hist(,i)}
for (i in 1:3){hist(data,i)}
data
for (i in 1:3){hist(data[,i])}
xx<-c(runif(10))
xx
for (i in 1:3){hist(xx,i])}
for (i in 1:3){hist(xx,i)}
for (i in 1:3){hist(xx[,i])}
pwd
getpwd
getwd
getwd()
library(Biobase)
data(sample.ExpressionSet)
sample.ExpressionSet
eval <- exprs(sample.ExpressionSet)
dim (eval)
evals <- exprs(sample.ExpressionSet)
dim (evals)
evals [1:4,1:3]
getwd
getwd()
setwd("~/Course_Materials/Day1/nki")
myfile <- file.choose()
myfile
dim(myfile)
file.exists(myfile)
expression <- read.delim(myfile)
dim(expression)
head(expression)
expression[1:5,1:10]
expression[1:10,]
expression[1:10,1]
subset <- expression [1:10,]
subset
#remove a row containing identifier that could be misinterpreted as numeric quantity
rownames(expression) <- expression[,1]
rownames(expression)
evals <- expression[,-1]
evals
evals <- read.delim(myfile,row.names=1,as.is=TRUE)
?as.is
?as
?is
dim(evals)
dim(expression)
head(evals)
evals[1:10,]
evals[1:5,1:10]
#methods other than : for generating row indices
expression[c(1,3,5,7),]
expression[seq(1, 10000,by=1000),]
?by
?seq
c(seq(1,10000,by=1000))
c(seq(0,10000,by=500))
c(seq(0,10000,5000))
dim(sample)
?sample
expression[sample(1:24481,10),]
?matrix
evals<- as.matrix(evals)
evals[1,]
#print gene expression values for genes in the first row
mean(evals[1,])
summary(evals[1,])
hist(evals[1,])
mean(evals[1,])
summary(evals[1,])
plot(evals[1,])
colors()
plot(evals[1,],col="steelblue")
rainbow(n=ncol(evals))
#vector of same length as number of cols
plot(evals[1,],col=rainbow)
plot(evals[1,],col=rainbow(n=ncol(evals)))
rainbow(n=ncol(evals))
plot(evals[1,],col=rainbow(n=ncol(evals)),pch=16)
plot(evals[1,],col=rainbow(n=ncol(evals)),pch=1:16)
plot(evals[1,],col="steelblue",pch=16,xlab="Patient",ylab="Expression value",main="Expression of gene X")
which(evals[1,]>4)
#expression greater than 4 for the first gene
values <- evals[1,]
values
plot(values,col="steelblue",pch=16)
?outl - #variable name not a function
outl <- which(abs(values) >4)
#abs takes - or +
outl
abline(h=c(-4,4))
points(outl,values[outl],col="red",pch=16)
limit <- 4
values <- evals[1,]
plot(values,col="steelblue",pch=16)
outl <- which(abs(values) > limit)
abline(h=c(-limit,limit))
points(outl,values[outl],col="red",pch=16)
?points
#points added to plot x,y co-ords, type
clinical <- file.choose()
clinical
clinical <- read.delim("NKI295.pdata.txt")
dim(clinical)
head(clinical)
features <- file.choose()
features <- read.delim("NKI295.fdata.txt")
head(features)
dim(features)
dim(clinical)
dim(expression)
dim(evals)
#Use Bioconductor
library(Biobase)
eset <- readExpressionSet("NKI295.exprs.txt", "NKI295.pdata.txt")
?readExpressionSet
eset
#shows summary
dim(exprs(eset))
head(exprs(eset))
dim(pData(eset))
head(pData(eset))
tail(pData(eset))
#esets designed to behave like other basic R data types
eset[1:1000,]
#gives 1000 features, all samples
eset[,1:10]
#gives all features, 10 samples
?sample
sampleset <- sample(eset,1000)
#no
?nrow
#number of rows
#subset random 1000 genes from first 50 samples
random.rows <- sample(1:nrow(eset),1000)
random.rows
eset.sub <- eset[random.rows,1:50]
eset.sub
dim(exprs(eset.sub))
head(exprs(eset.sub))
pData(eset.sub)
exprs(eset.sub)
#subset acc to sample annotation
sampleMat <- pData(eset)
sampleMat
head(sampleMat)
rownames(sampleMat)
colnames(sampleMat)
all(rownames(sampleMat) == colnames(exprs(eset)))
?all
table(sampleMat$ER)
#well curated dataset where no 'NEGATIVE' or 'negative' just 'Negative'
sampleMat$ER == "Negative"
erNegSet <- eset[,sampleMat$ER == "Negative"]
#all rows, cols which have ER == 'Negative'
erNegSet
View(erNegSet)
features <- read.delim ("NKI295.fdata.txt")
head(features)
eset<-readExpressionSet("NKI295.exprs.txt", "NKI295.pdata.txt")
eset
myrow <- match("ESR1", features[,2])
mygene <- exprs(eset)[myrow,]
mygene
plot(density(mygene))
?pData
ERStatus <- pData(eset)$ER
boxplot(mygene~ERStatus)
t.test(mygene~ERStatus)
library(genefilter)
tstats <- rowtests(exprs(eset), as.factor(ERStatus))
tstats <- rowttests(exprs(eset), as.factor(ERStatus))
#open rowttests
head(tstats)
ordt<-order(abs(tstats), decreasing=TRUE)
ordt
topGenes <- ordt[1:10]
features[topGenes,]
topGenes
pdf("TopGenes.pdf")
for(i in 1:10){
  rowind <- ordt[i]
  geneName <- features[rowind,2]
  boxplot(exprs(eset)[rowind,]~ERStatus,main=geneName)
}
dev.off()
### Mark refresh

?rownames
#use rownames not numbers e.g. rownames(expression) <- expression[,1]
#remember to give row and column index else output not as expected
seq(1,1000,(100-10))
?seq
seq(5,100)
example(hist)
?points
###
View(sampleMat)
#prefer to head fn?
sampleMat$ER
features 
features [1:10,]
match("ESR1", features[,2])
match("ESR1", features[,2])
features[match("ESR1", features[,2]),]
myrow <- match("ESR1", features[,2])
mygene <- exprs(eset[myrow,])
mygene
plot(density(mygene))
boxplot(mygene)
boxplot(mygene~myrow)
head(pData(eset))
ERStatus <- pData(eset)$ER
ERStatus
boxplot (mygene ~ ERStatus)
t.test(mygene ~ ERStatus)


#### 150727-PM ###


