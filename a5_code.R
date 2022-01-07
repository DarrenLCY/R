rm(list = ls())

#Loading the R packages
library(limma)
library(edgeR)

setwd("C:/Users/Owner/Dropbox/Memorial University - Year 3 Semester 1/BIOL 3951/Assignment 5")

## Reading the files
sampleInfo <- read.table("Sample_Info.csv", sep = ",", header = TRUE, stringsAsFactors = F, row.names =1)
sampleInfo$shortName <- paste(sampleInfo$MouseType, 1:3, sep = "")
sampleInfo

readCounts <- read.table("GSE37236_processed_data_file_normalized_counts.tab.txt", sep ="\t", 
                         header = T, stringsAsFactors = F, row.names =1)
dim(readCounts)
head(readCounts)
tail(readCounts)

counts <- readCounts[,-7] #don't want the last column
colnames(counts) <- sampleInfo$shortName
head(counts)
dim(counts)

nsamples <- ncol(counts)

## Processing the data

dge <- DGEList(counts = counts)

avgCPM <- aveLogCPM(dge)
head(avgCPM)
ToKeep <- avgCPM >= 1
sum(ToKeep)

lcpm <- cpm(dge, log = T)

par(mfrow=c(1,2))

# Plotting unfiltered data
plot(density(lcpm[,1]), lwd=2, ylim = c(0,0.225),las=2, main="Raw data", xlab="Log-cpm")

for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y,  lwd=2)
}
abline(v=1, lty=3)

#Filtering unexpressed genes
dge <- DGEList(counts = counts[ToKeep,])
lcpm <- cpm(dge, log = T)
#Plotting filtered data
plot(density(lcpm[,1]), lwd=2, las=2, ylim = c(0,0.225), main="Filtered data", xlab="Log-cpm")
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y,  lwd=2)
}
abline(v=1, lty=3)

#Summary filtered data
summary(lcpm)

#Visualizing what normalization does
par(mfrow=c(1,2))
boxplot(lcpm, las=2, main="Unnormalised data",ylab="Log-cpm", cex.axis = 0.5)

dge <- calcNormFactors(dge, method = "TMM") #Normalized using the method of trimmed mean of M-values
lcpm <- cpm(dge, log = T)
boxplot(lcpm, las=2, main="Normalised data",ylab="Log-cpm", cex.axis = 0.5)
dge$samples$norm.factors

#Creating the design and contrast matrix
design <- model.matrix(~0+sampleInfo$MouseType)
design
colnames(design) <- c("KnockOut", "WildType")
design

contr.matrix <- makeContrasts(
   KvsW = KnockOut-WildType, 
   levels = colnames(design))
contr.matrix

#Removing Mean-variance trend 
par(mfrow=c(1,2))
v <- voom(dge, design, plot=TRUE) # estimate the mean-variance relationship
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
#Differential expression with eBayes
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

head(decideTests(efit))
summary(decideTests(efit))

topTable(efit, coef = "KvsW")

counts["ENSMUSG00000027639",]

#DE visualization.
dt <- decideTests(efit)

par(mfrow=c(1,1))
plotMD(efit, column=1, status=dt[,1], main=colnames(efit)[1], xlim=c(0,13), ylim=c(-3,3))

library(gplots)

KvsWgenes <- topTable(efit, coef=1, n=100)
KvsW.topgenes <-  row.names(KvsWgenes)
mycol <- colorpanel(100,"blue","white","red")
heatmap.2(lcpm[KvsW.topgenes,], scale="row", labRow= KvsW.topgenes,  
          col=mycol, trace="none", density.info="none",  margin=c(8,10), 
          lhei=c(3,10), dendrogram="column")

#### To write all the results to a file
write.fit(efit, dt, file="results.txt")
