library(DESeq2)
library(GenomicFeatures)
library(Rsamtools)
library(GenomicAlignments)

#DESeq gets really angry when you have no replicates in your experiment, so...

#Read in sample coding table
sampleInfo <- read.table('sample_coding.txt')

#Read in count data
countdata <- read.table('RPE_counts.text', header=T)

#Remove extraneous rows in countdata
countdata <- countdata[,6:7]

#Duplicate data
countdata <- cbind(countdata, countdata)

#Create DESeq2 object
dds <- DeSeqDataSetFromMatrix(countData=countdata, colData=sampleInfo, design=~Treatment)

#Do Wald test because regular DESeq analysis will refuse to do the analysis
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

res <- results(dds)

write.table('DESeq_nbinomWaldTest.txt', res)