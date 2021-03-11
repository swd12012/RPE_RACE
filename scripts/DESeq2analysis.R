library(DESeq2)
library(GenomicFeatures)
library(Rsamtools)
library(GenomicAlignments)

#DESeq gets really angry when you have no replicates in your experiment, so...

#Read in sample coding table
sampleInfo <- read.table('RACEdesign.txt', header=T)

#Read in count data
countdata <- read.table('RACE_counts.txt', header=T, row.names=1)

#Remove extraneous rows in countdata
countdata <- countdata[,6:7]

#Duplicate data so DESeq2 doesn't throw a fit
countdata <- cbind(countdata, countdata)
colnames(countdata) <- c('P1-1','P2-1', 'P1-2','P2-2')

#Read in annotation file
annotations <- read.table('mouseENSEMBLID_genename.txt',header=T,sep=',')

annotated_countdata <- merge(countdata,annotations, by.x='row.names',by.y='Gene.stable.ID')

#Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=sampleInfo, design=~Condition)

#Do Wald test because regular DESeq analysis will refuse to do the analysis
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
dispersions(dds) <- mcols(dds)$dispGeneEst
dds <- nbinomWaldTest(dds)

res <- results(dds)

write.table('DESeq_nbinomWaldTest.txt', res)

#Exploratory data analysis
#plotMA( res, ylim = c(-1, 1) )
#plotDispEsts( dds )
hist( res$pvalue, breaks=20, col="grey" )
###  throw out lowly expressed genes?? ... I leave this as an exercise
###  add external annotation to "gene ids"
# log transform
rld = rlog( dds )
## this is where you could just extract the log transformed normalized data for each sample ID, and then roll your own analysis
head( assay(rld) )
mydata = assay(rld)

#sampleDists = dist( t( assay(rld) ) )
# heat map
#sampleDistMatrix = as.matrix( sampleDists )
#rownames(sampleDistMatrix) = rld$TissueCode
#colnames(sampleDistMatrix) = NULL
library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#heatmap.2( sampleDistMatrix, trace="none", col=colours)
# PCs
# wow you can sure tell tissue apart
#print( plotPCA( rld, intgroup = "TissueCode") )
# heat map with gene clustering
library( "genefilter" )
# these are the top genes (that tell tissue apart no doubt)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )

png('plots/heatmap.png', width=8, height=8, res=300, units='in')
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", trace="none", dendrogram="column", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), cexRow=0.5, cexCol=0.5, margins=c(12,8),srtCol=45)
dev.off()
# volcano plot this is an exercise

#res <- res[order(res$padj),]

#deg=subset(res,padj<0.05)

#write.table(as.data.frame(deg), file='DEGs.txt')

#Lifted this function off the internet to turn ENSEMBL IDs to gene symbols
library(org.Mm.eg.db)
convertIDs <- function( ids, fromKey, toKey, db, ifMultiple=c( "putNA", "useFirst" ) ) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select( 
  db, keys=ids, keytype=fromKey, columns=c(fromKey,toKey) ) )
  if( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]   
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ] }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}

heatmap_matrix <- assay(rld)[topVarGenes,]
heatmap_matrix$gene_symbol <- convertIDs(row.names(heatmap_matrix), "ENSEMBL","SYMBOL", org.Mm.eg.db)

#Make volcano plot
library(EnhancedVolcano)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')