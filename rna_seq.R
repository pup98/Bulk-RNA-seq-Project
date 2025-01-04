
library(DESeq2)
library(tidyverse)
library(airway)
library(EnhancedVolcano)
library(pheatmap)
library(Rsubread)
library(org.Hs.eg.db)
library(dplyr)

# Prepare counts data
counts_data = read.csv('counts_data_v1.csv',row.names=1)
head(counts_data)
colnames(counts_data)
# sample info
sample_metadata = read.csv('sample_metadata.csv',row.names=1)
sample_metadata
rownames(sample_metadata)
sample_metadata$experiment

level_factor = factor(sample_metadata$experiment)
level_factor
sample_metadata$experiment = level_factor
sample_metadata$experiment

# check for sample formats ! 
all(colnames(counts_data) %in% rownames(sample_metadata)) # all should be present
all(colnames(counts_data) == rownames(sample_metadata)) # order matters

# >>> Dseq object <<
dds = DESeqDataSetFromMatrix(countData = counts_data,
                             colData = sample_metadata,
                             design= ~experiment)

dds
# filter rows with alteast gene count 10
keep_genes = rowSums(counts(dds)) >= 10
keep_genes
dds = dds[keep_genes,]
dds

# define levels! standard/variable
dds$experiment = relevel(dds$experiment, ref = "control")


# Run DESeq!! 
dds = DESeq(dds)
res = results(dds)
res
summary(res)
names(res)
res$id = row.names(res)
names(res)
res = subset(res, select=c("id","padj","pvalue","lfcSE","stat","log2FoldChange","baseMean"))
names(res)
res = na.omit(res)

#ordered = res[order(res$padj),]
#ordered


#res = results(dds, alpha=0.01)
#summary(res)

resultsNames(dds) # gives level of comparison bw samples

write.csv(as.data.frame(res), file="deseq2_results.csv", row.names=F)

# >>>>>>>>> DEG's have pval < defalut 0.05 

deg = subset(res, padj<0.05 & abs(log2FoldChange) >= 1)
deg = deg[order(deg$padj),]
deg
dim(deg)
dim(res)
deg = as.data.frame(deg)
deg$geneSymbol = mapIds(org.Hs.eg.db, keys = rownames(deg), keytype = "ENSEMBL", column = "SYMBOL")
deg
write.csv(deg, file="deg_results.csv", row.names=F)

# >>>>>>>>> plotting ... 

plotMA(res, main='MA Plot')
#plotCounts(dds, gene="CD38", intgroup="cancer")

hist(res$padj, breaks = seq(0,1,length=20), col = "red", border = "white", main="Frequencies of padj-values")

# volcano
old.pal = palette(c("#00BFFF","#FF3030"))
par(mar=c(4,4,2,1), cex.main=1.5)

plot(res$log2FoldChange, -log10(res$padj),xlab="logfc",ylab="-log10(padj)",pch=20, cex=0.5)


with(subset(res, padj <0.05 & abs(log2FoldChange) >=1),
     points(log2FoldChange, -log10(padj), pch=20, col=(sign(log2FoldChange)+3)/2, cex=1))
legend("bottomright", title=paste("padj<",0.05,sep=""),
  legend=c("down","up"), pch=20,col=1:2)



EnhancedVolcano(res,
                lab = res$id,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Volcano Plot')


## >>>>>> pca
vsd = vst(dds,blind=FALSE)
plotPCA(vsd, intgroup=c('experiment'))


## heatmap
normalized_counts = counts(dds, normalized=T)
head(normalized_counts)
transformed_counts = log2(normalized_counts+1)
top_degs = row.names(deg[1:6,])
top_degs = transformed_counts[top_degs,]
pheatmap(top_degs, cluster_rows = FALSE, cluster_cols = FALSE, main="Heatmap")
