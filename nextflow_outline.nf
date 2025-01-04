process hisat {
conda "environment.yml"
publishDir params.outdir, mode: 'copy'

input:
path fastq_input
path reference

output:
path "${fastq_input.simpleName}.bam"

script:
"""
hisat -q -x ${reference} -U ${fastq_input} | samtools sort -o ${fastq_input.simpleName}.bam
"""
}
process featurecounts {
conda "environment.yml"
publishDir params.outdir, mode: 'copy'

input:
path bams
path gtf

output:
path "counts.out"

script:
"""
featureCounts -a ${gtf} -o counts.out ${bams}
"""
}

process deseq2 {
conda "environment.yml"
publishDir params.outdir, mode: 'copy'

input:
path counts
path metadata

output:
path "deseq_results.csv"

script:
"""
#!/usr/bin/env Rscript
library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(dplyr)

# >>>>>>>>> Prepare counts data
counts_data = read.csv("${counts}",row.names=1)
head(counts_data)
colnames(counts_data)

# >>>>>>>>> sample info
sample_metadata = read.csv("${metadata}",row.names=1)
rownames(sample_metadata)
sample_metadata$experiment

level_factor = factor(sample_metadata$experiment)
sample_metadata$experiment = level_factor

# >>>>>>>>> check for formats ! 
all(colnames(counts_data) %in% rownames(sample_metadata)) # all should be present
all(colnames(counts_data) == rownames(sample_metadata)) # order matters

# >>>>>>>>> Dseq object
dds = DESeqDataSetFromMatrix(countData = counts_data,
                             colData = sample_metadata,
                             design= ~experiment)

# >>>>>>>>> filter rows with alteast gene count 10
keep_genes = rowSums(counts(dds)) >= 10
dds = dds[keep_genes,]

# >>>>>>>>> define levels! standard/variable
dds$experiment = relevel(dds$experiment, ref = "control")

# >>>>>>>>> Run DESeq!! 
dds = DESeq(dds)
res = results(dds)
# summary(res)
# names(res)
res$id = row.names(res)
res = subset(res, select=c("id","padj","pvalue","lfcSE","stat","log2FoldChange","baseMean"))
res = na.omit(res)
write.csv(as.data.frame(res), file="deseq_results.csv", row.names=F)

# >>>>>>>>> DEG's
deg = subset(res, padj<0.05 & abs(log2FoldChange) >= 1)
deg = deg[order(deg$padj),]
deg = as.data.frame(deg)
deg$geneSymbol = mapIds(org.Hs.eg.db, keys = rownames(deg), keytype = "ENSEMBL", column = "SYMBOL")
write.csv(deg, file="top_DEGs.csv", row.names=F)

"""
}

process gct_file {
conda "environment.yml"
publishDir params.outdir, mode: 'copy'

input:
path deg_output

output:
path "final_dseq_data.gct"

script:
"""
#!/usr/bin/env python

import pandas as pd
raw_counts = pd.read_csv('${deg_output}', index_col=0)
genes, samples = raw_counts.shape
header = f"#1.2\n{genes}\t{samples}\n"
description_col=['metadata unavailable' for i in range(len(raw_counts)]
raw_counts.insert(0, 'Description', description_col)
gct_data = raw_counts.reset_index()
gct_data.columns = ['id'] + list(gct_data.columns[1:])
with open('final_deq_data.gct', 'w') as f:
    f.write(header)
    gct_data.to_csv(f, sep=',', index=False)
"""
}

########## main.nf
nextflow.enable.dsl = 2

include { hisat } from './modules/hisat'
include { featurecounts } from './modules/featurecounts'
include { hisat } from './modules/hisat'
include { gct_file } from './modules/gct_file'

params.reads = "path/*_{1,2}.fastq.gz"
params.reference = "path/genome_indices/"
params.gtf = "path/annotation.gtf"
params.metadata = "path/metadata.csv"
# params.bam = "path/*.bam"
params.outdir = "path/output_results"



workflow {
read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
aligned_reads = hisat(read_pairs_ch,params.reference)
# read_bam = Channel.fromPath(${aligned_reads.getParent()}, checkIfExists: true)
# if( Channel.of(read_pairs_ch).count().view() == number_of_samples ) {
#     counts = featurecounts(params.bam,params.gtf)
# }
counts = featurecounts(aligned_reads.getParent(),params.gtf)
deseq_result = deseq2(counts, params.metadata)
gct_file(deseq_result)
}




