#! /usr/bin/R

#===============================================================================
# RUVseq normalisation of sc transcript counts
# Copyright (C) 2017 Jens Kleinjung
#===============================================================================

library("RUVSeq");
library("RColorBrewer");
library("ggfortify");
library("org.Hs.eg.db");
library("rgl");

source("gene_name_ID.R");

#_______________________________________________________________________________
## transcript count data (genes x samples)
counts = read.table("genali.scnorm", header = TRUE);
gene_names = select(org.Hs.eg.db, keys = rownames(counts), columns = c("SYMBOL"), keytype="ENSEMBL");
gene_names.unique = gene_names[!duplicated(gene_names$ENSEMBL), ];

## RUVSeq object
cpe = newSeqExpressionSet(as.matrix(round(counts)),
	phenoData = data.frame(colnames(counts), row.names = colnames(counts)));

#colours = brewer.pal(8, "Set2");
pdf("plot.RLE.raw.pdf");
plotRLE(cpe, outline = FALSE, ylim = c(-4, 4));
dev.off();

#_______________________________________________________________________________
## between lane normalisation
cpe = betweenLaneNormalization(cpe, which = "full");

pdf("plot.RLE.norm.pdf");
plotRLE(cpe, outline = FALSE, ylim = c(-3, 3));
dev.off();

pdf("plot.PCA.norm.pdf");
plotPCA(cpe, cex = 1.2);
dev.off();

write.table(normCounts(cpe), file = "RUV.bln.counts");

#===============================================================================
 
