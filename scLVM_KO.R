#! /usr/bin/R

#____________________________________________________________________________
library("genefilter");
library("statmod");
library("ggplot2");
library("DESeq2");
library("scLVM");
library("org.Hs.eg.db");
library("Rtsne");
library("cluster");
library("rgl");
library("MASS");

#_______________________________________________________________________________
## gene names and IDs
source("gene_name_ID.R");

#_______________________________________________________________________________
## load data
dataHs.raw = read.table("RUV.bln.counts", header = TRUE);

gene_names = AnnotationDbi::select(org.Hs.eg.db, keys = rownames(dataHs.raw), columns = c("SYMBOL"), keytype="ENSEMBL");
## remove multiple mappings
gene_names.unique = gene_names[!duplicated(gene_names$ENSEMBL), ];

#_______________________________________________________________________________
## gene markers
markers_KO = read.table("deseq_KO_vs_uninjected_control_1.csv", header = TRUE);

#_______________________________________________________________________________
## cell meta data
samples = read.table("samples.dat", header = TRUE);
## select samples in count matrix
dataHs = dataHs.raw[ , which(colnames(dataHs.raw) %in% samples$sample_ID)];

#_______________________________________________________________________________
## get technical noise

## If no spike-ins are available, we can also use the endogenous read counts for
##   fitting the mean-CV2 relation using a log-linear fit in the log-space.
techNoiseLogFit = fitTechnicalNoise(dataHs, fit_type = 'log', use_ERCC = FALSE, plot = TRUE);

#_______________________________________________________________________________
## variable genes
is_het = getVariableGenes(dataHs, techNoiseLogFit$fit, method = "fit", threshold = 0.1, fit_type = "log");
table(is_het);

#_______________________________________________________________________________
## normalised pseudo-log-transformed read counts
Y = apply(asinh(dataHs) / 2, 1, function(x) scale(x, center = TRUE, scale = FALSE));
rownames(Y) = colnames(dataHs);

## variable genes
genes_het_bool = as.vector(is_het);
## gene IDs
geneID = rownames(dataHs);
## technical noise
tech_noise = as.vector(techNoiseLogFit$techNoiseLog);

#_______________________________________________________________________________
## subsetting on marker genes
#_______________________________________________________________________________
## KO genes 
markers_KO_in = markers_KO[as.character(markers_KO$Ensembl) %in% colnames(Y), ];

## 'padj' < 1e-2
markers_KO_ss = subset(markers_KO_in, padj < 1e-2); 
idx_KO_in = as.character(markers_KO_ss[, "Ensembl"]);

Yko = Y[ , idx_KO_in];

#_______________________________________________________________________________
##  colours
cell.types = unlist(sapply(rownames(Y), function(x) {
                samples[samples$sample_ID == x, "cell_type"]
			}));

crispr.types = unlist(sapply(rownames(Y), function(x) {
                samples[samples$sample_ID == x, "CRISPR"]
			}));

embryo.types = unlist(sapply(rownames(Y), function(x) {
                samples[samples$sample_ID == x, "embryo"]
			}));

sample.types = unlist(sapply(rownames(Y), function(x) {
                samples[samples$sample_ID == x, "smp_ID"]
			}));

cellID.types = unlist(sapply(rownames(Y), function(x) {
                samples[samples$sample_ID == x, "cell_ID"]
			}));

#colfunc = colorRampPalette(c("darkblue", "orange", "yellow", "blue", "lightblue", "black"));
colfunc.cell_types = colorRampPalette(c("lightblue", "blue", "orange", "black"));
colours.cell_types = colfunc.cell_types(4)[cell.types];

colfunc.crispr_types = colorRampPalette(c("blue", "orange"));
colours.crispr_types = colfunc.crispr_types(2)[crispr.types];

#_______________________________________________________________________________
## PCA
Y.pca = prcomp(Yko);
variance = sapply(Y.pca$sdev, function(x) { x^2 / sum(Y.pca$sdev^2) });
write.table(variance, "PCA_variance.dat");

#_______________________________________________________________________________
## plot1
pdf("PCA_KO_1.pdf");
ggplot(as.data.frame(cbind(Y.pca$x[ , 1], Y.pca$x[ , 2])),
	   aes(x = V1, y = V2,
			colour = cell.types, shape = embryo.types)) +
	geom_point(size = 4) +
	scale_colour_manual(values = c("yellow2", "orange", "blue", "grey40")) +
	scale_shape_manual(values = c(0, 1, 2, 5, 15, 16, 17, 18)) +
	xlab("PC1 (21% explained var.)") + ylab("PC2 (13% explained var.)") + 
	theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14)) +
	theme(legend.text = element_text(size = 12)) +
	theme(legend.title = element_text(size = 12)) +
	guides(shape = guide_legend(override.aes = list(size = 3))) +
	guides(colour = guide_legend(override.aes = list(size = 3)));
dev.off();

#_______________________________________________________________________________
## plot2, with sample labels
pdf("PCA_KO_2.pdf");
ggplot(as.data.frame(cbind(Y.pca$x[ , 1], Y.pca$x[ , 2])),
	   aes(x = V1, y = V2,
			colour = cell.types, shape = embryo.types, label = sample.types)) +
	geom_point(size = 4) +
	scale_colour_manual(values = c("yellow2", "orange", "blue", "grey40")) +
	scale_shape_manual(values = c(0, 1, 2, 5, 15, 16, 17, 18)) +
	geom_text(size = 3, hjust = -0.4, vjust = -0.4, colour = "black") +
	xlab("PC1 (21% explained var.)") + ylab("PC2 (13% explained var.)") + 
	theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14)) +
	theme(legend.text = element_text(size = 12)) +
	theme(legend.title = element_text(size = 12)) +
	guides(shape = guide_legend(override.aes = list(size = 3))) +
	guides(colour = guide_legend(override.aes = list(size = 3)));
dev.off();

#_______________________________________________________________________________
## plot3, with cellID labels
pdf("PCA_KO_3.pdf");
ggplot(as.data.frame(cbind(Y.pca$x[ , 1], Y.pca$x[ , 2])),
	   aes(x = V1, y = V2,
			colour = cell.types, shape = embryo.types, label = cellID.types)) +
	geom_point(size = 4) +
	scale_colour_manual(values = c("yellow2", "orange", "blue", "grey40")) +
	scale_shape_manual(values = c(0, 1, 2, 5, 15, 16, 17, 18)) +
	geom_text(size = 3, hjust = -0.2, vjust = -0.2, colour = "black") +
	xlab("PC1 (21% explained var.)") + ylab("PC2 (13% explained var.)") + 
	theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14)) +
	theme(legend.text = element_text(size = 12)) +
	theme(legend.title = element_text(size = 12)) +
	guides(shape = guide_legend(override.aes = list(size = 3))) +
	guides(colour = guide_legend(override.aes = list(size = 3)));
dev.off();

#_______________________________________________________________________________
# add read-depth (cellsums) and mapping rate
cellsums = read.table("cellsums.dat", header = TRUE);
mapping_rate = read.table("mapping_rate.dat", header = TRUE);
samples_rd = merge(samples, cellsums, by = c("sample_ID"), all = FALSE);
samples_rd_mr = merge(samples_rd, mapping_rate, by = c("sample_ID"), all = FALSE);
write.table(samples_rd_mr, "samples_stats.dat");

readCounts.types = as.character(floor(log2(unlist(sapply(rownames(Y), function(x) {
                samples_rd_mr[samples_rd_mr$sample_ID == x, "read_counts"]
			})))));

mappingRates.types = unlist(sapply(rownames(Y), function(x) {
                samples_rd_mr[samples_rd_mr$sample_ID == x, "mapping_rate"]
			}));

#_______________________________________________________________________________
## plot4, with read-depth labels
pdf("PCA_KO_4.pdf");
ggplot(as.data.frame(cbind(Y.pca$x[ , 1], Y.pca$x[ , 2])),
	   aes(x = V1, y = V2,
			colour = cell.types, shape = embryo.types, label = readCounts.types)) +
	geom_point(size = 4) +
	scale_colour_manual(values = c("yellow2", "orange", "blue", "grey40")) +
	scale_shape_manual(values = c(0, 1, 2, 5, 15, 16, 17, 18)) +
	geom_text(size = 4, hjust = -0.2, vjust = -0.2, colour = "black") +
	xlab("PC1 (21% explained var.)") + ylab("PC2 (13% explained var.)") + 
	theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14)) +
	theme(legend.text = element_text(size = 12)) +
	theme(legend.title = element_text(size = 12)) +
	guides(shape = guide_legend(override.aes = list(size = 3))) +
	guides(colour = guide_legend(override.aes = list(size = 3)));
dev.off();

#_______________________________________________________________________________
## plot5, with alignment rate labels
pdf("PCA_KO_5.pdf");
ggplot(as.data.frame(cbind(Y.pca$x[ , 1], Y.pca$x[ , 2])),
	   aes(x = V1, y = V2,
			colour = cell.types, shape = embryo.types, label = mappingRates.types)) +
	geom_point(size = 4) +
	scale_colour_manual(values = c("yellow2", "orange", "blue", "grey40")) +
	scale_shape_manual(values = c(0, 1, 2, 5, 15, 16, 17, 18)) +
	geom_text(size = 4, hjust = -0.2, vjust = -0.2, colour = "black") +
	xlab("PC1 (21% explained var.)") + ylab("PC2 (13% explained var.)") + 
	theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14)) +
	theme(legend.text = element_text(size = 12)) +
	theme(legend.title = element_text(size = 12)) +
	guides(shape = guide_legend(override.aes = list(size = 3))) +
	guides(colour = guide_legend(override.aes = list(size = 3)));
dev.off();

#_______________________________________________________________________________
## 3D PCA
## colouring by cell type
plot3d(Y.pca$x[ , 1:3], col = colours.cell_types, size = 8);
pp = par3d(no.readonly = TRUE);
rgl.postscript("tSNE_KO_celltype.pdf", "pdf");

## colouring by CRISPR type
plot3d(Y.pca$x[ , 1:3], col = colours.crispr_types, size = 8);
pp = par3d(no.readonly = TRUE);
rgl.postscript("tSNE_KO_crisprtype.pdf", "pdf");

#_______________________________________________________________________________
## 3D RTSNE
set.seed(2^4);
tsne = Rtsne(Yko, dims = 3, perplexity = 10, max_iter = 10000, check_duplicates = FALSE); 

## colouring by cell type
plot3d(tsne$Y[ , 1:3], col = colours.cell_types, size = 6);
rgl.postscript("tSNE_KO_celltype.pdf", "pdf")

## colouring by CRISPR type
plot3d(tsne$Y[ , 1:3], col = colours.crispr_types, size = 6);
## save view
pp = par3d(no.readonly = TRUE);
rgl.postscript("tSNE_two_groups_crisprtype.pdf", "pdf")

