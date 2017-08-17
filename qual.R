#! /usr/bin/R

#===============================================================================
# quality check of genes and cells 
# Copyright (C) 2017 Jens Kleinjung 
#===============================================================================

#_______________________________________________________________________________
## input
## transcript counts; genes x cells
dat.raw = read.table("genali.merge", header = TRUE);

dat = as.matrix(dat.raw[ , -1]);
dimnames(dat) = list(dat.raw[ , 1], colnames(dat.raw[ , -1]));

#_______________________________________________________________________________
## valid genes must have at least 5 counts in 2 samples
valid_genes = apply(dat, 1, function (x) { length(x[x > 5]) >= 5 });
write.table(valid_genes, file = "valid_genes.dat", quote = FALSE, col.names = FALSE);
vg = table(valid_genes);
write.table(vg, file = "valid_genes_table.dat");

#_______________________________________________________________________________
## valid cells (samples)
## vector of total transcript counts per cell
cellsum.v = apply(dat, 2, sum);
write.table(cellsum.v, file = "cellsums.dat", quote = FALSE, col.names = FALSE);
## valid cells: 
valid_cells = cellsum.v > .5e5 & cellsum.v < .5e8;
write.table(valid_cells, file = "valid_cells.dat", quote = FALSE, col.names = FALSE);
vc = table(valid_cells);
write.table(vc, file = "valid_cells_table.dat");
pdf("logcellcounts.pdf");
plot(log10(cellsum.v));
dev.off();

pdf("transcriptcounts_all.pdf");
boxplot(cellsum.v, xlab = c("all cell samples"), ylab = c("transcript count"));
dev.off();

pdf("transcriptcounts_selected.pdf");
boxplot(cellsum.v[valid_cells], xlab = c("selected cell samples"), ylab = c("transcript count"));
dev.off();

sink("invalid_cells.dat");
print(names(cellsum.v)[!valid_cells]);
sink();

#_______________________________________________________________________________
## output
write.table(dat[names(valid_genes)[valid_genes], names(valid_cells)[valid_cells]],
	file = "genali.quality", quote = FALSE);

