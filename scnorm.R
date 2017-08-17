#! /usr/bin/R

#===============================================================================
# SCnorm normalisation of single-cell transcript counts
# Copyright (C) 2017 Jens Kleinjung 
#===============================================================================

library(SCnorm);

## input data
counts.raw = read.table("genali.quality", header = TRUE);
counts = counts.raw[ , -1];
rownames(counts) = counts.raw[ , "gene_ID"];

## trivial design matrix
Conditions = rep(c(1), each = dim(counts)[2]);

## plot inital evaluation
checkCountDepth(Data = counts, Conditions = Conditions,
	OutputName = "checkData", FilterCellProportion = .1, NCores = 8);

## normalisation
DataNorm = SCnorm(counts, Conditions, OutputName = "normData",
	SavePDF = TRUE, FilterCellNum = 5, NCores = 8);

## normalised data
head(DataNorm$NormalizedData);

write.table(DataNorm$NormalizedData, file = "genali.scnorm");

