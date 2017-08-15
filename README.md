# Index of files
* genali.merge:
	raw transcript count data of all samples after alignment

* genali.quality:
	transcript count data after quality filter

* genali.scnorm:
	transcript count data after SCnorm normalisation

* RUV.bln.counts:
	transcript count data after RUVSeq between-lane normalisation
	using the 'full' quantile regression

* samples_meta.dat:
	sample meta data

* PCA_gene_list.dat:
	gene list to subset transcript count data for PCA

* PCA_KO_1.pdf:
	PCA of 'RUV.bln.counts' after subsetting to above gene list
	and log-transformation with per-gene centering

* PCA_KO_2.pdf:
	same plot labelled with sample identifier

* PCA_KO_3.pdf:
	same plot labelled with cell type

* PCA_KO_4.pdf:
	same plot labelled with sequencing depth (total transcript count per cell)
	expressed in rounded log2 format

* PCA_KO_5.pdf:
	same plot labelled with genome mapping rate in per-cent format

* program_versions.txt:
	list of programs and their versions used to process the data

