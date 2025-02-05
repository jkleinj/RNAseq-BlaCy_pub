# Index of files

The _*_ is a placeholder for the acronyms 'bln' or 'nobln',
standing for 'between lane normalisation' or 'no between lane normalisation'.

* genali.merge
: raw transcript count data of all samples after alignment

* genali.quality
: transcript count data after quality filter

* genali.scnorm
: transcript count data after SCnorm normalisation

* RUV.bln.counts
: transcript count data after RUVSeq between-lane normalisation
	using the 'full' quantile regression

* SupplementaryTable2_samples_meta.[csv,dat]
: sample meta data containing the columns
"sample_ID", "cell_ID", "cell_type", "smp_ID", "CRISPR", "embryo", "library_ID", "alt_cell_ID", "read_counts", "mapping_rate". "read_counts" is the total number of trancscript counts per cell, also referred to as 'sequencing depth'.

* PCA_gene_list.dat
: gene list used for sub-setting transcript count data to compute PCA plot

* PCA_KO_*_1.pdf
: PCA of 'RUV.bln.counts' after subsetting to above gene list
	and log-transformation with per-gene centering

* PCA_KO_*_2.pdf
: same plot labelled with sample identifier

* PCA_KO_*_3.pdf
: same plot labelled with cell type

* PCA_KO_*_4.pdf
: same plot labelled with sequencing depth (total transcript count per cell)
	expressed in rounded log2 format

* PCA_KO_*_5.pdf
: same plot labelled with genome mapping rate in per-cent format

* program_versions.txt
: list of programs and their versions used to process the data

* deseq_KO_vs_uninjected_control_1.csv:
: gene list for subsetting

* RUVSeq.R
: between-lane normalisation using the RUVSeq package

* gene_name_ID.R
: program to match gene names and IDss

* qual.R
: quality filter for genes and cells (samples)

* scLVM_KO_*.R
: subsetting of genes to the padj<1e-2 part of deseq_KO_vs_uninjected_control_1.csv

* scnorm.R
: within-sample normalisation using the SCnorm package

* heatmap_nature2017.r
: heatmap script

* pca_nature2017.r
: PCA script
