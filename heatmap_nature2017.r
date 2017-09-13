library('pheatmap')
library('RColorBrewer')

#increment function
new<-function(x){
  x+1
}

#insert file name here
fname<-""

#read rpkm table (2 columns for gene names and identintifier )
rpkm<-read.table(fname, header=TRUE, sep=",", check.names=FALSE)
rowlen=ncol(rpkm)

#set rownames
rpkm2<-as.matrix(rpkm[,3:rowlen])
rownames(rpkm2)<-rpkm[,2]

#convert to matrix and convert non-zero values
rpkm3<-as.matrix(apply(rpkm2,2,new))
#convert to log scale
rpkml2<-log2(rpkm3)

#generate heatmap
pheatmap(rpkml2, cellwidth=5, 
         cluster_rows=FALSE,
         cluster_cols=TRUE,
         scale='row',
         fontsize_col=6,
         fontsize_row=6,
         cellheight = 5,
         color = colorRampPalette(rev(brewer.pal(n = 4, name = "PiYG")))(4),
         breaks=c(-2,-1,1,2)
)