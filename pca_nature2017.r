install.packages('stats')
library('stats')
library('matrixStats')

#increment function
new<-function(x){
  x+1
}

#read csv file containing RPKM values (first 2 columns for gene names and ensembl ids) 
y<-read.csv("", header=TRUE)

#set boundary of additional data points (test set vs training set)
end=87
start=89

#read colour file (column 1=sample names; column 2=colour)
SampleCol<-read.csv("", header=TRUE)

#read colour file for additional data points (column 1=sample names; column 2=colour)
CrisprCol<-read.csv("")

#read legend file (column 1=colour; column 2=sample)
Legend<-read.csv("")

#prepare data
colnum<-ncol(y)
colnum2=colnum-2
y2<-as.matrix(y[,3:colnum])
rownames(y2)<-paste(y[,1],y[,2])
colour<-as.vector(SampleCol[,2])
crisprcolour<-as.vector(CrisprCol[,2])
shape<-as.vector(CrisprCol[,3])
legendlab<-as.vector(Legend[,2])
legendcol<-as.vector(Legend[,1])

#filter genes with rpkm more than 5 and take top 8000 most variably expressed genes
y2 <- cbind(y2, yr.above = rowSums(y2 > 5))
y2<-y2[y2[, "yr.above"] >10,]
y3<-y2[,1:colnum2-1]
y4<-as.matrix(apply(y3,2,new))
y4<-log2(y4)
y4 <- cbind(y4, yr.above = rowMaxs(y4)-rowMins(y4))
y5<-y4[order(y4[,colnum-2]),]
y6<-tail(y5, n = 8000)
crispr<-y6[,1:end]
yanonly<-y6[,start:colnum2-1]

#transpose main points
tscy<-t(yanonly)
#transpose additional points
crispr_tscy<-t(crispr)

#PCA
yanpca<-prcomp(tscy,center=T,scale=T)

#calculate PC percentages
pc1 <- round(yanpca$sdev[1]^2/sum(yanpca$sdev^2)*100,2)
pc2 <- round(yanpca$sdev[2]^2/sum(yanpca$sdev^2)*100,2)
pc3 <- round(yanpca$sdev[3]^2/sum(yanpca$sdev^2)*100,2)

#create labels and project additional points onto PCA
plabs=yanpca$x[,0]
plabs.df<-data.frame(rownames(plabs))
pblabs2<-as.vector(t(plabs.df))
predict_yanpca<-predict(yanpca, crispr_tscy)
plabs_pred=predict_yanpca[,0]
pred_plabs.df<-data.frame(rownames(plabs_pred))
pred_pblabs2<-as.vector(t(pred_plabs.df))

#plot PCA and save to file
pdf("test.pdf")
plot(yanpca$x[,2],yanpca$x[,3],xlab=paste('PC2', pc2, '%', sep=" "),ylab=paste('PC3', pc3, '%', sep=" "),col=colour,pch=19)
points(predict_yanpca[,2],predict_yanpca[,3], pch=shape,col=crisprcolour, text(predict_yanpca[,2],y=predict_yanpca[,3],pred_pblabs2,cex=0.25,adj=1.3))
legend(34,112,col=legendcol,legend=legendlab,cex=0.5,lwd=5)
dev.off()