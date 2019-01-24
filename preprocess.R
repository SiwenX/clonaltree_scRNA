#
library(Matrix)
matrix_dir = "/Users/xsw/bio/Multiple_myeloma/data/scRNA/2/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv")
features.path <- paste0(matrix_dir, "genes.tsv")
matrix.path <- paste0(matrix_dir, "matrix.mtx")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1
#
csv<-read.csv("/Users/xsw/bio/Multiple_myeloma/data/scRNA/2/matrix.csv",header = T)
csv[,-1]->csv
colnames(csv) = barcode.names$V1
rownames(csv) = feature.names$V1
#remove non-expression gene
csv[-(which(rowSums(csv)==0)),]->csv2
#statistics the number of cells have expression
expression_cells<-c()
f<-function(x){length(which(x!=0))}
apply(csv2,1,f)->expression_cells

#在%多少的cell中有molecule(6166*0.5 = 3083; 6166*0.6 = 3700; 6166*0.7 = 4317; 6166*0.8 = 4933; 6166*0.9 = 5550)
aa<-c()
for(i in 1:nrow(csv2)){
  if(length(which(csv2[i,]!=0))>=5550){c(aa,i)->aa}
}
csv2[aa,]->csv3
count1<-c()
count10<-c()
count20<-c()
count30<-c()
for(i in 1:nrow(csv3)){
  length(which(csv3[i,]!=0))->expression_cells
  if(sum(csv3[i,])/expression_cells >= 1){c(count1,i)->count1}
  if(sum(csv3[i,])/expression_cells >= 10){c(count10,i)->count10}
  if(sum(csv3[i,])/expression_cells >= 20){c(count20,i)->count20}
  if(sum(csv3[i,])/expression_cells >= 30){c(count30,i)->count30}
} 
rownames(csv3)->ROW; t(colnames(csv3))->COL;
write.table(csv3,"/Users/xsw/bio/Multiple_myeloma/data/scRNA/2/90%cells.csv",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(ROW,"/Users/xsw/bio/Multiple_myeloma/data/scRNA/2/rowname.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(COL,"/Users/xsw/bio/Multiple_myeloma/data/scRNA/2/colname.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
csv3[count10,]->csv4
rownames(csv4)->ROW; t(colnames(csv4))->COL;
write.table(csv4,"/Users/xsw/bio/Multiple_myeloma/data/scRNA/2/90%cellsUMI10.csv",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(ROW,"/Users/xsw/bio/Multiple_myeloma/data/scRNA/2/rowname.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(COL,"/Users/xsw/bio/Multiple_myeloma/data/scRNA/2/colname.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
##match gene ID with name
a1<-read.table("/Users/xsw/bio/Multiple_myeloma/data/scRNA/2/genes.tsv")
a2<-read.table("/Users/xsw/bio/Multiple_myeloma/data/scRNA/2/rowname.txt")
bb<-c()
for(i in 1:nrow(a2)){
  as.character(a1[which(as.character(a1$V1) == as.character(a2[i,1])),2])->b
  rbind(bb,b)->bb
}
write.table(bb,"/Users/xsw/bio/Multiple_myeloma/data/scRNA/2/temp.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
#match gene name with loci
genes<-read.table("/Users/xsw/bio/Multiple_myeloma/data/scRNA/2/rowname.txt")
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
aa<-c()
for(i in 1:nrow(genes)){
getBM(c("hgnc_symbol","chromosome_name","start_position","end_position","ensembl_gene_id"), "ensembl_gene_id", genes$V1[66], mart) -> a
  rbind(aa, a) -> aa
}
write.table(aa,"/Users/xsw/bio/Multiple_myeloma/data/scRNA/2/temp.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
#

