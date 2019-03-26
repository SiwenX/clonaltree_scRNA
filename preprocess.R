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
#genes in csv4 are the marker-gene
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
#
a1<-read.table("/Users/xsw/bio/scRNA/data/2/FTL",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
FTL<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(FTL,a)->FTL
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/SSR4",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
SSR4<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(SSR4,a)->SSR4
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/RPL34",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
RPL34<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(RPL34,a)->RPL34
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/IGLV2-18",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
IGLV2_18<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(IGLV2_18,a)->IGLV2_18
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/RPS27A",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
RPS27A<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(RPS27A,a)->RPS27A
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/RPL9",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
RPL9<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(RPL9,a)->RPL9
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/RPS28",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
RPS28<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(RPS28,a)->RPS28
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/RPS15",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
RPS15<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(RPS15,a)->RPS15
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/CYBA",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
CYBA<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(CYBA,a)->CYBA
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/RPL18A",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
RPL18A<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(RPL18A,a)->RPL18A
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/RPS19",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
RPS19<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(RPS19,a)->RPS19
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/RPS27",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
RPS27<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(RPS27,a)->RPS27
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/HLA-B",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
HLA_B<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(HLA_B,a)->HLA_B
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/RPS18",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
RPS18<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(RPS18,a)->RPS18
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/RPL13A",fill = T,header = F)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
RPL13A<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(RPL13A,a)->RPL13A
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/RPL8",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
RPL8<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(RPL8,a)->RPL8
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/MZB1",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
MZB1<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(MZB1,a)->MZB1
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/RPL3",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
RPL3<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(RPL3,a)->RPL3
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/EEF1A1",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
EEF1A1<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(EEF1A1,a)->EEF1A1
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/IGHA1",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
IGHA1<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(IGHA1,a)->IGHA1
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/RPL13",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
grep(UR,UMI)
RPL13<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(RPL13,a)->RPL13
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/RPLP1",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
RPLP1<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(RPLP1,a)->RPLP1
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/RPL10",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
RPL10<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(RPL10,a)->RPL10
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/HLA-C",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
HLA_C<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(HLA_C,a)->HLA_C
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/MALAT1",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
MALAT1<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(MALAT1,a)->MALAT1
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/IGHA2",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
IGHA2<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(IGHA2,a)->IGHA2
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/B2M",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
B2M<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(B2M,a)->B2M
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/JCHAIN",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
JCHAIN<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(JCHAIN,a)->JCHAIN
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/IGLC2",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
IGLC2<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(IGLC2,a)->IGLC2
}
a1<-read.table("/Users/xsw/bio/scRNA/data/2/RPL39",fill = T)
unique(a1$V25)->UMI
UMI[grep('UR',UMI)]->UMI
RPL39<-c()
for(i in 1:length(UMI)){
  length(which(a1$V25 == UMI[i]))->a
  c(RPL39,a)->RPL39
}

par(mfcol=c(2,3))
hist(IGLC2, plot = T)
hist(IGHA1, plot = T)
hist(JCHAIN, plot = T)
hist(B2M, plot = T)
hist(IGHA2, plot = T)
hist(MALAT1, plot = T)
par(mfcol=c(2,3))
hist(HLA_C, plot = T)
hist(RPL10, plot = T)
hist(RPLP1, plot = T)
hist(RPL13, plot = T)
hist(EEF1A1, plot = T)
hist(RPL3, plot = T)
par(mfcol=c(2,3))
hist(MZB1, plot = T)
hist(RPL8, plot = T)
hist(RPL13A, plot = T)
hist(RPS18, plot = T)
hist(SSR4, plot = T)
hist(FTL, plot = T)
par(mfcol=c(2,3))
hist(HLA_B, plot = T)
hist(RPS27, plot = T)
hist(RPS19, plot = T)
hist(RPL18A, plot = T)
hist(CYBA, plot = T)
hist(RPS15, plot = T)
par(mfcol=c(2,3))
hist(RPS28, plot = T)
hist(RPS27A, plot = T)
hist(RPL3, plot = T)
hist(RPL34, plot = T)
hist(IGLV2_18, plot = T)
hist(RPL39, plot = T)

RPL13->a
length(which(a==1))->a1
length(which(a==2))->a2
length(which(a==3))->a3
length(which(a==4))->a4
length(which(a==5))->a5
length(which(a==6))->a6
length(which(a==7))->a7
length(which(a==8))->a8
length(which(a==9))->a9
length(which(a==10))->a10
length(which(a>10))->a11
c(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11)->aa
