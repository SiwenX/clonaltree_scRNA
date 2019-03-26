Cell<-read.table("/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/barcodes.tsv",header=F)
Gene<-read.table("/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/gene_region.bed",header=F)
RM<-read.table('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/IGLC2_variant_matrix2.txt', header=F, col.names = c((as.numeric(Gene[3])):as.numeric(Gene[4])))
row.names(RM)<-as.character(Cell[,1])
write.table(RM,'/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/IGLC2_variant_matrix2.txt',quote=FALSE,row.names=T,col.names=T)


