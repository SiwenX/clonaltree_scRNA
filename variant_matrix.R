arg <- commandArgs(T)
library("GenomicRanges")
library(stringr)
gene<-read.table("/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/gene_region.bed",header=FALSE)
gene[1,1]->GeneName
as.numeric(gene[4])-as.numeric(gene[3])+1 -> length
rep(0,length) -> aaa
tryCatch({cell<-read.table("/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/temp.txt",header=F)} ,error = function(e)
{
  cat(aaa, file = paste('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/',GeneName,'_variant_matrix2.txt',sep = ''), append = TRUE)
  cat("\n", file = paste('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/',GeneName,'_variant_matrix2.txt',sep = ''), append = TRUE)
})
cell$V2->cell$V3
cell[,c(1,2,3)]->SNP

  gene[1,c(2,3,4)]->gene
  matrix(0, nrow = 1, ncol = (as.numeric(gene[3])-as.numeric(gene[2])+1))->variant_matrix
  colnames(variant_matrix)<-c(as.numeric(gene[2]):as.numeric(gene[3]))
  colnames(SNP)<-c("chr","start","end"); colnames(gene)<-c("chr","start","end");
  p1Range<-with(SNP,GRanges(seqnames=chr,ranges=IRanges(start=start,end=end),strand="*"));
  p2Range<-with(gene,GRanges(seqnames=chr,ranges=IRanges(start=start,end=end),strand="*"));
  intersect(p1Range,p2Range)->intersection; 
  as.data.frame(intersection)->intersection
  intersection[,1:3]->peak
  gene_leng<-rep(0,as.numeric(gene[3])-as.numeric(gene[2])+1)

  for(i in 1:nrow(peak)){ 
    cell[which(cell[,2]==peak[i,2]),]->temp
    temp$V2->temp2
    if(length(temp2)==0){break}
    str_split(temp$V10,':')[[1]][3]->temp3
    str_split(temp3,',')[[1]]->genotype
    if(as.numeric(genotype[1])!=0){
      gene_leng[which(colnames(variant_matrix)==temp2[1])]<-1}
    if(as.numeric(genotype[1])==0){
      gene_leng[which(colnames(variant_matrix)==temp2[1])]<-2}

  }
  
cat(gene_leng, file = paste('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/',GeneName,'_variant_matrix2.txt',sep = ''), append = TRUE)
cat("\n", file = paste('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/',GeneName,'_variant_matrix2.txt',sep = ''), append = TRUE)

