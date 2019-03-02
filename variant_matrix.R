arg <- commandArgs(T)
library("GenomicRanges")

cell<-read.table("/Users/xsw/bio/scRNA/data/2/6166cell/vcf/temp.txt",header=FALSE)
Gene<-read.table("/Users/xsw/bio/scRNA/data/2/gene_region.bed",header=F)
cell$V2->cell$V3
cell[,c(1,2,3)]->SNP
for(i in 1:100){
  Gene[i,1]->GeneName
  Gene[i,c(2,3,4)]->gene
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
    gene_leng[which(colnames(variant_matrix)==temp2[1])]<-1
  }

  rownames(variant_matrix)[i]<-as.character(arg[1])
  cat(gene_leng, file = paste('/Users/xsw/bio/scRNA/data/2/6166cell/vcf/',GeneName,'_variant_matrix.txt',sep = ''), append = TRUE)
  cat("\n", file = paste('/Users/xsw/bio/scRNA/data/2/6166cell/vcf/',GeneName,'_variant_matrix.txt',sep = ''), append = TRUE)
}
