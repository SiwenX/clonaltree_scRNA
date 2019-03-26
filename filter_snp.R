library(stringr)
tryCatch({snp<-read.table("/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/SNP.bed",header=FALSE)} ,error = function(e)
{ matrix(0,1,3)->aa
  cat(aa, file = "/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/SNP.bed")
})
a<-c()
for(i in 1:nrow(snp)){
  aa<-c()
  if(length(str_split(snp[i,3],'')[[1]])!=1){i->aa}
  c(a,aa)->a
}
if(!is.null(a)){
snp[-a,]->filter_snp
write.table(filter_snp,"/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/SNP.bed",quote=FALSE,row.names=FALSE,col.names=F)
}else{
    write.table(snp,"/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/SNP.bed",quote=FALSE,row.names=FALSE,col.names=F)
}
