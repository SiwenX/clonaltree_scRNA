arg <- commandArgs(T)
arg[1] -> gene
arg[2] -> cell
library(stringr)
tryCatch({snp<-read.table(paste0('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/', gene, '/cell_files/', cell, '/SNP.bed'), header=FALSE)
a<-c()
for(i in 1:nrow(snp)){
  aa<-c()
  if(length(str_split(snp[i,3],'')[[1]])!=1){i->aa}
  c(a,aa)->a
}
if(!is.null(a)){
snp[-a,]->filter_snp
write.table(filter_snp, paste0('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/', gene, '/cell_files/', cell, '/SNP.bed'), quote=FALSE, row.names=FALSE, col.names=F)
}else{
    write.table(snp, paste0('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/', gene, '/cell_files/', cell, '/SNP.bed'), quote=FALSE, row.names=FALSE, col.names=F)
}
}, error = function(e)
{ matrix(0,1,3)->aa
    cat(aa, file = paste0('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/', gene, '/cell_files/', cell, '/SNP.bed'))
})
