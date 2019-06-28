arg <- commandArgs(T)
arg[1] -> gene
arg[2] -> cell
setwd(paste0('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/', gene, '/cell_files/', cell))
tryCatch({temp<-read.table("temp.txt",header=F)} ,error = function(e)
{
  cat('no input')})
tryCatch({filter<-read.table("filter_variant.txt",header=F)
    temp[as.numeric(filter[,1]),]->temp
    write.table(temp,paste0('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/', gene, '/cell_files/', cell, '/temp2.txt'), quote=FALSE, row.names=F, col.names=F)
} ,error = function(e)
{
    cat('no snv')
    temp <- c()
    write.table(temp,paste0('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/', gene, '/cell_files/', cell, '/temp2.txt'), quote=FALSE, row.names=F, col.names=F)

})
