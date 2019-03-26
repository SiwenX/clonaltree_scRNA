setwd('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/')
tryCatch({temp<-read.table("temp.txt",header=F)} ,error = function(e)
{
  cat('no input')})
tryCatch({filter<-read.table("filter_variant.txt",header=F)} ,error = function(e)
{
    cat('no snv')})
temp[as.numeric(filter[,1]),]->temp
write.table(temp,'/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/temp.txt',quote=FALSE,row.names=F,col.names=F)
