arg <- commandArgs(T)
arg[1] -> gene_name
arg[2] -> gene_id
Gene<-read.table(paste0("/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/gene_list/",gene_id),header=FALSE)
RM<-read.table(paste0('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/', gene_name, '/',gene_name,'_variant_matrix_rule_4.txt'))
rownames(RM) <- RM[,1]
RM[,-1] -> RM
paste(Gene[[1]],c((as.numeric(Gene[3])):as.numeric(Gene[4])),sep = '') -> colnames(RM)
write.table(RM,paste0('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/', gene_name, '/',gene_name,'_variant_matrix_rule_4_withname.txt'),quote=FALSE,row.names=T,col.names=T)
