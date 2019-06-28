arg <- commandArgs(T)
arg[1] -> cell_name
arg[2] -> gene_name
arg[3] -> gene_id

gene<-read.table(paste0('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/gene_list/', gene_id),header=FALSE)
gene[1,1]->GeneName
as.numeric(gene[4])-as.numeric(gene[3])+1 -> length
rep(0,length) -> aaa
t(aaa)->aaa
tryCatch({temp2<-read.table(paste0('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/', gene_name, '/cell_files/', cell_name, '/temp2.txt'),header=FALSE)
    t(temp2)->temp2
write.table(temp2, file = paste0('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/', gene_name, '/cell_files/', cell_name, '/genodepth_matrix.txt'), col.names = F, row.names = F)
    
    
    
} ,error = function(e)
{   cat('e')
    write.table(aaa, file = paste0('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/', gene_name, '/cell_files/', cell_name, '/genodepth_matrix.txt'), col.names = F, row.names = F)
})
