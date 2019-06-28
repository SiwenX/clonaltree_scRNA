arg <- commandArgs(T)
arg[1] -> cell_name
arg[2] -> gene_name
#arg[3] -> gene_pos
#arg[4] -> gene_sta
#arg[5] -> gene_end
#c(1,as.numeric(gene_pos),as.numeric(gene_sta),as.numeric(gene_end))->gene
library("GenomicRanges")
library(stringr)
gene<-read.table(paste0('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/gene_list/',gene_name),header=FALSE)
gene[1,1]->GeneName
#gene_name -> GeneName
as.numeric(gene[4])-as.numeric(gene[3])+1 -> length
rep(0,length) -> aaa
tryCatch({cell<-read.table(paste0('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/', GeneName, '/cell_files/', cell_name, '/temp2.txt'),header=F)
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
    cat(cell_name, file = paste('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/', GeneName, '/cell_files/',cell_name, '/', GeneName,'_variant_matrix_rule_4.txt',sep = ''), append = TRUE)
    cat("\t", file = paste('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/', GeneName, '/cell_files/',cell_name, '/', GeneName,'_variant_matrix_rule_4.txt',sep = ''), append = TRUE)
    cat(gene_leng, file = paste('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/', GeneName, '/cell_files/',cell_name, '/', GeneName,'_variant_matrix_rule_4.txt',sep = ''), append = TRUE)
    cat("\n", file = paste('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/', GeneName, '/cell_files/',cell_name, '/', GeneName,'_variant_matrix_rule_4.txt',sep = ''), append = TRUE)
} ,error = function(e)
{
    cat(cell_name, file = paste('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/', GeneName, '/cell_files/',cell_name, '/', GeneName,'_variant_matrix_rule_4.txt',sep = ''), append = TRUE)
    cat("\t", file = paste('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/', GeneName, '/cell_files/',cell_name, '/', GeneName,'_variant_matrix_rule_4.txt',sep = ''), append = TRUE)
    cat(aaa, file = paste('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/', GeneName, '/cell_files/',cell_name, '/', GeneName,'_variant_matrix_rule_4.txt',sep = ''), append = TRUE)
    cat("\n", file = paste('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/', GeneName, '/cell_files/',cell_name, '/', GeneName,'_variant_matrix_rule_4.txt',sep = ''), append = TRUE)
})
