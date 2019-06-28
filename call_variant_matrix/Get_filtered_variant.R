arg <- commandArgs(T)
arg[1] -> gene
tryCatch({temp1 <- read.table(paste0('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/', gene, '/', gene, '_variant_matrix_rule_4_withname.txt'))
temp2 <- read.table(paste0('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/', gene, '/', gene, '_depth_matrix.txt'))

temp1[,-which(colSums(temp1)==0)]->temp11
temp2[,-which(colSums(temp1)==0)]->temp22

#
temp11+temp22->temp3   #x#

for(i in 1:ncol(temp11)){
  temp11[which(temp3[,i] < 4),i]<-3
}
#
a<-c();b<-c();c<-c();d<-c()
for(i in 1:ncol(temp11)){
  length(which(temp11[,i]==0))->a0
  length(which(temp11[,i]==1))->a1
  length(which(temp11[,i]==2))->a2
  length(which(temp11[,i]==3))->a3
  c(a,a0)->a
  c(b,a1)->b
  c(c,a2)->c
  c(d,a3)->d
}
t(a)->a; t(b)->b; t(c)->c; t(d)->d
rbind(colnames(temp11), a, b, c, d)->temp
t(temp)->variant
#temp[which(as.numeric(temp[,3]) > 1 | as.numeric(temp[,4]) > 1),]->variant

write.table(variant, file = paste0('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/filtered_variant_rule_4.txt'), row.names = F, col.names = F, quote = FALSE, append = TRUE)
} ,error = function(e)
{ 
cat('no snv')
})
