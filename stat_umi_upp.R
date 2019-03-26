arg <- commandArgs(T)
library(stringr)
setwd('/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/')
sam<-read.table("fragments.bed",header=FALSE)
SNP<-read.table("SNP.bed",header=F)
temp<-read.table("temp.txt",header=F)
tryCatch({temp<-read.table("temp.txt",header=F)} ,error = function(e)
{ aa<-c()
  cat(aa, file = '/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/filter_variant.txt', append = TRUE)
})

SNP[arg[1],]->snp;
filter_loci<-c()

str_split(snp$V4,',')[[1]][1]->snp[4]
colnames(snp) <- c('chr','pos','ref','alt')
unique(sam$V6)->UMI
adjust_cigar<-function(start,pos,cigar){
    len<-as.numeric(pos)-as.numeric(start)+1;
    
    num<-as.numeric(str_split(cigar[[1]],"D|I|M|N|S")[[1]]);
    num<-num[!is.na(num)]
    label<-(str_split(cigar[[1]],"\\d")[[1]] );
    label<-label[nchar(label)>0]
    
    cur_len<-0;insert_len<-0;del_len<-0;N_len<-0;S_len<-0
    for(i in 1:length(num)){
        
        if(label[i]=="I"){
            insert_len<-insert_len+num[i];
        }
        
        if(label[i]=="D"){
            del_len<-del_len+num[i];
        }
        
        if(label[i]=="M" || label[i]=="D"){
            cur_len<-cur_len+num[i]
        }
        
        if(label[i]=="N"){
            N_len<-N_len+num[i];
            cur_len<-cur_len+num[i]
        }
        
        if(label[i]=="S"){
            cur_len<-cur_len;
            S_len<-S_len+num[i];
        }
        
        
        if(cur_len>=len){
            break;
        }
        
        
    }
    
    return(as.numeric(pos)+as.numeric(insert_len)-as.numeric(del_len)-as.numeric(N_len)+as.numeric(S_len))
    
}

myfunction2 <- function(x){
    loci<-adjust_cigar(x[3],x[7],x[4])
    strsplit(as.character(x[5][[1]]), split = "", fixed = TRUE)[[1]]->seq
    seq[as.numeric(loci)-as.numeric(x[3])+1]->SNP
    if(length(SNP) > 1){
        "NA" -> SNP
    }
    return(SNP)
}

a <- c()
for(i in 1:length(UMI)){
    sam[which(sam$V6 == as.character(UMI[i])),]->subsam
    snp$pos->subsam[,7]
    apply(subsam, 1, myfunction2)->SNP
    as.numeric(length(which(SNP == snp$ref)))->ref
    length(which(SNP == snp$alt))->alt
    c((UMI[i]), alt/(ref+alt), ref, alt) -> aa     #为了得到numeric形式，UMI现在是数；如果想得到字符形式，改为as.character(UMI[i])
    rbind(a, aa) -> a
}

sum(a[,3])->ref_number
sum(a[,4])->alt_number

#vote function
##alt_ratio > 0.5 & number of alt_UMI > 2
if(ref_number<alt_number){
  which(a[,3]!=0)->ref_UMI
  if(length(ref_UMI) > 1 & all(as.numeric(a[ref_UMI,2]) > 0.4) ){arg[1]->keep_loci}
}
if(ref_number>alt_number){
  which(a[,4]!=0)->alt_UMI
  if(length(alt_UMI) > 1 & all(as.numeric(a[alt_UMI,2]) > 0.4) ){arg[1]->keep_loci}
}


if(!is.null(keep_loci)){
  cat(keep_loci, file = '/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/filter_variant.txt', append = TRUE)
  cat("\n", file = '/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/filter_variant.txt', append = TRUE)
}
