#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=48:00:00
#PBS -l vmem=16gb
#PBS -M siwxu@iu.edu
#PBS -m abe

module load samtools/1.9
module load vcftools/gnu/0.1.13
module load r/3.3.1

n=0
cat /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/barcodes.tsv | while read line
do
line_arr=( $line )
cell=${line_arr[0]}
grep -v '#' /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/$cell.recode.vcf > /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/temp.txt
cat /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/temp.txt | awk '{print $1 "\t" $2 "\t" $4 "\t" $5}' > /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/SNP.bed
Rscript  /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/filter_snp.R
m=0
cat /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/SNP.bed | while read line
do
line_arr2=( $line )
pos=${line_arr2[0]}
sta=${line_arr2[1]}
echo $pos
echo $sta
samtools view /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/$cell.sort.bam $pos:$sta-$sta | awk '{print $2 "\t" $3 "\t" $4 "\t" $6 "\t" $10 "\t" $(NF-1)}' > /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/fragments.bed
((m=m+1))
echo $m
Rscript /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/stat_umi_upp.R $m
done
((n=n+1))
echo $n
Rscript /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/filter_variant.R
Rscript /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/variant_matrix.R $cell

rm -f /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/filter_variant.txt
rm -f /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/temp.txt
rm -f /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/SNP.bed
rm -f /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/fragments.bed

done
Rscript /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/IGLC2/vcf/add_rowcol_name.R


