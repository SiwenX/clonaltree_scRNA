#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=48:00:00
#PBS -l vmem=16gb
#PBS -M siwxu@iu.edu
#PBS -m abe

module load samtools/1.9
module load r/3.3.1

touch /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/RPL10/vcf/genodepth_matrix_0.bed
n=0
m=0
cat /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/barcodes.tsv | while read line
do
line_arr=( $line )
cell=${line_arr[0]}
samtools depth -a /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/RPL10/$cell.sort.bam -r X:154389955-154409168 > /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/RPL10/vcf/temp.txt
cat /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/RPL10/vcf/temp.txt | awk '{print $3}' > /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/RPL10/vcf/temp2.txt
((m=m+1))
paste /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/RPL10/vcf/genodepth_matrix_$n.bed /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/RPL10/vcf/temp2.txt > /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/RPL10/vcf/genodepth_matrix_$m.bed
((n=n+1))
done
