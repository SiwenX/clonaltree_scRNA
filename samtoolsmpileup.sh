#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=48:00:00
#PBS -l vmem=16gb
#PBS -M siwxu@iu.edu
#PBS -m abe

module load samtools/1.9

n=0
cat /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/100GENES/gene_id.bed | while read line
do
line_arr=( $line )
gen=${line_arr[0]}

cat /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/possorted_genome_bam_header.bed $gen  \
>> /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/100GENES/$gen.sam
samtools view -b $gen.sam > $gen.bam
samtools sort $gen.bam > $gen.sort.bam
samtools index $gen.sort.bam 
samtools mpileup -l /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/100GENES/gene_region.bed --ff UNMAP,SECONDARY,QCFAIL -uvf /N/dc2/projects/ngs/users/siwxu/database/10X/genome.fa -t DP,AD \
/N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/100GENES/$gen.sort.bam > /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/100GENES/$gen.vcf
((n=n+1))
echo $n
done

