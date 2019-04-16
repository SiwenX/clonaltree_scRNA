#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=120:00:00
#PBS -l vmem=32gb
#PBS -M siwxu@iu.edu
#PBS -m abe

module load samtools/1.9
module load vcftools/gnu/0.1.13

n=0
cat /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/gene_region.bed | while read line
do
line_arr=( $line )
gene=${line_arr[0]}
mkdir /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/$gene
mkdir /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/$gene/vcf
    cat /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cells2/barcodes/barcode.tsv | while read line
    do
    line_arr2=( $line )
    id=${line_arr2[0]}
    grep "$id" /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/100GENES/$gene.sam > /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/$gene/$id
    cat /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/possorted_genome_bam_header.bed /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/$gene/$id > /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/$gene/$id.sam
    samtools view -b /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/$gene/$id.sam > /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/$gene/$id.bam
    samtools sort /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/$gene/$id.bam > /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/$gene/$id.sort.bam
    samtools index /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/$gene/$id.sort.bam
    samtools mpileup --ff UNMAP,SECONDARY,QCFAIL -uvf /N/dc2/projects/ngs/users/siwxu/database/10X/genome.fa -t DP,AD \
    /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/$gene/$id.sort.bam > /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/$gene/vcf/$id.vcf
    vcftools --vcf /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/$gene/vcf/$id.vcf --min-alleles 3 --recode --out /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/$gene/vcf/$id
    done
((n=n+1))
echo $n
done


