#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -M siwxu@iu.edu
#PBS -m abe
#PBS -l walltime=12:00:00

gene_id=$gene_id

module load samtools/1.3
module load vcftools/0.1.13
module load r/3.3.1

cat /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/gene_list/$gene_id | while read line
do
    line_arr3=( $line )
    gene=${line_arr3[0]}
    gene_pos=${line_arr3[1]}
    gene_sta=${line_arr3[2]}
    gene_end=${line_arr3[3]}
    mkdir /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/cell_files
    rm /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/$gene'_variant_matrix_rule_4.txt'
    j=0
    pids=""
    for i in $(cat /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/barcodes.tsv); do
        (
        mkdir /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/cell_files/$i
        rm /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/cell_files/$i/$gene'_variant_matrix_rule_4.txt'
        rm /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/cell_files/$i/temp.txt
        rm /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/cell_files/$i/SNP.bed
        rm /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/cell_files/$i/fragments.bed
        rm /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/cell_files/$i/filter_variant.txt
        grep -v '#' /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/vcf/$i.recode.vcf > /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/cell_files/$i/temp.txt
        cat /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/cell_files/$i/temp.txt | awk '{print $1 "\t" $2 "\t" $4 "\t" $5}' > /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/cell_files/$i/SNP.bed
        Rscript  /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/src/filter_snp.R $gene $i
        m=0
        cat /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/cell_files/$i/SNP.bed | while read line
        do
            line_arr2=( $line )
            pos=${line_arr2[0]}
            sta=${line_arr2[1]}
            samtools view /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/$i.sort.bam $pos:$sta-$sta | awk '{print $2 "\t" $3 "\t" $4 "\t" $6 "\t" $10 "\t" $(NF-1)}' > /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/cell_files/$i/fragments.bed
            ((m=m+1))
            Rscript /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/src/stat_umi_upp.R $m $gene $i
        done
        Rscript /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/src/filter_variant.R $gene $i
        Rscript /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/src/variant_matrix.R $i $gene_id
        ) &
        pids=`echo $! $pids`
        j=$(($j+1))
       if (( $j % 16 == 0 )); then wait $pids; pids=""; fi
    done
    wait $pids
    for i in $(cat /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/barcodes.tsv); do
        cat /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/cell_files/$i/$gene'_variant_matrix_rule_4.txt' >> /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/$gene'_variant_matrix_rule_4.txt'
    done
    Rscript /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/src/add_rowcol_name.R $gene $gene_id
    Rscript /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/src/Get_filtered_variant.R $gene
done
