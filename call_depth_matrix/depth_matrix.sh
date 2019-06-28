#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -M siwxu@iu.edu
#PBS -m abe
#PBS -l walltime=12:00:00

gene_id=$gene_id

module load samtools/1.3
module load r/3.3.1

cat /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/gene_list/$gene_id | while read line
do
line_arr3=( $line )
gene=${line_arr3[0]}
gene_pos=${line_arr3[1]}
gene_sta=${line_arr3[2]}
gene_end=${line_arr3[3]}
mkdir /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/cell_files
j=0
pids=""
for i in $(cat /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/barcodes.tsv); do
(
mkdir /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/cell_files/$i
#rm /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/Cell3/$gene/cell_files/$i/$gene'_variant_matrix_rule_t.txt'
samtools depth -a /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/$i.sort.bam -r $gene_pos:$gene_sta-$gene_end > /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/cell_files/$i/temp.txt
cat /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/cell_files/$i/temp.txt | awk '{print $3}' > /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/cell_files/$i/temp2.txt
Rscript  /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/src/genodepth_matrix.R $i $gene $gene_id

) &
pids=`echo $! $pids`
j=$(($j+1))
if (( $j % 6 == 0 )); then wait $pids; pids=""; fi
done
wait $pids
for i in $(cat /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/barcodes.tsv); do
cat /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/cell_files/$i/genodepth_matrix.txt >> /N/dc2/projects/ngs/users/siwxu/projects/scRNA/data/4/CD138plus/HighGene/$gene/$gene'_depth_matrix.txt'
done
done

