#!/bin/bash

n=0
cat /Users/xsw/bio/scRNA/data/2/6166cell/vcf/Cell_id2.txt | while read line
do
line_arr=( $line )
cell=${line_arr[0]}
mkdir /Users/xsw/bio/scRNA/data/2/6166cell/vcf/$cell
grep -v '#' $cell.recode.vcf > /Users/xsw/bio/scRNA/data/2/6166cell/vcf/$cell/temp.txt

((n=n+1))
echo $n
echo $cell
Rscript /Users/xsw/bio/scRNA/data/2/6166cell/variant_matrix.R $cell

done


