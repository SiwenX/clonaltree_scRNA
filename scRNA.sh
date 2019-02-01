#!/bin/bash

n=0
cat SNP.bed | while read line
do
line_arr=( $line )
pos=${line_arr[0]}
sta=${line_arr[1]}

samtools view /Users/xsw/bio/scRNA/data/2/mutation_call/AAACCTGAGAAGGTGA-1.sort.bam $pos:$sta-$sta | awk '{print $2 "\t" $3 "\t" $4 "\t" $6 "\t" $10 "\t" $(NF-1)}' > /Users/xsw/bio/scRNA/data/2/mutation_call/fragments.bed
grep -v "S" /Users/xsw/bio/scRNA/data/2/mutation_call/fragments.bed > /Users/xsw/bio/scRNA/data/2/mutation_call/fragments_noS.bed
((n=n+1))
echo $n
Rscript scRNA.R $pos $sta $n

done
