# scRNA-mutation
This project helps us to call variants from single cell exactly, then we use these variants to generate a clonal tree.
# DATA FROM
/N/dc2/projects/ngs/projects/iucmg/10X_singleCell/Nova_S2_008/Analysis/0661-2_CD138Plus/outs/filtered_gene_bc_matrices/GRCh38/
# workflow 
1-preprocess.R


# vatrix
'./vartrix --scoring-method coverage --umi --threads 8 --bam /home/siwxu/bio/clonalTree/data/4/possorted_genome_bam.bam \
--cell-barcodes /home/siwxu/bio/clonalTree/data/4/barcodes.tsv --fasta /home/siwxu/database/10X/genome.fa --out-matrix \
/home/siwxu/bio/clonalTree/data/4/out_matrix_umi --vcf /home/siwxu/bio/clonalTree/data/4/SNP_freebayes_patient4_QUAL20_0728.vcf'
