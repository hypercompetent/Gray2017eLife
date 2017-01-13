#!/bin/bash

## Raw fastq files from SRA
## Only necessary if you want to redo alignment and peak calling.

mkdir fastq

# sratools fastq-dump location
FQD=/data/mct-t200/T502/atac_analysis/sratools/sratoolkit.2.8.1-centos_linux64/bin/fastq-dump

# SRA accessions for Gray, et al. (2017)
# SRA Study SRP090744
# GEO Accession GSE87548
ACCLIST=(SRR4340343 SRR4340344 SRR4340345 SRR4340346 SRR4340347 SRR4340348 SRR4340349 SRR4340359 SRR4340362 SRR4340366 SRR4340364 SRR4340363 SRR4340361 SRR4340360 SRR4340358 SRR4340357 SRR4340356 SRR4340355 SRR4340354 SRR4340353 SRR4340352 SRR4340351 SRR4340350 SRR4340370 SRR4340369 SRR4340368 SRR4340367 SRR4340365)

# File names for FASTQ files
NAMELIST=(cux2_rep1 cux2_rep2 cux2_rep3_L1 cux2_rep3_L2 gad2_rep1 gad2_rep2 gad2_rep3 rbp4_rep3_L2 scnn1a_rep3_L2 mes_rep1 mes_rep2 scnn1a_rep3_L1 scnn1a_rep2 scnn1a_rep1 rbp4_rep3_L1 rbp4_rep1_L1 rbp4_rep1_L2 rbp4_rep2_L1 rbp4_rep2_L2 ntsr1_rep3_L1 ntsr1_rep3_L2 ntsr1_rep1 ntsr1_rep2 genomic_rep2_L1 genomic_rep2_L2 genomic_rep1_L1 genomic_rep1_L2 mes_rep3)

# Retrieve FASTQ files directly using fastq-dump
# See https://github.com/ncbi/sra-tools/wiki/Download-On-Demand
N=${#ACCLIST[@]}
for ((i = 0; i < $N; i++))
do
 	ACC=${ACCLIST[$i]}
 	BASE=${NAMELIST[$i]}
	$FQD --split-files --outdir fastq --gzip $ACC
	mv fastq/${ACC}_1.fastq.gz fastq/${BASE}_R1.fastq.gz
	mv fastq/${ACC}_2.fastq.gz fastq/${BASE}_R2.fastq.gz
done

# Merge multiple lanes for a few samples to simplify alignment
MERGES=(cux2_rep3 ntsr1_rep3 rbp4_rep1 rbp4_rep2 rbp4_rep3 scnn1a_rep3 genomic_rep1 genomic_rep2)
for BASE in ${MERGES[*]}
do
	zcat fastq/${BASE}_L1_R1.fastq.gz fastq/${BASE}_L2_R1.fastq.gz | gzip -c > fastq/${BASE}_R1.fastq.gz
	zcat fastq/${BASE}_L1_R2.fastq.gz fastq/${BASE}_L2_R2.fastq.gz | gzip -c > fastq/${BASE}_R2.fastq.gz
	rm fastq/${BASE}_L1_R1.fastq.gz
	rm fastq/${BASE}_L2_R1.fastq.gz
	rm fastq/${BASE}_L1_R2.fastq.gz
	rm fastq/${BASE}_L2_R2.fastq.gz
done


## Mo 2015 files
## Used for comparisons to previous ATAC-Seq study.
MD=mo_2015_files
mkdir $MD

# Histone modification ChIP-Seq from Camk2a-Cre cells for comparisons to ATAC-seq data
wget -O ${MD}/H3K27ac_peaks.bed.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63137&format=file&file=GSE63137%5FChIP%2Dseq%5FH3K27ac%5Fexcitatory%5Fneurons%5FSICER%5Fpeaks%2Ebed%2Egz"
wget -O ${MD}/H3K27me3_peaks.bed.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63137&format=file&file=GSE63137%5FChIP%2Dseq%5FH3K27me3%5Fexcitatory%5Fneurons%5FSICER%5Fpeaks%2Ebed%2Egz"
wget -O ${MD}/H3K4me1_peaks.bed.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63137&format=file&file=GSE63137%5FChIP%2Dseq%5FH3K4me1%5Fexcitatory%5Fneurons%5FSICER%5Fpeaks%2Ebed%2Egz"
wget -O ${MD}/H3K4me3_peaks.bed.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63137&format=file&file=GSE63137%5FChIP%2Dseq%5FH3K4me3%5Fexcitatory%5Fneurons%5FSICER%5Fpeaks%2Ebed%2Egz"

# JASPAR Motif Database for AME and FIMO
mkdir meme_db
wget -o meme_db/motif_databases.12.15.tgz "http://meme-suite.org/meme-software/Databases/motifs/motif_databases.12.15.tgz"
tar -zxvf meme_db/motif_databases.12.15.tgz motif_databases/JASPAR/JASPAR_CORE_2016.meme
cat motif_databases/JASPAR/JASPAR_CORE_2016.meme neurod2_motif.meme > temp.meme
mv temp.meme motif_databases/JASPAR/JASPAR_CORE_2016.meme
