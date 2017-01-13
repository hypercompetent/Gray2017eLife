#!/bin/bash

## Directories
# Fastq Directory
FD=fastq
# BAM Directory
BD=bam
# Stats Directory
SD=stats
# Preseq Directory
PD=preseq
# Insert Sizes Directory
ID=inserts

## Applications
# Bowtie Location
BT=bowtie
# Trim Galore location
TG=/data/mct-t200/T502/Lab_Notebook/tools/trim_galore
# Samtools location
ST=samtools
# Preseq location
PS=preseq
# bamToBed location (from bedtools)
BB=bamToBed
# Collect Insert Size Metrics location (from Picard Tools)
CI="java -jar /usr/bin/bioapps/picard/picard-tools-1.110/CollectInsertSizeMetrics.jar"

#Bowtie settings
BS="-p 6 -m 1 -S -X 2000 --chunkmbs 256 --un unaligned.fastq"
# Genome Index Location
GL=/data/rnaseqanalysis/RNAseq/indexes/mm10/genome

## Make directories
mkdir $BD
mkdir $SD
mkdir $PD
mkdir $ID

## Perform alignments
R1=($FD/*_R1*)
R2=($FD/*_R2*)

PAIRS=${#R1[@]}

echo $PAIRS

for ((i = 0; i < $PAIRS; i++))
do
	FASTQ=$(basename ${R1[i]})
	BASE=`echo "$FASTQ" | cut -d'.' -f1`
	OUTBASE=${BASE%_L1*}
	echo $OUTBASE
	# Initial bowtie alignment
	$BT $BS $GL -1 <( zcat ${R1[i]} ) -2 <( zcat ${R2[i]} ) aligned.sam
	# Convert SAM results to BAM
	$ST view -Sb aligned.sam > aligned.bam
	# Trim Nextera primer sequences from unaligned reads
	$TG --nextera --paired --gzip --three_prime_clip_R1 1 --three_prime_clip_R2 1 unaligned_1.fastq unaligned_2.fastq
	# Align trimmed reads
	$BT $BS $GL -1 <( zcat unaligned_1_val_1.fq.gz ) -2 <( zcat unaligned_2_val_2.fq.gz ) aligned.sam
	# Convert trimmed read results to BAM
	$ST view -Sb aligned.sam > trimmed.bam
	# Concatenate initial and trimmed results
	$ST cat -o merged.bam aligned.bam trimmed.bam
	# Retain only aligned pairs based on BAM flags
	$ST view -b -F 4 merged.bam > mapped.bam
	# Sort mapped reads
	$ST sort mapped.bam ${BD}/${OUTBASE}.srt
	# Remove duplicates
	$ST rmdup ${BD}/${OUTBASE}.srt.bam ${BD}/${OUTBASE}.rmd.srt.bam
	# Index the unique, sorted reads
	$ST index ${BD}/${OUTBASE}.rmd.srt.bam
	# Compute stats for all mapped reads and for reads after removing duplicates
	$ST flagstat merged.bam > ${SD}/${OUTBASE}.stats
	$ST flagstat ${BD}/${OUTBASE}.rmd.srt.bam >> ${SD}/${OUTBASE}.stats
	# Remove intermediate files
	rm aligned.sam
	rm aligned.bam
	rm mapped.bam
	rm trimmed.bam
	rm merged.bam
	rm ./unaligned*
	# Perform preseq on sorted reads to assess library diversity
	$PS lc_extrap -P -o ${PD}/${OUTBASE}.txt <( $BB -i ${BD}/${OUTBASE}.srt.bam )
	# Perform insert size analysis on unique, sorted reads
	$CI I=${BD}/${OUTBASE}.rmd.srt.bam O=${ID}/${OUTBASE}_insert_sizes.txt H=${ID}/${OUTBASE}.pdf
done
