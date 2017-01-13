#!/bin/bash

## Downsample BAM Files
# millions of fragments (read pairs) to retain
DS=3.2

# BAM directory
BD=bam
# Stats directory
SD=stats
#Downsample directory
DD=bam_${DS}M

# Make downsample directory
mkdir $DD

# Run downsampling scripts
for F in ${BD}/*.rmd.srt.bam
do
	NODIR=$(basename ${F})
	BASE=`echo "$NODIR" | cut -d'.' -f1`
	
	# Sort the bame file by name to group pairs together
	#samtools sort -n ${BD}/${BASE}.rmd.srt.bam ${BD}/${BASE}_byname
	# Randomly select pairs to retain
	Rscript select_downsample_lines.R ${SD}/${BASE}.stats $DS keep_lines.txt
	# Pull header from the bam file to start sam output
	samtools view -H ${BD}/${BASE}_byname.bam > temp.sam
	# Filter using perl
	samtools view ${BD}/${BASE}_byname.bam | ./sam_subsample_filter_stdio.pl keep_lines.txt >> temp.sam
	# Convert back to bam for downstream use
	samtools view -bS temp.sam > temp.bam
	# Sort filtered BAM file
	samtools sort temp.bam ${DD}/${BASE}
	# Remove intermediate files
	#rm keep_lines.txt
	#rm temp.sam
	#rm temp.bam
done

