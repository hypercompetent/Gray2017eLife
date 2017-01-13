library(DiffBind)
library(dplyr)

# Build sample table from BAM files
sample_files <- list.files("bam")
control_files <- sample_files[grepl("genomic",sample_files)]
sample_files <- sample_files[!grepl("genomic",sample_files)]

samples <- data.frame(bamReads = paste0("bam/",sample_files)) %>%
  mutate(SampleID = sub("bam/","",sub("_3.2M.bam","",bamReads))) %>%
  mutate(Tissue = sub("_.+","",SampleID)) %>%
  mutate(Tissue = ifelse(Tissue == "wt","mes",Tissue)) %>%
  mutate(Replicate = sub(".+rep","",SampleID)) %>%
  mutate(bamControl = paste0("bam/",control_files[1])) %>%
  mutate(Peaks = sub(".bam",".narrowPeak",sub("bam/","hotspot/",bamReads))) %>%
  mutate(peakCaller = "narrowPeak")

# Build the Diffbind Object (dbo)
dbo <- dba(sampleSheet=samples)

#  Run the count function to get results based on affinity (TMM), 
#  not just peak overlaps
dbo <- dba.count(dbo, minOverlap = 3)

# Write the affinity-based histogram to a PDF for use in Figure 3A
pdf("figure_elements/Figure_3A.pdf", height = 8, width = 8)
dba.plotHeatmap(dbo, colScheme = "Blues")
dev.off()

## DESeq2 Contrast analysis for each pair of tissues
dbo <- dba.contrast(dbo, categories = DBA_TISSUE, minMembers = 3)
dbo <- dba.analyze(dbo, method = DBA_DESEQ2)

save(dbo,samples, file = "diffbind_results.rda")
