library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(gtools)
library(ggplot2)
options(stringsAsFactors = F)

source("common_functions.R")

## Positions relative to TSS annotations
# get_tss_regions() and bed_to_GRanges() are found in common_functions.R

tss_points <- bed_to_GRanges(get_tss_regions(expand = 1))


peak_files <- paste0(peak_dir,list.files(peak_dir))


cre_order <- data.frame(order = 1:7, cre=c("camk2a","ntsr1","rbp4","scnn1a","cux2","gad2","mes"))

samples <- data.frame(filename = peak_files,
                      cre = c("scATAC",rep("camk2a",2),rep("pvalb",2),rep("vip",2),rep("cux2",3),rep("gad2",3),rep("ntsr1",3),rep("rbp4",3),rep("scnn1a",3),rep("mes",3)))

samples <- samples %>%
  filter(cre %in% cre_order$cre)


bin.df <- data.frame(binmin = seq(0,7,by=7/100)) %>%
  mutate(binmax = lead(binmin,default=7))

results <- bin.df

for(i in 1:nrow(cre_order)) {
  
  cre_samples <- samples %>% filter(cre == cre_order$cre[i])
  
  cre.freq <- matrix(nrow = nrow(bin.df),ncol = nrow(cre_samples))
  
  for(j in 1:nrow(cre_samples)) {
    sample.peaks <- import(cre_samples$filename[j],format="bed")
    sample.dist <- as.data.frame(mcols(distanceToNearest(sample.peaks,tss_points)))
    
    sample.bins <- bin.df %>%
      rowwise() %>%
      mutate(freq = sum(log10(sample.dist$distance) > binmin & log10(sample.dist$distance) <= binmax)/nrow(sample.dist))
    
    sample.freq <- unlist(sample.bins$freq)
    
    cre.freq[,j] <- sample.freq
  }
  
  cre.means <- rowMeans(cre.freq)
  cre.sd <- unlist(apply(cre.freq,1,sd))
  
  results <- cbind(results,cre.means)
  names(results)[ncol(results)] <- paste0(cre_order$cre[i],".m")
  results <- cbind(results,cre.sd)
  names(results)[ncol(results)] <- paste0(cre_order$cre[i],".sd")
  
}
