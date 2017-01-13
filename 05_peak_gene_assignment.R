library(GenomicRanges)
library(dplyr)
library(readxl)
library(DiffBind)
library(ggplot2)
source("common_functions.R")
options(stringsAsFactors = F)

# diffbind_results were computed in 01_DiffBind_Analysis.R
load("diffbind_results.rda")

cell.classes <- c("cux2","gad2","ntsr1","rbp4","scnn1a","mes")

# Fetch TSS locations from UCSC refGene table
# get_tss_regions() and bed_to_GRanges() are in common_functions.R
tss <- get_tss_regions(expand = 0) %>%
  unique() %>%
  filter(!grepl("_",chr)) %>%
  bed_to_GRanges()

save(tss, file = "tss.rda")

# Pull samples, peaks, and p-values from the DiffBind object (dbo)
samples <- dbo$samples
peaks <- dbo$peaks
peak.sig <- dbo$binding[,-c(1:3)]

# Convert the peak locations to a GRanges object
peak_chr <- dbo$chrmap[dbo$binding[,1]]
peak_data <- as.data.frame(dbo$binding)
all.peaks <- GRanges(IRanges(start = peak_data$START, 
                             end = peak_data$END), 
                     seqnames = peak_chr)

# Calculate mean TMM values for each cell class
peak.cre <- do.call("cbind", 
                    tapply(samples$SampleID, 
                           samples$Tissue, 
                           function(x) rowMeans(peak.sig[,x])))

peak.cre <- peak.cre[,-6]

# Load gene expression data from supplemental table 7
counts <- read_excel("Supp Table 7 Gene Expression.xlsx")
gene_names <- counts$gene
counts <- as.matrix(counts[,-1])
rownames(counts) <- gene_names
expr <- log2(counts + 1)

# Gene expression sample info
expr.sample.df <- data.frame(sample_id = colnames(expr),
                             cre = sub("_.+","",colnames(expr)))

# Get means from scRNA-Seq for each cell class.
expr.cre = do.call("cbind",
                   tapply(expr.sample.df$sample_id, 
                          expr.sample.df$cre,
                          function(x) rowMeans(counts[,x])))

# Filter expression and TSS locations for genes that are found both in the expression data
# and in the TSS table
common.genes <- intersect(row.names(expr.cre),tss$mcols.name)
expr.cre <- expr.cre[common.genes,]
expr.cre <- log2(expr.cre + 1)
tss <- tss[tss$mcols.name %in% common.genes]
counts <- counts[common.genes,]

# Get distance from peaks to the nearest TSS
# getNearestDist() is in common_functions.R
tmp <- getNearestDist(all.peaks,tss,return.nearest=T)
all.peaks$gene.id  <- tmp$mcols.name
all.peaks$dist2tss <- tmp$dist

# Compute correlations between peak accessibility and gene expression

# numerical vector correlation function, vecCor
vecCor <- function(A,B)
{
  cA <- A - rowMeans(A)
  cB <- B - rowMeans(B)
  sA <- sqrt(rowMeans(cA^2))
  sB <- sqrt(rowMeans(cB^2))
  rowMeans(cA * cB) / (sA * sB)
}

# correlation between peak accessibility and gene expression across all cell classes
gene.cor <- vecCor(peak.cre, expr.cre[all.peaks$gene.id,])
gene.cor[is.na(gene.cor)] <- 0
gene.cor.df <- data.frame(cor = gene.cor)

# Permutation of correlations
perm.cor <- sapply(1:10, function(i){
  # Randomly reassign gene names to randomize associations
  p <- setNames(sample(common.genes), common.genes)
  # Calculate permuted correlations
  gene.cor <- vecCor(peak.cre, expr.cre[p[all.peaks$gene.id],])
})

# Plotting actual vs permuted correlations
perm.cor[is.na(perm.cor)] <- 0
colnames(perm.cor) <- paste0("perm", 1:10)
perm.cor.df <- as.data.frame(as.table(perm.cor))
colnames(perm.cor.df) <- c("site", "cat", "cor")

cor.df <- rbind(perm.cor.df, 
                data.frame(site = names(gene.cor), 
                           cat = rep("obs", length(gene.cor)), 
                           cor = gene.cor))
# Density Plot
perm.cor.density <- ggplot(cor.df,
                        aes(x = cor,
                            color = cat, 
                            fill = cat,
                            alpha = 0.5)) + 
  geom_density() + 
  scale_fill_manual(values = c(rep("blue",10), "red"))

# Box Plot
perm.core.boxplot <- ggplot(cor.df,
                            aes(x = rev(cat),
                                y = cor)) +
  geom_boxplot() +
  scale_fill_manual(values = c(rep("gray80",10), "skyblue")) +
  theme_classic()

# Write results for the next step.
save(all.peaks, peak.sig, peak.cre, gene.cor, file = "atac_data.rda")
save(counts, expr.sample.df, expr, expr.cre, file = "expr_data.rda")



