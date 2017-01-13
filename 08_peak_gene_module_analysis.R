library(GenomicRanges)
library(dplyr)
library(readxl)
library(DiffBind)
library(ggplot2)
library(matrixStats)
library(reshape2)
source("common_functions.R")
options(stringsAsFactors = F)

load("atac_data.rda")
load("expr_data.rda")
load("all_peaks.rda")
load("de_gene_expr.rda")

all.peaks$peak_id <- 1:length(all.peaks)

## Initial Peak accessibility k-means clustering

# select columns for glutamatergic peaks (start with p, exclude Gad2 and mES)
p.cols <- grep("^p", colnames(values(all.peaks)), value = T)
p.cols <- setdiff(p.cols, c(grep("gad2", p.cols, value=T), c(grep("wt", p.cols, value = T))))
p.cols <- grep("_",p.cols, value=T)

# filter peaks for minimum p-value < 0.01, and > 4-fold differential accessibility ( > 2 on log2 scale)
ex.peak.cre <- peak.cre[,-2]
select <- all.peaks$ex.pPadj < 0.01 & rowMaxs(ex.peak.cre) - rowMins(ex.peak.cre) > 2

# normalize peak accessibility values by dividing by subtracting minimum accessibility score and 
# dividing by maximum mean accessibility score per peak
ex.peak.norm <- ex.peak.cre[select,] - rowMins(ex.peak.cre[select,])
ex.peak.norm <- ex.peak.norm/rowMaxs(ex.peak.norm)

## Initial k-means clustering.
# will give different results each run.
tmp <- kmeans(ex.peak.norm, 12, iter.max = 100)
cl <- tmp$cluster
ex.peak.cl <- cl
ex.peak.cl.means <- round(do.call("rbind", 
                                  tapply(1:length(cl), 
                                         cl, 
                                         function(x){
                                           colMeans(ex.peak.norm[x,])
                                          })
                                  ),
                          digits = 2)


# heatmap of peaks in k-means
ex.peak.cor <- cor(t(ex.peak.norm), t(ex.peak.cl.means))
ex.peak.cl <- apply(ex.peak.cor, 1, which.max)

# too many peaks to plot, so choose the top 500 from each cluster
ex.peak.plot <- data.frame(ex.peak.norm, cl = ex.peak.cl, max = rowMaxs(ex.peak.cor))
ex.peak.plot <- ex.peak.plot %>%
  group_by(cl) %>%
  arrange(-max) %>%
  filter(row_number() <= 500) %>%
  arrange(cl, -max)

ex.peak.plot2 <- melt(ex.peak.plot %>% select(-max), "cl")

ggplot(ex.peak.plot2) +
  geom_tile(aes(x = rep(1:6000,4), y = variable, fill = value))


## Initial Gene expression k-means clustering

# Select columns for glutamatergic gene expression (exclude Gad2)
tmp <- names(de.expr)
tmp <- setdiff(tmp, grep("gad",tmp, value=T))

# filter genes for minimum p-value < 0.05 and 4-fold differential expression ( > 1 on log2 scale)
ex.de.g <- unique(unlist(sapply(de.expr[tmp], 
                                function(df) {
                                  row.names(df)[with(df, which(abs(lfc) > 1 & padj < 0.05))]
                                },
                                simplify = F)
                         )
                  )

ex.de.g <- intersect(ex.de.g, row.names(expr.cre))
tmp.dat <- expr.cre[ex.de.g, -2]

# normalize gene expression values by dividing by subtracting minimum mean expression values and 
# dividing by maximum mean expression values per peak
tmp.dat <- tmp.dat - rowMins(tmp.dat)
ex.expr.norm <- tmp.dat / rowMaxs(tmp.dat)

## Initial k-means clustering.
# will give different results each run.
tmp <- kmeans(ex.expr.norm, 12, iter.max = 100)
cl <- tmp$cluster
ex.expr.cl.means <- round(do.call("rbind",
                                  tapply(1:length(cl), 
                                         cl, 
                                         function(x){
                                           colMeans(ex.expr.norm[x,])
                                          })
                                  ),
                          digits = 2)

# heatmap of genes in k-means
ex.gene.cor <- cor(t(ex.expr.norm), t(ex.expr.cl.means))
ex.gene.cl <- apply(ex.gene.cor, 1, which.max)

ex.gene.plot <- data.frame(ex.expr.norm, cl = ex.gene.cl, max = rowMaxs(ex.gene.cor))
ex.gene.plot <- ex.gene.plot %>%
  arrange(cl, -max)

ex.gene.plot2 <- melt(ex.gene.plot %>% select(-max), "cl")

ggplot(ex.gene.plot2) +
  geom_tile(aes(x = rep(1:nrow(ex.gene.plot),4), y = variable, fill = value))



## K-means clustering with common centroid seeds
# Will provide cleaner and more consistent results
centroid.seeds <- rbind(c(1,0,0,0), # Cux2 High
                        c(0,0,0,1), # Scnn1a High
                        c(0,0,1,0), # Rbp4
                        c(0,1,0,0), # Ntsr1
                        c(1,0,0,1), # Upper  
                        c(0,1,1,0), # Lower
                        c(1,1,1,0), # L4 -  
                        c(1,0,1,1)) # L6 -

# Seeded peak k-means clustering
tmp <- kmeans(ex.peak.norm, centers = centroid.seeds, iter.max = 100)
cl <- tmp$cluster
ex.peak.cl <- cl
ex.peak.cl.means <- round(do.call("rbind", 
                                  tapply(1:length(cl), 
                                         cl, 
                                         function(x){
                                           colMeans(ex.peak.norm[x,])
                                         })
                                  ),
                          digits = 2)


# heatmap of peaks in k-means
ex.peak.cor <- cor(t(ex.peak.norm), t(ex.peak.cl.means))
ex.peak.cl <- apply(ex.peak.cor, 1, which.max)

# add clusters to all.peaks
all.peaks$ex.peak.cl <- 0
all.peaks$ex.peak.cl[as.numeric(names(ex.peak.cl))] <- ex.peak.cl

# too many peaks to plot, so choose the top 500 from each cluster
ex.peak.plot <- data.frame(ex.peak.norm, cl = ex.peak.cl, max = rowMaxs(ex.peak.cor))
ex.peak.plot <- ex.peak.plot %>%
  group_by(cl) %>%
  arrange(-max) %>%
  filter(row_number() <= 500) %>%
  arrange(cl, -max)

ex.peak.plot2 <- melt(ex.peak.plot %>% select(-max), "cl")

ggplot(ex.peak.plot2) +
  geom_tile(aes(x = rep(1:nrow(ex.peak.plot),4), y = variable, fill = value))


# Seeded gene k-means clustering
tmp <- kmeans(ex.expr.norm, centroid.seeds, iter.max = 100)
cl <- tmp$cluster
ex.expr.cl.means <- round(do.call("rbind",
                                  tapply(1:length(cl), 
                                         cl, 
                                         function(x){
                                           colMeans(ex.expr.norm[x,])
                                         })
                                  ),
                          digits = 2)

# heatmap of genes in k-means
ex.gene.cor <- cor(t(ex.expr.norm), t(ex.expr.cl.means))
ex.gene.cl <- apply(ex.gene.cor, 1, which.max)

# add clusters to all.peaks
gene_to_cl <- data.frame(gene.id = names(ex.gene.cl),
                         ex.gene.cl = ex.gene.cl)
tmp <- data.frame(gene.id = all.peaks$gene.id) %>%
  left_join(gene_to_cl)
tmp[is.na(tmp)] <- 0
all.peaks$ex.gene.cl <- tmp$ex.gene.cl

ex.gene.plot <- data.frame(ex.expr.norm, cl = ex.gene.cl, max = rowMaxs(ex.gene.cor))
ex.gene.plot <- ex.gene.plot %>%
  arrange(cl, -max)

ex.gene.plot2 <- melt(ex.gene.plot %>% select(-max), "cl")

ggplot(ex.gene.plot2) +
  geom_tile(aes(x = rep(1:nrow(ex.gene.plot),4), y = variable, fill = value))




## Fisher's exact tests for enrichment of peaks in each module near genes in each module

# Select peaks that are near clustered genes
select.peak <- names(ex.peak.cl)[all.peaks[as.integer(names(ex.peak.cl))]$gene.id %in% names(ex.gene.cl)]

# Count number of peaks in each cluster that are associated with each clustered gene
tb <- table(ex.peak.cl[select.peak], ex.gene.cl[all.peaks[as.integer(select.peak)]$gene.id])

# Compute Fisher's exact test p-values and log-odds ratios (estimate)
tb.pval <- matrix(1, nrow = nrow(tb), ncol = ncol(tb))
tb.odds <- matrix(1, nrow = nrow(tb), ncol = ncol(tb))
for(i in 1:nrow(tb)){
  for(j in 1:ncol(tb)){
    # build contingency matrix (tmp.tb)
    a <- sum(tb[i,])
    b <- sum(tb[,j])
    tmp.tb <- matrix(c(sum(tb) - a - b + tb[i,j],  a - tb[i,j], b - tb[i,j], tb[i,j]), nrow = 2, ncol = 2)
    # run test
    tmp <- fisher.test(tmp.tb)
    # extract adjusted p-value from results
    pval <- tmp$p.value
    tb.pval[i,j] <- pval
    # if(tmp$estimate < 1) {
    #   tb.pval[i,j] <- -1 * tb.pval[i,j]
    # }
    tb.odds[i,j] <- tmp$estimate
  }
}

# Benjamini and Hochberg multiple hypothesis correction
tb.adj.pval <- matrix(p.adjust(tb.pval, method = "BH"),
                      nrow = nrow(tb.pval),
                      ncol = ncol(tb.pval))
# Log scale and log-odds direction
tb.log.pval <- matrix(ifelse(tb.odds > 1,-log10(tb.adj.pval),log10(tb.adj.pval)),
                      nrow = nrow(tb.pval),
                      ncol = ncol(tb.pval))

colnames(tb.log.pval) <- colnames(tb)
row.names(tb.log.pval) <- row.names(tb)

colnames(tb.odds) <- colnames(tb)
row.names(tb.odds) <- row.names(tb)

# building the heatmap for the peak-gene enrichment
pval.df <- as.data.frame(as.table(tb.log.pval))
colnames(pval.df) <- c("peak.cl", "gene.cl", "logP")
pval.df <- pval.df %>% mutate_if(is.factor, as.numeric)
#pval.df$logP[pval.df$logP > 30]  <- 30
#pval.df$logP[pval.df$logP < -20] <- -20

low_colors <- c("darkblue","blue","skyblue")
high_colors <- c("orange","orangered","red")

pval.df <- pval.df %>%
  rowwise() %>%
  mutate(fill = ifelse(abs(logP) < 2,
                       "#000000",
                       ifelse(logP > 0,
                       values_to_colors(logP, 2, 120, high_colors),
                       values_to_colors(logP, -35, 2, low_colors))))

module_pval_plot <- ggplot(pval.df) + 
  geom_tile(aes(x = gene.cl, y = peak.cl, fill = fill)) + 
  scale_y_reverse() +
  scale_fill_identity() +
  theme_classic()

ggsave("module_pval_plot.pdf",module_pval_plot,width = 4, height = 4)

odds.df <- as.data.frame(as.table(tb.odds))
colnames(odds.df) <- c("peak.cl", "gene.cl", "odds")
odds.df <- odds.df %>% mutate_if(is.factor, as.numeric) %>%
  rowwise() %>%
  mutate(fill = ifelse(log(odds,2) > 0,
                       values_to_colors(log(odds,2),0,4,high_colors),
                       values_to_colors(log(odds,2),-4,0,low_colors)))

module_odds_plot <- ggplot(odds.df) + 
  geom_tile(aes(x = gene.cl, y = peak.cl, fill = fill)) + 
  scale_fill_identity() +
  scale_y_reverse() +
  theme_classic()

ggsave("module_logodds_plot.pdf",module_odds_plot,width = 6, height = 4)

save(all.peaks,file="all_peaks.rda")

cl.comparison.pvals <- pval.df
save(ex.peak.cl, ex.gene.cl, ex.expr.cl.means, ex.peak.cl.means, cl.comparison.pvals, file = "ex_peak_gene_modules.rda")

