library(GenomicRanges)
library(dplyr)
library(readxl)
library(DiffBind)
library(ggplot2)
library(matrixStats)
source("common_functions.R")
options(stringsAsFactors = F)

load("diffbind_results.rda")
load("atac_data.rda")
load("expr_data.rda")


# Get differential peaks from diffbind contrasts
de.peaks = sapply(dbo$contrasts, 
                  function(x) {
                    tmp = x$DESeq2$de
                    tmp[which(tmp$padj < 0.05),]
                  }, 
                  simplify = F)

names(de.peaks) <- sapply(dbo$contrasts, function(x){ paste(x$name1, x$name2, sep = "_")})

# compute differentially expressed genes between each pair of cell classes (cre lines).
# This may take a while (overnight using 4 cores).
# DESeq.genes.pw() can be found in common_functions.R
de.expr <- DESeq.genes.pw(counts, expr.sample.df$cre)

# change names to match peak data
names(de.expr) <- tolower(names(de.expr))
names(de.expr) <- sub(".tg3","",names(de.expr))

save(de.expr, file = "de_gene_expr.rda")

# Get differentially expressed genes with lfc > 1 and p < 0.05 (de.g)
de.g <- unique(unlist(sapply(de.expr, 
                             function(df) {
                               row.names(df)[with(df, which(abs(lfc) > 1 & padj < 0.05))]
                             },
                             simplify = F)
                      )
               )

# Get differentially-accessible peaks (de.p)
de.p <- unique(unlist(sapply(de.peaks, 
                             function(x) {
                               x$id
                             },
                             simplify = F)
                      )
               )

# Building the all.peaks table with gene and peak values and p values

# mean peak and gene values for each cell class (gene values start with g., peak values start with p.)
tmp1 <- expr.cre
colnames(tmp1) <- paste0("g.", sub(".tg3","",tolower(colnames(tmp1))))
tmp2 <- peak.cre
colnames(tmp2) <- paste0("p.", colnames(tmp2))

values(all.peaks) <- cbind(values(all.peaks), 
                           gene.cor = gene.cor, 
                           as.data.frame(tmp1[as.character(all.peaks$gene.id),]), 
                           as.data.frame(tmp2))

# adjusted differential gene expression p-values for each pairwise comparison
gPadj <- sapply(de.expr, 
                function(x) {
                  x[all.peaks$gene.id, "padj"]
                })

gPadj[is.na(gPadj)] <- 1

# Get minimum p-value across all comparisons (gPadj)
colnames(gPadj) <- paste0("g.", colnames(gPadj))
all.peaks$gPadj <- rowMins(gPadj)

# Get minimum p-value across glutamatergic cell types (ex.gPadj)
tmp <- setdiff(colnames(gPadj), grep("gad", colnames(gPadj), value = T))
all.peaks$ex.gPadj <- rowMins(gPadj[,tmp])

for(x in colnames(gPadj)) {
  values(all.peaks)[, x] <- gPadj[,x]
}

# adjusted differential peak accessibility p-values for each pairwise comparison
pPadj <- sapply(de.peaks, 
                function(x){
                  x[match(1:length(all.peaks), x$id), "padj"]
                })

pPadj[is.na(pPadj)] <- 1

# Get minimum p-value across all comparisons (pPadj)
colnames(pPadj) <- paste0("p.", colnames(pPadj))
all.peaks$pPadj <- rowMin(pPadj)

# Get minimum p-value across glutamatergic cell types (ex.pPadj)
tmp <- setdiff(colnames(pPadj), grep("gad", colnames(pPadj), value = T))
tmp <- setdiff(tmp, grep("wt", tmp, value=T))
all.peaks$ex.pPadj <- rowMins(pPadj[,tmp])

for(x in colnames(pPadj)) {
  values(all.peaks)[, x] <- pPadj[,x]
}

save(all.peaks, file="all_peaks.rda")



