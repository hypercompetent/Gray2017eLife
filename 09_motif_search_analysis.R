library(GenomicRanges)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)
library(matrixStats)
library(dplyr)
source("common_functions.R")
options(stringsAsFactors = F)
load("tss.rda")
load("expr_data.rda")
load("all_peaks.rda")
load("ex_peak_gene_modules.rda")

all.peaks$peak_id <- 1:length(all.peaks)

# computing AME enrichment between clusters 

dir.create("motif_seq")
dir.create("ame_out")

# Retrieve sequenes for each clustered peak
# getSequence() is in common_functions.R
gr <- all.peaks[as.integer(names(ex.peak.cl))]
seq <- getSequence(gr,genome=BSgenome.Mmusculus.UCSC.mm10)
names(seq) <- paste0(seqnames(gr),":",start(gr),"_",end(gr))


# filtering the peaks based on their similarity to each cluster 
# (remove other clusters from the negative set that are too similar)
for(i in unique(ex.peak.cl)) {
  # Remove peaks whose merged width is >= 400 bp
  select <- as.vector(ex.peak.cl == i) & nchar(seq) < 400
  
  # Peak module foreground sequences
  file <- file.path("motif_seq", paste0("cl.",i,".fasta"))
  writeXStringSet(seq[select], file)
  
  pos <- ex.peak.cl.means[i,] > 0.5
  neg <- which(rowMaxs(ex.peak.cl.means[, pos, drop = F]) < 0.5)
  
  not.select <- as.vector(ex.peak.cl %in% neg) & nchar(seq) < 400
  
  # Peak module background sequences
  file <- file.path("motif_seq", paste0("cl.ex", i, ".fasta"))
  writeXStringSet(seq[not.select],file)
  
}

tmp <- sapply(sort(unique(ex.peak.cl)), 
              function(i){
                pos <- ex.peak.cl.means[i,] > 0.5
                rowMaxs(ex.peak.cl.means[,pos,drop=F]) < 0.5
              })
row.names(tmp) <- 1:8
colnames(tmp) <- row.names(tmp)
ame_contrast_groups <- tmp

# This function requires the external tool AME as well as a motif database.
ame.path <- "ame"
motif.db.path <- "motif_databases/JASPAR/JASPAR_CORE_2016.meme"

findMotif <- function(cl, ame.path, motif.db.path)
{
  file1=file.path("motif_seq", paste0("cl.",cl,".fasta"))
  file2=file.path("motif_seq", paste0("cl.ex",cl,".fasta"))
  ame <- ame.path
  cmd = paste(ame, "-o", file.path("ame_out",cl), "-control", file2, file1, motif.db.path)
  print(cmd)
  system(cmd)
  cmd = paste(ame, "-o", file.path("ame_out",paste0("ex",cl)), "-control", file1, file2, motif.db.path)
  print(cmd)
  system(cmd)
}

for(cl in 1:8) {
  findMotif(cl, ame.path, motif.db.path)
}

# parse AME results and get enrichment using AME results for comparisons between contrast groups
parse_ame <- function(fn)
{
  lines <- readLines(fn)[-(1:12)]
  tmp <- do.call("rbind",strsplit(lines, "[ \\(\\)]"))
  tb <- data.frame(tmp[,c(8,9,13,17)], stringsAsFactors = F)
  colnames(tb) <- c("ID", "TF", "pval", "padj")
  tb$pval <- as.numeric(tb$pval)
  tb$padj <- as.numeric(tb$padj)
  tb
}

ame_result <- do.call("rbind",
                      sapply(sort(unique(ex.peak.cl)), 
                             function(i) {
                               tb1 = parse_ame(file.path("ame_out", i, "ame.txt"))
                               tb2 = parse_ame(file.path("ame_out", paste0("ex",i), "ame.txt"))
                               tb1$cat = "enriched"
                               tb2$cat = "depleted"
                               tb=rbind(tb1, tb2)
                               tb$peak.cl = i
                               tb
                             },
                           simplify = F)
                      )

ame_result$logP <- -log10(ame_result$padj)
select <- ame_result$cat == "depleted"
ame_result$logP[select] <- - ame_result$logP[select]

# Selected motifs to retain for downstream analysis
select.motif = list(DLX    = "MA0882.1", 
                    NEUROD = "MA0000.1", 
                    TBR1   = c("MA0802.1", "MA0800.1"),
                    FOXP   = c("MA0593.1", "MA0481.1"), 
                    NFIA   = "MA0670.1", 
                    RORB   = c("MA0071.1", "MA0072.1"), 
                    POU3F  = c("MA0787.1", "MA0789.1", "MA0786.1", "MA0788.1"), 
                    MEF2   = c("MA0487.1", "MA0773.1", "MA0660.1", "MA0052.3"),
                    FOS    = c("MA0477.1", "MA0478.1", "MA0476.1", "MA0099.2"), 
                    EGR1   = c("MA0472.2", "MA0732.1", "MA0162.2", "MA0733.1"), 
                    MEIS   = c("MA0775.1", "MA0774.1"), 
                    RFX3   = c("MA0798.1", "MA0600.2", "MA0799.1", "MA0365.1", "MA510.2"),
                    CUX    = c("MA0755.1", "MA0679.1", "MA0757.1", "MA0756.1"))

motif_families <- data.frame(motif_id = unlist(select.motif),
                             motif_family = rep(names(select.motif),unlist(lapply(select.motif,length))))

n.cl <- sort(unique(ex.peak.cl))
ex.ame.score <- sapply(select.motif, 
                       function(x) {
                         select <- ame_result$ID %in% x
                         tmp <- with(ame_result[select,], tapply(logP, peak.cl, function(x) x[which.max(abs(x))]))
                         score <- setNames(rep(0, length(n.cl)), n.cl)
                         score[names(tmp)] <- tmp
                         score
                       })

# Some RORB comparisons returning null.
ex.ame.score[["6","RORB"]] <- 0
ex.ame.score[["8","RORB"]] <- 0


ex.ame.score <- as.matrix(ex.ame.score)
ex.ame.df <- data.frame(ex.ame.score)
ex.ame.df <- as.data.frame(lapply(ex.ame.df,unlist))

# Heatmap for enrichment of motifs per peak module
tmp.df <- data.frame(peak.cl = rownames(ex.ame.df), ex.ame.df)
tmp.df <- melt(tmp.df,"peak.cl")
colnames(tmp.df) <- c("peak.cl","motif","logP")

min.p <- min(tmp.df$logP)
max.p <- max(tmp.df$logP)

plot.df <- tmp.df %>%
  mutate_if(is.factor,as.character) %>%
  arrange(motif,peak.cl) %>%
  rowwise() %>%
  mutate(fill = values_to_colors(logP, min.p, max.p,colorset = c("darkblue","blue","skyblue","white","orange","orangered","red"))) %>%
  ungroup() %>%
  group_by(motif)

plot.df <- plot.df %>%
  ungroup() %>%
  mutate(ypos = group_indices(plot.df))

p2 <- ggplot(plot.df) + 
  geom_tile(aes(x = peak.cl, y = ypos, fill = fill),  colour = "white") + 
  scale_fill_identity() +
  scale_y_reverse()
 
p2


ame.comparison.pvals <- tmp.df

save(ame_result, ex.ame.df, ame.comparison.pvals, file = "ame_results.rda")

# FIMO to identify motif locations within all peaks

# Output all peak sequences for FIMO
dir.create("fimo")
all.seq <- getSequence(all.peaks)
names(all.seq) <- paste0(seqnames(all.peaks),":",start(all.peaks),"_",end(all.peaks))

writeXStringSet(all.seq, "fimo/all_peaks.fasta")

# This section requires the external tool FIMO as well as a motif database.
fimo.path <- "fimo"
motif.db.path <- "motif_databases/JASPAR/JASPAR_CORE_2016.meme"

motifs <- unlist(select.motif)

for(motif in motifs) {
  cmd <- paste(fimo.path, "--o", paste0("fimo/",motif), "-motif", motif, motif.db.path, "fimo/all_peaks.fasta")
  system(cmd)
}

# Parse the FIMO results as a GRanges object

all_fimo_gr <- GRanges()

for(motif in unlist(select.motif)) {
  fimo_file <- paste0("fimo/",motif,"/fimo.txt")
  fimo_out <- data.frame()
  try(fimo_out <- read.table(fimo_file))
  if(nrow(fimo_out) > 0) {
    names(fimo_out) <- c("motif_id","fasta_peak_name","rel_start","rel_end","strand","score","pval","X","motif_seq")
    fimo_out <- fimo_out %>%
      mutate(chr = sub(":.+","",fasta_peak_name),
             peak_start = as.numeric(sub(".+:","",sub("_.+","",fasta_peak_name))),
             peak_end = as.numeric(sub(".+_","",fasta_peak_name))) %>%
      mutate(start = ifelse(rel_start < rel_end,
                            peak_start + rel_start,
                            peak_start + rel_end),
             end = ifelse(rel_start < rel_end,
                          peak_start + rel_end,
                          peak_start + rel_start))
    
    fimo_gr <- GRanges(seqnames = fimo_out$chr, IRanges(start = fimo_out$start, end = fimo_out$end))
    fimo_gr$motif_id <- fimo_out$motif_id
    fimo_gr$motif_seq <- fimo_out$motif_seq
    fimo_gr$motif_pval <- fimo_out$pval
    fimo_gr$motif_family <- motif_families$motif_family[motif_families$motif_id == motif]
    fimo_gr$peak_location <- fimo_out$fasta_peak_name
    all_fimo_gr <- c(all_fimo_gr, fimo_gr)
  }
}

# Match the motifs to peaks

all.peaks_motif.hits <- findOverlaps(all.peaks, all_fimo_gr)
all.peaks_motif.hits <- as.data.frame(all.peaks_motif.hits)

missed_peaks <- all.peaks[-unique(all.peaks_motif.hits$queryHits),]

# Fully joined GRanges object

all.peaks_all.motifs <- all.peaks[all.peaks_motif.hits$queryHits]
all.peaks_all.motifs$motif_id <- all_fimo_gr$motif_id[all.peaks_motif.hits$subjectHits]
all.peaks_all.motifs$motif_pval <- all_fimo_gr$motif_pval[all.peaks_motif.hits$subjectHits]
all.peaks_all.motifs$motif_seq <- all_fimo_gr$motif_seq[all.peaks_motif.hits$subjectHits]
all.peaks_all.motifs$motif_family <- all_fimo_gr$motif_family[all.peaks_motif.hits$subjectHits]

# Summarize the results
all.peaks_all.motifs_df <- as.data.frame(all.peaks_all.motifs)

# Get the number of motifs for each family found in each peak
all.peaks_all.motifs_counts <- all.peaks_all.motifs_df %>%
  group_by(peak_id, motif_family) %>%
  summarise(n_motifs = n()) %>%
  ungroup() %>%
  mutate_if(is.factor,as.character)

family_summary <- data.frame(peak_id = 1:length(all.peaks),
                             DLX    = 0, 
                             NEUROD = 0, 
                             TBR1   = 0,
                             FOXP   = 0, 
                             NFIA   = 0, 
                             RORB   = 0, 
                             POU3F  = 0, 
                             MEF2   = 0,
                             FOS    = 0, 
                             EGR1   = 0, 
                             MEIS   = 0, 
                             RFX3   = 0,
                             CUX    = 0)

for(i in 1:nrow(all.peaks_all.motifs_counts)) {
  counts_row <- all.peaks_all.motifs_counts[i,]
  family_summary[counts_row$peak_id, counts_row$motif_family] <- counts_row$n_motifs
}

# add peak ids to the all_fimo results
all_fimo_gr$peak_id <- 0
all_fimo_gr$peak_id[all.peaks_motif.hits$subjectHits] <- all.peaks_motif.hits$queryHits

save(all_fimo_gr, all.peaks_all.motifs, family_summary, file = "fimo_results.rda")
