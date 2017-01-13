library(GenomicRanges)
library(dplyr)
library(reshape2)
library(XML)
options(StringsAsFactors = F)

load("all_peaks.rda")
load("ex_peak_gene_modules.rda")
load("ame_results.rda")
load("fimo_results.rda")
load("expr_data.rda")

all <- as.data.frame(all.peaks) %>%
  mutate(peak_gene_comb = paste(ex.peak.cl,ex.gene.cl))

# TF annotations from AnimalTFDB
tf_url <- "http://www.bioguo.org/AnimalTFDB/BrowseAllTF.php?spe=Mus_musculus"
tf <- readHTMLTable(tf_url)$table1

## NODE SELECTION
# Manually assigned from TreeFam and literature searches

nodes <- data.frame(motif_family = c("CUX","FOXP","MEF2","MEIS","NEUROD","NFIA","POU3F","RFX3","RORB","TBR1"),
                    gene.id = c("Cux1","Foxp2","Mef2c","Meis2","Neurod6","Nfia","Pou3f2","Rfx3","Rorb","Tbr1"),
                    is.repressor = c(T,T,F,F,F,F,F,F,F,F))

# For key regulators, assign each to the most strongly correlated gene cluster
ex.expr.cre <- expr.cre[,-2]
nodes.ex.expr.cre <- ex.expr.cre[nodes$gene.id,]
best_match <- numeric()
for(i in 1:nrow(nodes.ex.expr.cre)) {
  best <- which(cor(nodes.ex.expr.cre[i,],t(ex.expr.cl.means)) == max(cor(nodes.ex.expr.cre[i,],t(ex.expr.cl.means))))
  best_match <- c(best_match,best)
}

nodes <- nodes %>%
  mutate(ex.gene.cl = best_match)

# Propagate assignments to the peaks
for(i in 1:nrow(nodes)) {
  all$ex.gene.cl[all$gene.id == nodes$gene.id[i]] <- nodes$ex.gene.cl[i]
}

## EDGE SELECTION
# Significance Cutoffs
motif_peak_module_cut <- 0.01
peak_gene_module_cut <- 0.01
peak_diffbind_cut <- 0.01
gene_deseq_cut <- 0.01
activator_cor_cut <- 0.3
repressor_cor_cut <- 0.3

edge.motifs <- as.data.frame(all_fimo_gr) %>%
  left_join(all, by = "peak_id")

# Motif is siginificantly associated with a peak module?
log_mpm_cut <- -log(motif_peak_module_cut, 10)

keep_peak_ame <- ame.comparison.pvals %>%
  filter(abs(logP) > log_mpm_cut) %>%
  mutate(peak_ame_comb = paste(peak.cl,motif))

edge.motifs <- edge.motifs %>%
  mutate(peak_ame_comb = paste(ex.peak.cl, motif_family)) %>%
  filter(peak_ame_comb %in% keep_peak_ame$peak_ame_comb)

# Peak module is significantly associated with a gene module that the nearest gene is a member of?
log_pgm_cut <- -log(peak_gene_module_cut,10)

keep_peak_gene_modules <- cl.comparison.pvals %>%
  filter(logP > log_pgm_cut) %>%
  mutate(peak_gene_comb = paste(peak.cl,gene.cl))

edge.motifs <- edge.motifs %>%
  filter(peak_gene_comb %in% keep_peak_gene_modules$peak_gene_comb)

# Peak is differentially accessible among excitatory cell types?
edge.motifs <- edge.motifs %>%
  filter(ex.pPadj < peak_diffbind_cut)

# Nearest gene is differentially expressed among excitatory cell types?
edge.motifs <- edge.motifs %>%
  filter(ex.gPadj < gene_deseq_cut)

# Nearest gene is a transcription factor?
edge.motifs <- edge.motifs %>%
  filter(gene.id %in% tf$Symbol)

# node TF is a repressor?
activator_families <- nodes$motif_family[nodes$is.repressor == F]
activator.motifs <- edge.motifs %>% 
  filter(edge.motifs$motif_family %in% activator_families)

repressor_families <- nodes$motif_family[nodes$is.repressor == T]
repressor.motifs <- edge.motifs %>%
  filter(edge.motifs$motif_family %in% repressor_families)

# Accessibility and node TF are correlated/anticorrelated?
tmp <- data.frame(motif_family = activator.motifs$motif_family) %>% left_join(nodes)
activator.motifs$motif_gene <- tmp$gene.id

tmp <- data.frame(motif_family = repressor.motifs$motif_family) %>% left_join(nodes)
repressor.motifs$motif_gene <- tmp$gene.id

check_cor <- function(peaks,genes) {
  
  cors <- numeric()
  
  for(i in 1:length(peaks)) {
    gene_expr <- all %>%
      filter(gene.id == genes[i]) %>%
      select(one_of("g.cux2","g.ntsr1","g.rbp4","g.scnn1a")) %>%
      unique() %>% unlist()
    
    scaled_gene_expr <- (gene_expr - min(gene_expr))/max(gene_expr - min(gene_expr))
    
    peak_acc <- all %>%
      filter(peak_id == peaks[i]) %>%
      select(one_of("p.cux2","p.ntsr1","p.rbp4","p.scnn1a")) %>%
      unlist()
    scaled_peak_acc <- (peak_acc - min(peak_acc))/max(peak_acc - min(peak_acc))
    cors <- c(cors,cor(scaled_gene_expr,scaled_peak_acc))
  }
  
  cors
}

activator.cor <- check_cor(activator.motifs$peak_id, activator.motifs$motif_gene)
activator.motifs <- activator.motifs[activator.cor > activator_cor_cut,]

repressor.cor <- check_cor(repressor.motifs$peak_id, repressor.motifs$motif_gene)
repressor.motifs <- repressor.motifs[repressor.cor < repressor_cor_cut,]

# Collapse motifs into edges by counting the number of peaks that contain these motifs
activator.edges <- activator.motifs %>%
  group_by(motif_gene, gene.id) %>%
  summarise(n = length(unique(peak_id))) %>%
  mutate(type = "activator")

repressor.edges <- repressor.motifs %>%
  group_by(motif_gene, gene.id) %>%
  summarise(n = length(unique(peak_id))) %>%
  mutate(type = "repressor")

all.edges <- rbind(activator.edges,repressor.edges) %>% as.data.frame




# for each gene set, make sure the peaks are in clusters that match the gene clusters
select.peak.gene.df = cl.comparison.pvals[cl.comparison.pvals$logP > -log10(peak_gene_module_cut), ]
select.peak.gene.df$peak.gene.id = paste(select.peak.gene.df$peak.cl, select.peak.gene.df$gene.cl)

# make sure motif enrichment is in the correct group
tmp.df = cbind(rep(1:8,ncol(ex.ame.df)),melt(ex.ame.df))
colnames(tmp.df) = c("peak.cl","motif","logP")
select.peak.motif.df = tmp.df[tmp.df$logP > -log10(motif_peak_module_cut),]
select.peak.motif.df$peak.motif.id = paste(select.peak.motif.df$peak.cl, select.peak.motif.df$motif)

# manual setting of cluster to keep these from dropping out
ex.gene.cl["Neurod6"] = 3
ex.gene.cl["Mef2c"] = 6
ex.gene.cl["Tbr1"] = 1

# linking motifs to peaks and genes that were selected above
tmp=load("select.motif.rda")
tmp=match(select.gr$gene.id, names(ex.gene.cl))
select.gr$ex.gene.cl[!is.na(tmp)] = ex.gene.cl[tmp[!is.na(tmp)]]
tmp = setNames(rep(names(select.motif),sapply(select.motif,length)),unlist(select.motif))
select.gr$motif.TF = tmp[as.character(select.gr$motif)]
select.gr$peak.gene.id = with(values(select.gr), paste(ex.peak.cl, ex.gene.cl))
select.gr$peak.motif.id = with(values(select.gr), paste(ex.peak.cl, motif.TF))

select1 = select.gr$peak.gene.id %in% select.peak.gene.df$peak.gene.id 
select2 = select.gr$peak.motif.id %in% select.peak.motif.df$peak.motif.id
# removed in the final network: FOS and EGR1 because they have many targets
# keep only most significant peaks and genes so that the network doesn't get too cluttered
select3=!select.gr$motif.TF %in% c("FOS","EGR1") & select.gr$ex.pPadj < 0.05 & select.gr$ex.gPadj < 10^-5
# assigning motifs to real genes
# DLX/Emx1 may not appear in network because the peaks might not have passed other criteria
tmp.tf = c(CUX="Cux1",FOXP="Foxp2", MEF2="Mef2c",MEIS="Meis2",NEUROD="Neurod6",NFIA="Nfia",POU3F="Pou3f2",RFX3="Rfx3",RORB="Rorb",TBR1="Tbr1","DLX"="Emx1")
tmp.df = values(select.gr)[select1&select2&select3,c("gene.id","dist2tss","motif.TF","score","score","pvalue", "gPadj","pPadj", "ex.peak.cl","ex.gene.cl","ex.gPadj")]
tmp.df$tf = tmp.tf[as.character(tmp.df$motif.TF)]
# FOXP and CUX likely repressors 
tmp.df$repress = tmp.df$motif.TF %in% c("FOXP","CUX")
tmp.cor  =cor(t(expr.cre[tmp.df$gene.id,-2]),t(expr.cre[tmp.df$tf,-2]))
tmp1.cor  =cor(t(expr[tmp.df$gene.id,-2]),t(expr[tmp.df$tf,-2]))
tmp.df$gene.cor1 = diag(tmp.cor)
tmp.df$gene.cor2 = diag(tmp1.cor)

# segregation of Neurod1 and Neurod6: peaks that were better correlated with 1 were assigned to 1 and 6 to 6
select = tmp.df$motif.TF == "NEUROD"
neurod.cor = cor(t(expr[tmp.df$gene.id, ]),t(expr[c("Neurod1","Neurod6"),]))
select1 = select & neurod.cor[,"Neurod1"] > neurod.cor[,"Neurod6"]
tmp.df[select1,"tf"]="Neurod1"
select2 = select & neurod.cor[,"Neurod1"] < neurod.cor[,"Neurod6"]
tmp.df[select2,"tf"]="Neurod6"

tmp.df = with(tmp.df, tmp.df[repress& pmin(gene.cor1,gene.cor2) < -0.3 | !repress & pmax(gene.cor1, gene.cor2) > 0.3,])
tmp.df = tmp.df[order(tmp.df$ex.gPadj),]

gene.df = unique(values(select.gr)[,c("gene.id","ex.gene.cl","ex.gPadj")])
gene.df$is.TF = gene.df$gene.id %in% tf$Symbol
row.names(gene.df)=gene.df$gene.id
gene.df$gene.id=NULL

tmp=tapply(1:nrow(tmp.df), list(tmp.df$gene.id, tmp.df$tf), function(x){
  length(x)
})
tmp = as.data.frame(as.table(tmp))
tmp=tmp[!is.na(tmp$Freq),]
tmp = data.frame(tmp, gene.df[as.character(tmp[,1]),])
tmp.df=tmp
colnames(tmp.df)[1:2]=c("gene.id","tf")
tmp.df$repress = tmp.df$tf %in% c("Foxp2","Cux1")

# refiltering of genes: more relaxed filtering if the gene is a TF, but keep some non-TFs that are very significant.
tmp2.df = tmp.df[tmp.df$ex.gPadj < 10^-35 | tmp.df$is.TF & tmp.df$ex.gPadj < 10^-5| tmp.df$gene.id %in% tmp.tf ,]
write.csv(tmp2.df, file="select.peak.gene.csv",quote=F)

# scaling for the plots in cytoscape
tmp1.df = gene.df[row.names(gene.df)%in% c(as.character(tmp2.df$gene.id),as.character(tmp2.df$tf)),] 
tmp1.df$g.score = -log10(tmp1.df$ex.gPadj)
tmp1.df$font.size = pmax(tmp1.df$g.score , 30 * tmp1.df$is.TF)
write.csv(tmp1.df, file="select.gene.csv",quote=F)
