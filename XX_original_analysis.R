library(GenomicRanges)
library(limma)
#library(affy)
#library(genefilter)
library(gplots)
library(dplyr)
source("common_functions.R")
options(stringsAsFactors = F)
# 
# 
# jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow",   "#FF7F00", "red", "#7F0000"))
# blue.red <-colorRampPalette(c("blue", "white", "red"))
# green <-colorRampPalette(c("white", "green"))
# red <-colorRampPalette(c("white", "red"))
# blue <-colorRampPalette(c("white", "blue"))

# tmp=load("~/zizhen/mm10/refseq.exons.rda")
# tmp1= tapply(start(refseq.exons), refseq.exons$transcript_id, min)
# tmp2= tapply(end(refseq.exons), refseq.exons$transcript_id, max)
# dup = duplicated(refseq.exons$transcript_id)
# strand = setNames(as.vector(strand(refseq.exons)[!dup]), refseq.exons$transcript_id[!dup])
# chr = setNames(as.vector(seqnames(refseq.exons)[!dup]), refseq.exons$transcript_id[!dup])
# gene.name = setNames(refseq.exons$gene_id[!dup], refseq.exons$transcript_id[!dup])
# id = names(tmp1)
# transcript = GRanges(IRanges(tmp1[id], tmp2[id]), seqnames=chr[id], strand=strand[id], transcript.id=id, gene.id = gene.name[id],seqinfo=Seqinfo(paste0("chr",c(1:19,"X","Y"))))
# tss = getTss(transcript, field=c("transcript.id","gene.id"))
# save(tss, file="tss.rda")

cell.classes <- c("cux2","gad2","ntsr1","rbp4","scnn1a","mes")

# Fetch TSS locations from UCSC refGene table
# get_tss_regions() and bed_to_GRanges() are in common_functions.R
tss <- get_tss_regions(expand = 0) %>%
  unique() %>%
  filter(!grepl("_",chr)) %>%
  bed_to_GRanges()

save(tss,file = "tss.rda")

# diffbind_results were computed in 01_DiffBind_Analysis.R
load("diffbind_results.rda")

# Pull samples, peaks, and p-values from the DiffBind object (dbo)
samples <- dbo$samples
peaks <- dbo$peaks
peak.sig <- dbo$binding[,-c(1:3)]

# Convert the peak locations to a GRanges object
# chr <- factor(dbo$chrmap[dbo$binding[,1]], seqlevels(tss))
# levels(chr) <- paste0("chr",c(1:19,"X","Y"))
all.peaks <- with(dbo$binding, 
                  GRanges(IRanges(START, END), seqnames = chr))

# Calculate mean TMM values for each cell class
peak.cre <- do.call("cbind", 
                    tapply(samples$SampleID, 
                           samples$Tissue, 
                           function(x) rowMeans(peak.sig[,x])))

# Calculate mean expression values for each cell class
# Data from Tasic, et al. (2016) is retrieved in 00_retrieve_source_files.R
load("expr.rda")
expr.sample.df$cre <- tolower(as.character(expr.sample.df$cre))
expr.sample.df$cre[expr.sample.df$cre=="scnn1a-tg3"] <- "scnn1a"
select <- expr.sample.df$cre %in% colnames(peak.cre)
expr.sample.df <- expr.sample.df[select,]
counts < counts[,select]
expr <- log2(counts + 1)

# Get means from scRNA-Seq for each cell class.
expr.cre = do.call("cbind",
                   tapply(row.names(expr.sample.df), 
                          expr.sample.df$cre,
                          function(x) rowMeans(counts[,x])))

# Filter expression and TSS locations for genes that are found both in the expression data
# and in the TSS table
common.genes <- intersect(row.names(expr.cre),tss$gene.id)
expr.cre <- expr.cre[common.genes,]
expr.cre <- log2(expr.cre + 1)
tss <- tss[tss$gene.id %in% common.genes]

# Get distance from peaks to the nearest TSS
tmp = getNearestDist(all.peaks,tss,return.nearest=T)
all.peaks$gene.id = tmp$gene.id
all.peaks$dist2tss = tmp$dist
save(all.peaks, peak.sig, sample.df, file="atac.data.rda")

# correlation of each peak to nearest gene across the Cre lines.
vecCor <- function(A,B)
  {
    cA <- A - rowMeans(A)
    cB <- B - rowMeans(B)
    sA <- sqrt(rowMeans(cA^2))
    sB <- sqrt(rowMeans(cB^2))
    rowMeans(cA * cB) / (sA * sB)
  }
###correlation with expression
gene.cor <- vecCor(peak.cre, expr.cre[all.peaks$gene.id,])
gene.cor[is.na(gene.cor)] <- 0
gene.cor.df <- data.frame(cor = gene.cor)

# Permutation of correlation to find cutoffs
perm.cor= sapply(1:10, function(i){
  p = setNames(sample(common.genes), common.genes)
  gene.cor = vecCor(peak.cre, expr.cre[p[all.peaks$gene.id],])
})
perm.cor[is.na(perm.cor)] = 0
colnames(perm.cor) = paste0("perm", 1:10)
perm.cor.df = as.data.frame(as.table(perm.cor))
colnames(perm.cor.df) = c("site", "cat", "cor")
cor.df = rbind(perm.cor.df, data.frame(site=names(gene.cor), cat=rep("obs", length(gene.cor)), cor=gene.cor))
g=ggplot(cor.df,aes(x=cor,color=cat, fill=cat,alpha=0.5))+ geom_density() + scale_fill_manual(values=c(rep("blue",10), "red"))
pdf("cor.density.pdf")
g
dev.off()
# Cutoff of 0.5 looks reasonable

# Get differential peaks from diffbind contrasts
de.peaks = sapply(dbo$contrasts, function(x){
  tmp=x$DESeq2$de
  tmp[which(tmp$padj < 0.05),]
},simplify=F)
names(de.peaks)  =  sapply(dbo$contrasts, function(x){paste(x$name1, x$name2,sep="_")})

# de genes from DESeq analysis of V1 dataset (Tasic 2016)
source("~/zizhen/My_R/de.genes.R")
counts = counts[common.genes,]
de.expr = DESeq.genes.pw(counts, expr.sample.df$cre)
tmp = DESeq.genes.pw(counts, as.factor(expr.sample.df$cre!="Gad2"))
# Gad2 vs excitatory
names(tmp) = c("Gad2_.Gad2")
tmp1= c(de.expr, tmp)
de.expr = tmp1
save(de.expr, file="DESeq.de.df.rda")

# Get differentially expressed genes with lfc > 1 and p < 0.05
de.df = de.expr
#de.df = DE.genes.pw(expr[common.genes,], expr.sample.df$cre)
de.g = unique(unlist(sapply(de.expr, function(df){
  row.names(df)[with(df, which(abs(lfc) > 1 & padj < 0.05))]
},simplify=F)))
# Get list of unique DE peaks across all contrasts
de.p = unique(unlist(sapply(de.peaks, function(x)x$id,simplify=F)))

# Building the all.peaks table with gene and peak values and p values
tmp1 = expr.cre
colnames(tmp1)=paste0("g", colnames(tmp1))
tmp2 = peak.cre
colnames(tmp2)=paste0("p", colnames(tmp2))
values(all.peaks) = cbind(values(all.peaks),gene.cor= gene.cor, as.data.frame(tmp1[as.character(all.peaks$gene.id),]), as.data.frame(tmp2))
gPadj = sapply(de.df, function(x){
  x[all.peaks$gene.id, "padj"]
})
gPadj[is.na(gPadj)]=1
colnames(gPadj)=paste0("g", colnames(gPadj))
all.peaks$gPadj = rowMins(gPadj)
tmp = setdiff(colnames(gPadj), grep("Gad",colnames(gPadj),value=T))
all.peaks$ex.gPadj = rowMins(gPadj[,tmp])

for(x in colnames(gPadj)){
  values(all.peaks)[, x] = gPadj[,x]
}

pPadj = sapply(de.peaks, function(x){
  x[match(1:length(all.peaks), x$id), "padj"]
})
pPadj[is.na(pPadj)] = 1
colnames(pPadj)=paste0("p", colnames(pPadj))
all.peaks$pPadj = rowMin(pPadj)
values(all.peaks) = cbind(values(all.peaks),as.data.frame(gPadj), as.data.frame(pPadj))
save(all.peaks, file="all.peaks.rda")
# CSV of peaks that are significant for both differential accessibility and expression
# pPadj and gPadj are most significant value from all comparisons
tmp.peaks = all.peaks[all.peaks$pPadj < 0.01 & all.peaks$gPadj < 0.01]
write.table(as.data.frame(tmp.peaks), file="de.peaks.csv", sep=",", quote=F,row.names=F)


###findMotifOccurences
# First pass with Zizhen's code.

tmp=load("~/tool/motif_databases/jaspar.rda")
jaspar = all.pwms
grep("FOXP2", names(jaspar),value=T)
grep("RFX", names(jaspar),value=T)
grep("RORA", names(jaspar))
select.motif = c("MA0658.1 LHX6","MA0048.2 NHLH1","MA0461.2 Atoh1","MA0802.1 TBR1","MA0593.1 FOXP2","MA0670.1 NFIA", "MA0071.1 RORA", "MA0787.1 POU3F2","MA0052.3 MEF2A","MA0476.1 FOS", "MA0162.2 EGR1","MA0798.1 RFX3")
library(BSgenome.Mmusculus.UCSC.mm10)
comb.seq = getSequence(all.peaks)
source("~/zizhen/My_R/motifRG/R/pwm.model.R")

pm <- sapply(select.motif, function(x){
  m=matchPWMCoor(logPWM(jaspar[[x]]), comb.seq, all.peaks,min.score=10)
  m$nscore = m$score/max(m$score)
  m
},simplify=F)

for(x in names(pm)){
  pm[[x]]$tf = x
}
load("comb.cov.rda")
all.peaks.summit = findSummits(all.peaks, combCoverage(comb.cov))
save(all.peaks.summit, file="all.peaks.summit.rda")

tmp = unlist(GRangesList(pm))
tmp$dist2summits = getNearestDist(tmp, all.peaks.summit)
values(tmp) = cbind(values(tmp), values(all.peaks)[tmp$seq.id,])

motif.match = tmp
row.names(motif.match)= NULL
names(motif.match)= NULL
save(motif.match, file="motif.match.rda")
tmp = as.data.frame(motif.match)
tmp$tf = gsub("M[^ ]* (.*)$", "\\1", tmp$tf)
###only select the top candidates
select = with(tmp, gene.cor > 0.5 & gPadj < 10^-3 & pPadj < 10^-3  & abs(dist2summits) < 100 & score > 0.7)
write.table(tmp[select,], file="motif.match.csv", sep=",", row.names=F,quote=F)

write.table(tmp[,c(1:3, 10, 7,5)], sep="\t", col.names=F, row.names=F, quote=F, file="motif.match.bed")



###co-cluster peaks and genes
# May not have been used
select = with(values(all.peaks), pPadj < 10^-3 & rowMaxs(peak.cre) - rowMins(peak.cre) > 3)
cols = c("Gad2", "Ntsr1","Rbp4" ,"Scnn1a", "Cux2")
peak.cre= peak.cre[,cols]
peak.norm = peak.cre[select,] - rowMins(peak.cre[select,])
peak.norm = peak.norm/rowMaxs(peak.norm)
hc  = hclust(dist(t(peak.cre)),method="average")
pdf("peak.hc.pdf")
plot(hc)
dev.off()

hc  = hclust(dist(t(expr.cre)),method="average")
pdf("expr.hc.pdf")
plot(hc)
dev.off()


 
tmp=kmeans(peak.norm, 13)
#tmp=kmeans(peak.norm, rbind(peak.cl.means, rbind(c(0,0,1,0,0),c(0,0,0,1,1))))
tmp=kmeans(peak.norm, rbind(peak.cl.means, rbind(c(0,0,1,1,0))),iter.max=100)
cl = tmp$cluster
peak.cl = cl
peak.cl.means = round(do.call("rbind",tapply(1:length(cl), cl, function(x){
  colMeans(peak.norm[x,])
})),digits=2)
n.cl = length(unique(cl))
cl.de <- matrix(0, nrow=n.cl, ncol=n.cl)
for(i in 2:n.cl){
  for(j in 1:(i-1)){
    cl.de[i,j] = sum(abs(peak.cl.means[i,] - peak.cl.means[j,]) > 0.5)
    cl.de[j,i] = cl.de[i,j]
  }
}
rowSums(cl.de==0)
peak.cl.means = do.call("rbind",tapply(1:length(cl), cl, function(x){
  colMeans(peak.norm[x,])
}))
peak.cl.means = round(peak.cl.means[,cols],digits=2)





cl.df = as.data.frame(table(cl))
cl.df$cl = as.integer(cl.df$cl)
g=ggplot(cl.df,aes(x=cl, y=Freq)) + geom_bar(stat="identity") + scale_x_continuous(trans="reverse", breaks=unique(tmp.df$cl))+ coord_flip()
pdf("cl.size.pdf")
g
dev.off()
peak.cl = cl
save(peak.cl, file="peak.cl.rda")

# Pulling genes near the peaks in the peak clusters to see if they follow similar patterns to peaks?
gene.cl.df = do.call("rbind", sapply(levels(cl), function(i){
  g = unique(all.peaks$gene.id[as.integer(names(cl)[cl==i])])
  g = intersect(g, de.g)
  tmp.dat = expr.cre[g,]
  tmp.dat = tmp.dat - rowMeans(tmp.dat)
  tmp.dat = tmp.dat / rowMax(tmp.dat)
  tmp.df = as.data.frame(as.table(tmp.dat))
  tmp.df$cl = i
  tmp.df
},simplify=F))


gene.cl.df$cl = factor(gene.cl.df$cl, levels=levels(cl))
gene.cl.df$Var2 = factor(as.character(gene.cl.df$Var2), levels=cols)
g= ggplot(gene.cl.df, aes(Var2, Freq, color=cl)) +geom_boxplot() + facet_grid(cl~.) + ylim(0.5, 1)
g = g+   theme( axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none",legend.background=element_blank(),panel.grid.major=element_blank(),panel.background=element_blank())
pdf("peak.cl.gene.pdf")
g
dev.off()

####cluster genes
# only DEGenes for clustering
de.g = intersect(de.g, row.names(expr.cre))
tmp.dat = expr.cre[de.g,cols]
tmp.dat = tmp.dat - rowMins(tmp.dat)
expr.norm = tmp.dat / rowMaxs(tmp.dat)

# checks clustering of genes (adding 01110 to see if that pattern can be found among genes)
# done after separate clustering to see if patterns that seem to be missing can be obtained by explicitly defining these means
tmp=kmeans(expr.norm, centers=rbind(expr.cl.means,c(0,1,1,1,0)),iter.max=100)
#tmp=kmeans(expr.norm, centers=expr.cl.means[-6,],iter.max=100)
cl = tmp$cluster
expr.cl.means = round(do.call("rbind",tapply(1:length(cl), cl, function(x){
  colMeans(expr.norm[x,])
})),digits=2)
gene.cl=cl
n.cl = length(unique(cl))
cl.de <- matrix(0, nrow=n.cl, ncol=n.cl)
# defining how different each cluster is from each other cluster to build contrast groups for comparison
for(i in 2:n.cl){
  for(j in 1:(i-1)){
    cl.de[i,j] = sum(abs(expr.cl.means[i,] - expr.cl.means[j,]) > 0.45)
    cl.de[j,i] = cl.de[i,j]
  }
}
rowSums(cl.de==0)
cl.hc = hclust(dist(expr.cl.means))
pdf("tmp.pdf")
plot(cl.hc)
dev.off()

gene.cl = cl
save(gene.cl, peak.cl, file="gene.peak.cl.rda")
# Reordering for consistency
peak.cl = setNames(factor(as.integer(peak.cl), levels=c(1,2,4,3,7,5,6,8,10, 9, 11,12)), names(peak.cl))
gene.cl = setNames(factor(as.integer(gene.cl), levels=c(8,9,11,4,7,6,1,3,2,10,5)), names(gene.cl))
levels(gene.cl) = 1:length(levels(gene.cl))
levels(peak.cl) = 1:length(levels(peak.cl))
save(gene.cl, peak.cl, file="gene.peak.cl.rda")

###compute the overlap of de
# comparisons of peak clusters and gene clusters to see if they match
# Fisher test for enrichment of each peak cluster near genes in each gene cluster
select.peak = names(peak.cl)[all.peaks[as.integer(names(peak.cl))]$gene.id %in% names(gene.cl)]
tb= table(peak.cl[select.peak], gene.cl[all.peaks[as.integer(select.peak)]$gene.id])
tb.pval = matrix(1, nrow=nrow(tb), ncol=ncol(tb))
for(i in 1:nrow(tb)){
  for(j in 1:ncol(tb)){
    a = sum(tb[i,])
    b = sum(tb[,j])
    tmp.tb = matrix(c(sum(tb)- a -b + tb[i,j],  a- tb[i,j], b-tb[i,j], tb[i,j]), nrow=2, ncol=2)
    tmp=fisher.test(tmp.tb)
    tb.pval[i,j] = -log10(tmp$p.value)
    if(tmp$estimate < 1){
      tb.pval[i,j] = - tb.pval[i,j]
    }
  }
}

colnames(tb.pval) = levels(gene.cl)
row.names(tb.pval) =row.names(tb)

tmp.df = as.data.frame(as.table(tb.pval))
colnames(tmp.df) = c("peak.cl", "gene.cl", "logP")
tmp.df$gene.cl = factor(as.integer(tmp.df$gene.cl),levels=length(levels(gene.cl)):1 )
tmp.df$logP[tmp.df$logP > 30] = 30
tmp.df$logP[tmp.df$logP < -20] = -20
p <- ggplot(tmp.df, aes(peak.cl, gene.cl)) + geom_tile(aes(fill = logP),  colour = "white") + scale_fill_gradient2(low = "blue",mid="white",high = "red") 
pdf("peak.gene.cl.overlap.pdf")
p
dev.off()



# Showing gene profiles for each cluster
cl = gene.cl
tmp.df= as.data.frame(as.table(expr.norm[names(cl),]))
tmp.df$cl = cl[as.character(tmp.df$Var1)]
tmp.df$Var2 = factor(as.character(tmp.df$Var2), levels=cols)
g= ggplot(tmp.df, aes(Var2, Freq,group= Var1, alpha=0.1, color=cl)) +geom_line() + facet_grid(cl~.) + scale_color_manual(values=jet.colors(length(levels(gene.cl))))
g = g+   theme( axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none",legend.background=element_blank(),panel.grid.major=element_blank(),panel.background=element_blank())
g1 = g + stat_summary(aes(y=Freq, group=Var2),fun.y="mean", colour="black", geom="point") 
pdf("gene.cl.13.pdf",height=10,width=4)
g1
dev.off()

# showing peak profiles for each cluster (top 500 most strongly correlated peaks with the mean of the cluster)
cl = peak.cl
peak.cor = cor(t(peak.norm), t(peak.cl.means))
select = row.names(peak.norm)[rowMaxs(peak.cor) > 0.8 & apply(peak.cor, 1, which.max) == cl]
n.cl = length(unique(cl))
select.peak=sapply(1:n.cl, function(i){
  tmp.select= select[cl[select] == i]
  tmp.select=tmp.select[head(order(peak.cor[tmp.select, i],decreasing=T), 500)]
},simplify=F)
select.peak = unique(unlist(select.peak))
tmp.df= as.data.frame(as.table(peak.norm[select.peak,]))
tmp.df$cl = factor(cl[as.character(tmp.df$Var1)])
tmp.df$Var2 = factor(as.character(tmp.df$Var2), levels=cols)
g= ggplot(tmp.df, aes(Var2, Freq,group= Var1, alpha=0.1, color=cl)) +geom_line() + facet_grid(cl~.) + scale_color_manual(values=jet.colors(length(levels(peak.cl))))
g = g+   theme( axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none",legend.background=element_blank(),panel.grid.major=element_blank(),panel.background=element_blank())
g1 = g + stat_summary(aes(y=Freq, group=Var2),fun.y="mean", colour="black", geom="point")
pdf("peak.cl.14.pdf",height=10,width=4)
g1
dev.off()

# plotting heatmap of the peaks that are in the peak clusters
peak.cor = cor(t(peak.norm), t(peak.cl.means))
peak.cl = apply(peak.cor, 1, which.max)
peak.cl = setNames(factor(peak.cl, level=c(1:4,13, 6,8:12,14,5,7)),names(peak.cl))
ord = order(peak.cl, order(rowMaxs(peak.cor[names(peak.cl),]),decreasing=T))
peak.col = jet.colors(length(unique(peak.cl)))[peak.cl]
pdf("peak.cl.heatmap.pdf",height=12, width=5)
heatmap.2(peak.norm[ord,], Rowv=NULL, Colv=NULL, cexRow=0.1,trace="none",RowSideColors=peak.col[ord],col=blue.red(100))
dev.off()

# plotting heatmap of genes that are in the gene clusters
gene.cor = cor(t(expr.norm), t(expr.cl.means))
gene.cl = apply(gene.cor, 1, which.max)
ord = order(gene.cl, order(rowMaxs(gene.cor[names(gene.cl),]),decreasing=T))
gene.col = jet.colors(length(unique(gene.cl)))[gene.cl]
pdf("gene.cl.heatmap.pdf",height=12, width=5)
heatmap.2(expr.norm[ord,], Rowv=NULL, Colv=NULL, cexRow=0.1,trace="none",RowSideColors=gene.col[ord],col=blue.red(100))
dev.off()

save(peak.cl, gene.cl, expr.cl.means, peak.cl.means, file="gene.peak.cl.rda")


# replacing previous results with FIMO searches for motifs
# searched using FIMO with TOMTOM database (may be JASPAR CORE 2014)
load("comb.cov.rda")
source("~/zizhen/My_R/misc.R")
tb=read.table("fimo_out/motif.csv",sep=",",colClasses=c("character","character","integer","integer","integer","integer","numeric","character","character","character"))
tb[,8] = as.numeric(tb[,8])
tb[,9] = as.numeric(tb[,9])
tmp1 = as.integer(tb[,3]) + tb[,5] - 1
tmp2 = as.integer(tb[,3]) + tb[,6] - 1
start = pmin(tmp1, tmp2)
end = pmax(tmp1, tmp2)
strand = rep("+", nrow(tb))
strand[tmp2 < tmp1] = "-"
gr = GRanges(IRanges(start,end), seqnames=tb[,2],strand=strand)
values(gr)= data.frame(motif=tb[,1], score = tb[,7], pvalue=tb[,8],qvalue=tb[,9],pattern=tb[,10])
tmp=do.call("rbind", strsplit(names(jaspar), " "))
tf = setNames(tmp[,2], tmp[,1])
gr$tf = tf[as.character(gr$motif)]
gr$pattern = as.character(gr$pattern)
save(gr, file="motif.match.gr.rda")

# These motifs came from pairwise comparisons in the diffbind peaks using AME (see findPairMotif)
# collapsing similar motifs by family
select.motif = list(DLX="MA0882.1", NEUROD="MA0000.1", TBR1=c("MA0802.1","MA0800.1"),FOXP=c("MA0593.1","MA0481.1"), NFIA="MA0670.1", RORB=c("MA0071.1", "MA0072.1"), POU3F=c("MA0787.1","MA0789.1", "MA0786.1", "MA0788.1"), MEF2=c("MA0487.1","MA0773.1","MA0660.1","MA0052.3"),FOS=c("MA0477.1","MA0478.1","MA0476.1","MA0099.2"), EGR1=c("MA0472.2","MA0732.1","MA0162.2","MA0733.1"), MEIS = c("MA0775.1", "MA0774.1"), RFX3=c("MA0798.1","MA0600.2","MA0799.1","MA0365.1","MA510.2"),CUX=c("MA0755.1","MA0679.1","MA0757.1","MA0756.1"))


merged.gr = sapply(names(select.motif), function(x){
  print(x)
  select.gr = gr[as.character(gr$motif) %in% select.motif[[x]]]
  merged.gr=reduce(select.gr)
  tmp=findOverlaps(select.gr, merged.gr, select="first")
  gr.padj = values(select.gr)$pvalue
  select = tapply(1:length(select.gr), tmp, function(x){x[which.min(gr.padj[x])]})
  values(merged.gr) = values(select.gr)[select,]
  merged.gr$merged.tf = x
  merged.gr
},simplify=F)
select.gr = unlist(GRangesList(merged.gr))

# assigning motifs to the peaks
tmp = getNearestDist(select.gr, all.peaks.summit,return.nearest=T)
select.gr$peak.id = tmp$map.id
select.gr$dist2summits=tmp$dist
values(select.gr) = cbind(values(select.gr), values(all.peaks)[select.gr$peak.id,])
save(select.gr, file="select.motif.rda")


###compare limma and DESeq
# Lucas can ignore
load("all.peaks.deseq.rda")
tmp = all.peaks
load("all.peaks.rda")
xyplot(-log10(tmp$gPadj) ~ (-log10(all.peaks$gPadj)))


load("/data/mct-t200/T502/atac_analysis/diffbind_downsampled_3.2M/peaks_all_samples/deseq2_padj.RData")
for(x in colnames(deseq2_results[-c(1:3)])){
  values(all.peaks)[,x] = deseq2_results[,x]
}

for(x in colnames(gPadj)){
  values(all.peaks)[,x] = gPadj[,x]
}
# pairwise comparisons of peaks to genes and gene to peaks
processPair <- function(p)
{
  print(p)
  pair= unlist(strsplit(p, "_"))
  peak.df = de.peaks[[p]]
  peak.df = cbind(peak.df, peak.cre[peak.df$id,pair])
  df = de.df[[p]]
  df = de.df[[p]]
  df$padj[is.na(df$padj)]=1
  df$lfc[is.na(df$lfc)]=0
  df = df[order(df$padj),]
  g1 = row.names(df)[with(df, which(lfc > 1 & padj < 0.05))]
  g2 = row.names(df)[with(df, which(lfc < -1 & padj < 0.05))]
    
  tmp = match(all.peaks$gene.id, g1)
  ord = order(tmp)
  p1 = ord[all.peaks$gene.id[ord] %in% g1]
  tmp = match(all.peaks$gene.id, g2)
  ord = order(tmp)
  p2 = ord[all.peaks$gene.id[ord] %in% g2]
  p.gcat = rep("", length(all.peaks))
  p.gcat[p1] = "p1"
  p.gcat[p2] = "p2"
  peak.df = peak.df[order(peak.df$padj),]
  de.p1 = with(peak.df, id[padj < 0.01 & peak.df[[pair[1]]] > peak.df[[pair[2]]]+1])
  de.p2 = with(peak.df, id[padj < 0.01 & peak.df[[pair[1]]]+1 < peak.df[[pair[2]]]]) 
  pcat = rep("", length(all.peaks))
  pcat[de.p1] = "p1"
  pcat[de.p2] = "p2"

  table(p.gcat, pcat)
  ####compute the overlap
  df$gsig = -log10(df$padj)
  neg = df$lfc < 0
  df$gsig[neg]  = - df$gsig[neg]
  
  tmp.df = data.frame(gene.id=all.peaks$gene.id[peak.df$id], peak.df, stringsAsFactors=F)
  tmp.df= cbind(tmp.df, df[as.character(tmp.df$gene.id),])
  tmp.df$gsig[tmp.df$gsig > 10] = 10
  tmp.df$gsig[tmp.df$gsig < -10] = -10
  
  p.g1= sum(p.gcat=="p1")
  p.g2= sum(p.gcat=="p2")
  p1.g1= sum(pcat=="p1" & p.gcat=="p1")
  p2.g2= sum(pcat=="p2" & p.gcat=="p2")
  p.pval1 = fisher.test(table(pcat=="p1", p.gcat=="p1"))$p.value
  p.pval2 = fisher.test(table(pcat=="p2", p.gcat=="p2"))$p.value
  
  g.p1 = unique(all.peaks$gene.id[pcat=="p1"])
  g.pval1 = fisher.test(table(row.names(df) %in% g1, row.names(df) %in% g.p1))$p.value
  g1.p1 = length(intersect(g1, g.p1))
  g.p1 = length(g.p1)
  
  g.p2 = unique(all.peaks$gene.id[pcat=="p2"])
  g.pval2 = fisher.test(table(row.names(df) %in% g2, row.names(df) %in% g.p2))$p.value
  g2.p2 = length(intersect(g2, g.p2))
  g.p2 = length(g.p2)
  
  
  stats=list(g1=length(g1), p1=length(de.p1), p.g1=p.g1, p1.g1=p1.g1,p.pval1=p.pval1, g2=length(g2),p2=length(de.p2), p.g2=p.g2,p2.g2=p2.g2,p.pval2=p.pval2, g.p1=g.p1, g1.p1 = g1.p1, g.pval1=g.pval1,  g.p2 = g.p2, g2.p2=g2.p2, g.pval2=g.pval2)
  g=ggplot(tmp.df, aes_string(x=pair[1], y=pair[2])) + geom_point(aes(color=gsig),size=1)+ scale_colour_gradient2(high="red", low="blue",midpoint=0)+xlim(0,8) + ylim(0,8)
  label1 = paste0("p1 =", stats$p1,"\n",
    "%p1.g1 =", format(stats$p1.g1/stats$p1,digits=2),"\n",
    "p.pval1 =", format(stats$p.pval1, digits=2),"\n",
    "g1 =", stats$g1,"\n",
    "%g1.p1 =", format(stats$g1.p1/stats$g1,digits=2),"\n",
    "g.pval1 =", format(stats$g.pval1, digits=2),"\n")
  
  label2 = paste0("p2 =", stats$p2,"\n",
    "%p2.g2 =", format(stats$p2.g2/stats$p2,digits=2),"\n",
    "p.pval2 =", format(stats$p.pval2, digits=2),"\n",
    "g2 =", stats$g2,"\n",
    "%g2.p2 =", format(stats$g2.p2/stats$g2,digits=2),"\n",
    "g.pval2 =", format(stats$g.pval2, digits=2),"\n")
  
  g = g+ annotate("text", x=5.8, y=1.2,  label=label1, hjust=0)
  g = g+ annotate("text", x=0, y=7,  label=label2, hjust=0)
  pdf(paste0(p, ".pw.pdf"))
  print(g)
  dev.off()
  return(list(de=tmp.df, p.gcat =p.gcat, pcat=pcat, stats=stats))
}

de.result <- sapply(names(de.df), processPair,simplify=F)
de.result.df = t(sapply(de.result, function(x)x$stats))
save(de.result.df , file="de.result.df.rda")
save(de.result, file="de.result.rda")
  
select.peak = as.integer(names(peak.cl))
gene.peak.df = data.frame(peak=select.peak, peak.cl = peak.cl, gene = all.peaks$gene.id[select.peak])
gene.peak.df$gene.cl = gene.cl[as.character(gene.peak.df$gene)]

tmp = all.peaks[gene.peak.df$peak,]       
values(tmp) = cbind(values(tmp),gene.peak.df[,-3])
select.peak = tmp
save(select.peak, file="select.peak.rda")

# Getting the FIMO results and sequences of the motifs
library(BSgenome.Mmusculus.UCSC.mm10)
gr = all.peaks[select.peak$peak]
seq = getSequence(gr)
names(seq)= paste0(seqnames(gr),":",start(gr),"_",end(gr))
for(i in levels(peak.cl)){
  # exclude peaks that are too large
  select = as.vector(select.peak$peak.cl == i) & width(gr) < 400
  file=file.path("motif_seq", paste0("cl.",i,".fastq"))
  writeXStringSet(seq[select],file)
}

# pairwise comparisons of all peak clusters
write.table(do.call("rbind",sapply(2:length(levels(peak.cl)),function(i){
  cbind(1:(i-1),rep(i, i-1))
},simplify=F)),file="peak.cl.pairs",row.names=F)

load("select.motif.rda")
# compute motif cluster enrichment (not using the cluster contrast table)
tmp=findOverlaps(select.gr, select.peak, select="first")
tmp.tb=table(select.gr$tf, tmp)
tmp.mat = matrix(0, nrow=length(select.peak),ncol=nrow(tmp.tb))
colnames(tmp.mat) = row.names(tmp.tb)
tmp.mat[as.integer(colnames(tmp.tb)),]= t(tmp.tb)
values(select.peak) = cbind(values(select.peak), as.data.frame(tmp.mat))
write.table(as.data.frame(select.peak), file="peak.gene.motif.cl.csv",sep=",",row.names=F)


tmp = sort(subsetByOverlaps(gr, all.peaks[177110,]))
tmp = as.data.frame(tmp)[,c(1:3,11,7,10)]

tf = read.table("~/zizhen/mm10/mm10.tf.csv", header=T,sep=",")

# finding DE TFs
cre = setNames(factor(expr.sample.df$cre,levels=cols), row.names(expr.sample.df))
top.g = selectMarkers(expr, cre, de.df,  n.markers=300, q=0.65,lfc=1.5)
de.tf = intersect(top.g, tf$Symbol)
tmp.dat = as.matrix(expr[de.tf,])
tmp.dat = tmp.dat - rowMeans(tmp.dat)
cre.col = jet.colors(length(levels(cre)))[cre]
ord = order(cre, expr.sample.df$cl)
gene.hc = hclust(dist(tmp.dat),method="ward")
pdf("de.tf.cre.pdf")
heatmap.3(tmp.dat[,ord],Colv=NULL, Rowv= as.dendrogram(gene.hc),trace="none",ColSideColors=cre.col[ord],col=blue.red(100) )
dev.off()

# refining TF set
select.de.tf = c("Dlx5","Nr2f2","Maf","Tcf4","Pbx1","Neurod2","Satb1","Satb2","Tbr1","Klf10","Neurod6", "Foxp1","Bcl11b","Fezf2","Meis2","Foxp2","Nfia","Atf6","Nfib","Id2","Tox","Etv1","Rorb", "Bhlhe40","Pou3f2","Nr4a1", "Egr1","Fos","Cux1", "Tshz1", "Neurod1")
tmp.dat = as.matrix(expr[select.de.tf,])
tmp.dat = tmp.dat - rowMeans(tmp.dat)
cre.col = jet.colors(length(levels(cre)))[cre]
gene.hc = hclust(dist(tmp.dat),method="ward")
color.df = read.table("/data/mct-t200/Vignette/vignette_clusters_colors_v2.csv", sep=",", header=T,fill=T,comment.char="",stringsAsFactors=F)
tmp.col = setNames(color.df$vignette_color,color.df$cluster_id)
expr.sample.df$col = tmp.col[as.character(expr.sample.df$cl)]
cl.col = expr.sample.df$col
ord = order(cre, expr.sample.df$cl)
pdf("de.tf.cre.pdf")
heatmap.3(tmp.dat[,ord],Colv=NULL, Rowv= NULL,trace="none",ColSideColors=rbind(cl.col,cre.col)[,ord],col=blue.red(100) )
dev.off()
head(expr.sample.df[ord,c("cre","cl","col")])



# 1 generation without centers provided, then find cleaner trend using patterns

###only cluster excitatory peaks
# as before, but using only the excitatory data
center = rbind(c(1,0,0,0), c(0,1,0,0), c(0,0,1,0), c(0,0,0,1), c(1,1,0,0),  c(0,0,1,1), c(0,1,1,1), c(1,1,0,1))

p.cols = grep("^p",colnames(values(all.peaks)),value=T)
p.cols= setdiff(p.cols, c(grep("Gad2", p.cols, value=T),c(grep("mES", p.cols, value=T))))
p.cols = grep("_",p.cols, value=T)

ex.pPadj = rowMins(as.matrix(values(all.peaks)[,p.cols]))
all.peaks$ex.pPadj = ex.pPadj
# filtered to get significantly differentially accessible peaks between excitatory types
# four-fold differentially accessible
select = ex.pPadj < 10^-2 & rowMaxs(peak.cre[,-2]) - rowMins(peak.cre[,-2]) > 2
ex.peak.cre= peak.cre[,cols[-1]]
ex.peak.norm = ex.peak.cre[select,] - rowMins(ex.peak.cre[select,])
ex.peak.norm = ex.peak.norm/rowMaxs(ex.peak.norm)
 
tmp=kmeans(ex.peak.norm, center,iter.max=100)
tmp=kmeans(ex.peak.norm,12,iter.max=100)
cl = tmp$cluster
ex.peak.cl = cl
ex.peak.cl.means = round(do.call("rbind",tapply(1:length(cl), cl, function(x){
  colMeans(ex.peak.norm[x,])
})),digits=2)

####FINDME

# heatmap of peaks in ex.clusters
ex.peak.cor = cor(t(ex.peak.norm), t(ex.peak.cl.means))
ex.peak.cl = apply(ex.peak.cor, 1, which.max)
ord = order(ex.peak.cl, order(rowMaxs(ex.peak.cor[names(ex.peak.cl),]),decreasing=T))
ex.peak.col = jet.colors(length(unique(ex.peak.cl)))[ex.peak.cl]
pdf("ex.peak.cl.8.heatmap.pdf",height=12, width=5)
heatmap.2(ex.peak.norm[ord,], Rowv=NULL, Colv=NULL, cexRow=0.1,trace="none",RowSideColors=ex.peak.col[ord],col=blue.red(100))
dev.off()

tmp = names(de.expr)
tmp=setdiff(tmp, grep("Gad",tmp, value=T))
# picking significant genes, excluding gad2
ex.de.g = unique(unlist(sapply(de.expr[tmp], function(df){
  row.names(df)[with(df, which(abs(lfc) > 1 & padj < 0.05))]
},simplify=F)))
ex.de.g = intersect(ex.de.g, row.names(expr.cre))
tmp.dat = expr.cre[ex.de.g,cols[-1]]
tmp.dat = tmp.dat - rowMins(tmp.dat)
ex.expr.norm = tmp.dat / rowMaxs(tmp.dat)

# clustering of excitatory DEGenes and heatmap
tmp=kmeans(ex.expr.norm,center, iter.max=100)
cl = tmp$cluster
ex.expr.cl.means = round(do.call("rbind",tapply(1:length(cl), cl, function(x){
  colMeans(ex.expr.norm[x,])
})),digits=2)
ex.gene.cor = cor(t(ex.expr.norm), t(ex.expr.cl.means))
ex.gene.cl = apply(ex.gene.cor, 1, which.max)
ord = order(ex.gene.cl, order(rowMaxs(ex.gene.cor[names(ex.gene.cl),]),decreasing=T))
gene.col = jet.colors(length(unique(ex.gene.cl)))[ex.gene.cl]
pdf("ex.gene.cl.8.heatmap.pdf",height=12, width=5)
heatmap.2(ex.expr.norm[ord,], Rowv=NULL, Colv=NULL, cexRow=0.1,trace="none",RowSideColors=gene.col[ord],col=blue.red(100))
dev.off()

# testing enrichment of peak clusters in the gene clusters (again using Fisher test)
select.peak = names(ex.peak.cl)[all.peaks[as.integer(names(ex.peak.cl))]$gene.id %in% names(ex.gene.cl)]
tb= table(ex.peak.cl[select.peak], ex.gene.cl[all.peaks[as.integer(select.peak)]$gene.id])
tb.pval = matrix(1, nrow=nrow(tb), ncol=ncol(tb))
for(i in 1:nrow(tb)){
  for(j in 1:ncol(tb)){
    a = sum(tb[i,])
    b = sum(tb[,j])
    tmp.tb = matrix(c(sum(tb)- a -b + tb[i,j],  a- tb[i,j], b-tb[i,j], tb[i,j]), nrow=2, ncol=2)
    tmp=fisher.test(tmp.tb)
    tb.pval[i,j] = -log10(tmp$p.value)
    if(tmp$estimate < 1){
      tb.pval[i,j] = - tb.pval[i,j]
    }
  }
}

colnames(tb.pval) = colnames(tb)
row.names(tb.pval) =row.names(tb)

# building the heatmap for the peak-gene enrichment
tmp.df = as.data.frame(as.table(tb.pval))
colnames(tmp.df) = c("peak.cl", "gene.cl", "logP")
tmp = as.integer(tmp.df$gene.cl)
tmp.df$gene.cl = factor(tmp, levels=max(tmp):1)
tmp.df$logP[tmp.df$logP > 30] = 30
tmp.df$logP[tmp.df$logP < -20] = -20
write.table(tmp.df,"ex.peak.gene.cl.overlap.values.tsv",sep="\t",row.names=F,quote=F)
p <- ggplot(tmp.df, aes(peak.cl, gene.cl)) + geom_tile(aes(fill = logP),  colour = "white") + scale_fill_gradient2(low = "blue",mid="white",high = "red") 
pdf("ex.peak.gene.cl.overlap.pdf")
p
dev.off()

save(ex.peak.cl, ex.gene.cl, ex.expr.cl.means, ex.peak.cl.means, file="ex.gene.peak.cl.rda")

# Profiles of peak and gene clusters (again top 500 for peaks)
cl = ex.peak.cl
ex.peak.cor = cor(t(ex.peak.norm), t(ex.peak.cl.means))
select = row.names(ex.peak.norm)[rowMaxs(ex.peak.cor) > 0.8 & apply(ex.peak.cor, 1, which.max) == cl]
n.cl = length(unique(cl))
select.peak=sapply(1:n.cl, function(i){
  tmp.select= select[cl[select] == i]
  tmp.select=tmp.select[head(order(ex.peak.cor[tmp.select, i],decreasing=T), 500)]
},simplify=F)
select.peak = unique(unlist(select.peak))
tmp.df= as.data.frame(as.table(ex.peak.norm[select.peak,]))
tmp.df$cl = factor(cl[as.character(tmp.df$Var1)])
tmp.df$Var2 = factor(as.character(tmp.df$Var2), levels=cols)
g= ggplot(tmp.df, aes(Var2, Freq,group= Var1, alpha=0.1, color=cl)) +geom_line() + facet_grid(cl~.) + scale_color_manual(values=jet.colors(length(unique(ex.peak.cl))))
g = g+   theme( axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none",legend.background=element_blank(),panel.grid.major=element_blank(),panel.background=element_blank())
g1 = g + stat_summary(aes(y=Freq, group=Var2),fun.y="mean", colour="black", geom="point")
pdf("ex.peak.cl.pdf",height=10,width=4)
g1
dev.off()

cl = ex.gene.cl
tmp.df= as.data.frame(as.table(ex.expr.norm[names(cl),]))
tmp.df$cl = as.factor(cl[as.character(tmp.df$Var1)])
tmp.df$Var2 = factor(as.character(tmp.df$Var2), levels=cols)
g= ggplot(tmp.df, aes(Var2, Freq,group= Var1, alpha=0.1, color=cl)) +geom_line() + facet_grid(cl~.) + scale_color_manual(values=jet.colors(length(unique(ex.gene.cl))))
g = g+   theme( axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none",legend.background=element_blank(),panel.grid.major=element_blank(),panel.background=element_blank())
g1 = g + stat_summary(aes(y=Freq, group=Var2),fun.y="mean", colour="black", geom="point") 
pdf("ex.gene.cl.pdf",height=10,width=4)
g1
dev.off()


# computing AME enrichment between clusters 
gr = all.peaks[as.integer(names(ex.peak.cl))]
seq = getSequence(gr)
names(seq)= paste0(seqnames(gr),":",start(gr),"_",end(gr))
# filtering the peaks based on their similarity to each cluster (remove other clusters from the negative set that are too similar)
for(i in unique(ex.peak.cl)){
  select = as.vector(ex.peak.cl == i) & width(gr) < 400
  file=file.path("motif_seq", paste0("cl.",i,".fasta"))
  writeXStringSet(seq[select],file)
  pos = ex.peak.cl.means[i,]> 0.5
  neg=which(rowMax(ex.peak.cl.means[,pos,drop=F]) < 0.5)
  not.select = as.vector(ex.peak.cl %in% neg) & width(gr) < 400
  file=file.path("motif_seq", paste0("cl.ex",i,".fasta"))
  writeXStringSet(seq[not.select],file)
}
tmp = sapply(sort(unique(ex.peak.cl)), function(i){
  pos = ex.peak.cl.means[i,]> 0.5
  rowMax(ex.peak.cl.means[,pos,drop=F]) < 0.5
})
colnames(tmp) = row.names(tmp) = 1:8
write.table(tmp, file="ame.contrast.group.csv",sep=",")

# enrichment of motifs across clusters using all neural types (no contrast groups)
load("select.motif.rda")
load("gene.peak.cl.rda")
select = as.character(select.gr$peak.id) %in% names(peak.cl)  & width(select.gr)<400
select.gr  = select.gr[select]
select.gr$ex.pPadj = all.peaks$ex.pPadj[select.gr$peak.id.1]
levels(peak.cl)=1:length(levels(peak.cl))
tb=table(select.gr$tf, peak.cl[as.character(select.gr$peak.id)])
tb.pval = matrix(1, nrow=nrow(tb), ncol=ncol(tb))
for(i in 1:nrow(tb)){
  for(j in 1:ncol(tb)){
    a = sum(tb[i,])
    b = sum(tb[,j])
    tmp.tb = matrix(c(sum(tb)- a -b + tb[i,j],  a- tb[i,j], b-tb[i,j], tb[i,j]), nrow=2, ncol=2)
    tmp=fisher.test(tmp.tb)
    tb.pval[i,j] = -log10(tmp$p.value)
    if(tmp$estimate < 1){
      tb.pval[i,j] = - tb.pval[i,j]
    }
  }
}
row.names(tb.pval)=row.names(tb)
colnames(tb.pval)=colnames(tb)
tmp.df = as.data.frame(as.table(tb.pval))
colnames(tmp.df) = c("motif", "peak.cl", "logP")
p <- ggplot(tmp.df, aes(peak.cl, motif)) + geom_tile(aes(fill = logP),  colour = "white") + scale_fill_gradient2(low = "blue",mid="white",high = "red") 
pdf("motif.peak.cl.pdf")
p
dev.off()

# enrichment of motifs across excitatory clusters not using contrast groups
load("select.motif.rda")
load("ex.gene.peak.cl.rda")
select = as.character(select.gr$peak.id) %in% names(ex.peak.cl)  & width(select.gr)<400
select.gr  = select.gr[select]
#tmp = table(select.gr$tf, select.gr$peak.id)
#tb=do.call("cbind",tapply(colnames(tmp), ex.peak.cl[colnames(tmp)], function(x) rowSums(tmp[,x] > 0)))
tb=table(select.gr$tf, ex.peak.cl[as.character(select.gr$peak.id)])
tb.pval = matrix(1, nrow=nrow(tb), ncol=ncol(tb))
for(i in 1:nrow(tb)){
  for(j in 1:ncol(tb)){
    a = sum(tb[i,])
    b = sum(tb[,j])
    tmp.tb = matrix(c(sum(tb)- a -b + tb[i,j],  a- tb[i,j], b-tb[i,j], tb[i,j]), nrow=2, ncol=2)
    tmp=fisher.test(tmp.tb)
    tb.pval[i,j] = -log10(tmp$p.value)
    if(tmp$estimate < 1){
      tb.pval[i,j] = - tb.pval[i,j]
    }
  }
}
row.names(tb.pval)=row.names(tb)
colnames(tb.pval)=colnames(tb)
tmp.df = as.data.frame(as.table(tb.pval))
colnames(tmp.df) = c("motif", "peak.cl", "logP")
p <- ggplot(tmp.df, aes(peak.cl, motif)) + geom_tile(aes(fill = logP),  colour = "white") + scale_fill_gradient2(low = "blue",mid="white",high = "red") 
pdf("ex.motif.peak.cl.pdf")
p
dev.off()

# parse AME and get enrichment using AME results for comparisons between contrast groups (instead of manual Fisher test -- AME may use Fisher test)
parse_ame <- function(fn)
  {
    lines=readLines(fn)[-(1:12)]
    tmp = do.call("rbind",strsplit(lines, "[ \\(\\)]"))
    tb= data.frame(tmp[,c(8,9,13,17)],stringsAsFactors=F)
    colnames(tb)=c("ID","TF","pval","padj")
    tb$pval = as.numeric(tb$pval)
    tb$padj = as.numeric(tb$padj)
    tb
  }

ame_result=do.call("rbind",sapply(sort(unique(ex.peak.cl)), function(i){
  tb1 = parse_ame(file.path("ame_out", i, "ame.txt"))
  tb2 = parse_ame(file.path("ame_out", paste0("ex",i), "ame.txt"))
  tb1$cat = "enriched"
  tb2$cat = "depleted"
  tb=rbind(tb1, tb2)
  tb$peak.cl = i
  tb
},simplify=F))
ame_result$logP = -log10(ame_result$padj)
select=ame_result$cat == "depleted"
ame_result$logP[select] =  - ame_result$logP[select]
select.motif = list(DLX="MA0882.1", NEUROD="MA0000.1", TBR1=c("MA0802.1","MA0800.1"),FOXP=c("MA0593.1","MA0481.1"), NFIA="MA0670.1", RORB=c("MA0071.1", "MA0072.1"), POU3F=c("MA0787.1","MA0789.1", "MA0786.1", "MA0788.1"), MEF2=c("MA0487.1","MA0773.1","MA0660.1","MA0052.3"),FOS=c("MA0477.1","MA0478.1","MA0476.1","MA0099.2"), EGR1=c("MA0472.2","MA0732.1","MA0162.2","MA0733.1"), MEIS = c("MA0775.1", "MA0774.1"), RFX3=c("MA0798.1","MA0600.2","MA0799.1","MA0365.1","MA510.2"),CUX=c("MA0755.1","MA0679.1","MA0757.1","MA0756.1"))

n.cl = sort(unique(ex.peak.cl))
ex.ame.score=sapply(select.motif, function(x){
  select=ame_result$ID %in% x
  tmp = with(ame_result[select,], tapply(logP, peak.cl, function(x)x[which.max(abs(x))]))
  score= setNames(rep(0, length(n.cl)), n.cl)
  score[names(tmp)] = tmp
  score
})
save(ex.ame.score,file="ex.ame.score.rda")

tmp.df = as.data.frame(as.table(ex.ame.score))
colnames(tmp.df) = c("peak.cl","motif","logP")
p <- ggplot(tmp.df, aes(peak.cl, motif)) + geom_tile(aes(fill = logP),  colour = "white") + scale_fill_gradient2(low = "blue",mid="white",high = "red") 
pdf("ex.motif.peak.cl.ame.pdf")
p
dev.off()


#####
# Checking results
#####
all.peaks$ex.peak.cl=NA
all.peaks$ex.peak.cl[as.integer(names(ex.peak.cl))] = ex.peak.cl
all.peaks$ex.gene.cl=NA
tmp=match(all.peaks$gene.id, names(ex.gene.cl))
select = !is.na(tmp)
all.peaks$ex.gene.cl[select] = ex.gene.cl[tmp[select]]

tmp=load("gene.peak.cl.rda")
all.peaks$peak.cl=NA
all.peaks$peak.cl[as.integer(names(peak.cl))] = peak.cl
all.peaks$gene.cl=NA
tmp=match(all.peaks$gene.id, names(gene.cl))
select = !is.na(tmp)
all.peaks$gene.cl[select] = gene.cl[tmp[select]]

logP = -log10(select.gr$pvalue) - 4
tmp=tapply(logP, list(select.gr$peak.id, select.gr$merged.tf), sum, na.rm=T)
tmp[is.na(tmp)] = 0
peak.motif.score = matrix(0, nrow=length(all.peaks), ncol= ncol(tmp))
peak.motif.score[as.integer(row.names(tmp)),] = tmp
colnames(peak.motif.score) = colnames(tmp)
values(all.peaks) = cbind(values(all.peaks), as.data.frame(peak.motif.score))
save(all.peaks, file="all.peaks.rda")

for(x in colnames(peak.motif.score)){
  values(all.peaks)[,x] = peak.motif.score[,x]
}


select = with(values(all.peaks), which(ex.gene.cl == 4 & ex.peak.cl==4& RFX > 0 & gPadj < 0.001))
tmp = all.peaks[select,]
tmp = tmp[order(tmp$ex.gPadj),]
head(tmp[,c(1:15,64, 51:63)])


select = with(values(all.peaks), which(ex.gene.cl == 4 & ex.peak.cl==4& FOXP > 0 & gPadj < 0.001))
tmp = all.peaks[select,]
tmp = tmp[order(tmp$ex.gPadj),]
head(tmp[,c(1:15,64, 51:63,65)])

select.gr = select.gr[select.gr$score > 10]

###Stard8
subsetByOverlaps(select.gr, all.peaks[218354])[,c(1:6,10)]
###Sty17
subsetByOverlaps(select.gr, all.peaks[189903])[,c(1:6,10)]
subsetByOverlaps(select.gr, all.peaks[189906])[,c(1:6,10)]
##Enpp2
subsetByOverlaps(select.gr, all.peaks[75415])[,c(1:6,10)]
subsetByOverlaps(select.gr, all.peaks[75421])[,c(1:6,10)]

select = with(values(all.peaks), which(ex.gene.cl == 1 & ex.peak.cl==1& NFIA > 0 & gPadj < 0.001))
tmp = all.peaks[select,]
tmp = tmp[order(tmp$ex.gPadj),]
head(tmp[,c(1:15,64, 51:63,65)],10)

##Foxp2
subsetByOverlaps(select.gr, all.peaks[168954])[,c(1:6,9:10)]
##NFIA

#Prkcb
subsetByOverlaps(select.gr, all.peaks[136169])[,c(1:6,9:10)]
subsetByOverlaps(select.gr, all.peaks[190308])[,c(1:6,9:10)]
subsetByOverlaps(select.gr, all.peaks[190332])[,c(1:6,9:10)]


select = with(values(all.peaks), which(ex.gene.cl %in% c(3,6) & ex.peak.cl==3& gPadj < 0.001))
tmp = all.peaks[select,]
tmp = tmp[order(tmp$ex.gPadj),]
head(tmp[,c(1:15,64, 51:63,65)],10)


tmp = peak.motif.score > 0.5
tb = t(tmp) %*% tmp
motif.occur = colSums(tmp)
total = nrow(peak.motif.score)
tb.pval = matrix(1, nrow=nrow(tb), ncol=ncol(tb))
for(i in 1:nrow(tb)){
  for(j in 1:ncol(tb)){
    cat(i,j, "\n")
    a = motif.occur[i]
    b = motif.occur[j]
    tmp.tb = matrix(c(total- a -b + tb[i,j],  a- tb[i,j], b-tb[i,j], tb[i,j]), nrow=2, ncol=2)
    tmp=fisher.test(tmp.tb)
    tb.pval[i,j] = -log10(tmp$p.value)
    if(tmp$estimate < 1){
      tb.pval[i,j] = - tb.pval[i,j]
    }
  }
}
colnames(tb.pval)=row.names(tb.pval)= colnames(peak.motif.score)



####find gene-gene regulatory network
# find genes in ex clusters
tmp=load("ex.gene.peak.cl.rda")
select.peak = names(ex.peak.cl)[all.peaks[as.integer(names(ex.peak.cl))]$gene.id %in% names(ex.gene.cl)]
tb= table(ex.peak.cl[select.peak], ex.gene.cl[all.peaks[as.integer(select.peak)]$gene.id])
tb.pval = matrix(1, nrow=nrow(tb), ncol=ncol(tb))
for(i in 1:nrow(tb)){
  for(j in 1:ncol(tb)){
    a = sum(tb[i,])
    b = sum(tb[,j])
    tmp.tb = matrix(c(sum(tb)- a -b + tb[i,j],  a- tb[i,j], b-tb[i,j], tb[i,j]), nrow=2, ncol=2)
    tmp=fisher.test(tmp.tb)
    tb.pval[i,j] = -log10(tmp$p.value)
    if(tmp$estimate < 1){
      tb.pval[i,j] = - tb.pval[i,j]
    }
  }
}

colnames(tb.pval) = colnames(tb)
row.names(tb.pval) =row.names(tb)

# for each gene set, make sure the peaks are in clusters that match the gene clusters
tmp.df = as.data.frame(as.table(tb.pval))
colnames(tmp.df) = c("peak.cl", "gene.cl", "logP")
select.peak.gene.df = tmp.df[tmp.df$logP >5, ]
select.peak.gene.df$peak.gene.id = paste(select.peak.gene.df$peak.cl, select.peak.gene.df$gene.cl)

# make sure motif enrichment is in the correct group
load("ex.ame.score.rda")
tmp.df = as.data.frame(as.table(ex.ame.score))
colnames(tmp.df) = c("peak.cl","motif","logP")
select.peak.motif.df = tmp.df[tmp.df$logP>5,]
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





