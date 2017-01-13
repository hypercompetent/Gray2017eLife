# Retrieve TSS regions from the UCSC refGene table
get_tss_regions <- function(symbols=NULL,expand=5000,genome="mm10") {
  library(dplyr)
  library(rtracklayer)
  
  session <- browserSession("UCSC")
  genome(session) <- genome
  refgene <- getTable(ucscTableQuery(session,table="refGene"))
  
  if(length(expand) == 1) {
    expand5 <- expand
    expand3 <- expand
  } else if(length(expand) == 2) {
    expand5 <- expand[1]
    expand3 <- expand[2]
  }
  
  if(!is.null(symbols)) {
    
    refgene <- refgene %>%
      filter(name2 %in% symbols)
    
  }
  
  tss_windows <- refgene %>%
    mutate(tss = ifelse(strand == "+", txStart, txEnd)) %>%
    mutate(tss_up = ifelse(strand == "+", tss - expand5, tss - expand3),
           tss_dn = ifelse(strand == "+", tss + expand3, tss + expand5)) %>%
    select(chrom,tss_up,tss_dn,name2,score,strand)
  
  names(tss_windows) <- c("chr","start","end","name","score","strand")
  
  return(tss_windows)
  
}

# Convert data.frames in BED-like format to GRanges objects
bed_to_GRanges <- function(bed) {
  library(rtracklayer)
  
  gr <- GRanges(seqnames=bed$chr,
                IRanges(start=bed$start,
                        end=bed$end),
                strand=bed$strand,
                mcols=bed[,c("name","score")])
  
  return(gr)
}

# Get the distance to the nearest GRanges object.

getNearestDist <- function(gr1, gr2, return.nearest=F, r1.pos="mid", r2.pos="mid", strand=2)
{
  r1 <- ranges(gr1)
  r2 <- ranges(gr2)
  map.id <- rep(NA, length(r1))
  dist <- c()
  for(chr in sort(unique(seqnames(gr1)))){
    select1 <- seqnames(gr1)==chr
    gr1.select <- gr1[select1,]
    select2 <- seqnames(gr2)==chr
    gr2.select <- gr2[select2,]
    l <- nearest(ranges(gr1.select), ranges(gr2.select))
    map.id[which(select1)] <- which(select2)[l]
    
    start1=start(gr1.select)
    end1 = end(gr1.select)
    end5.1 <- start1
    end3.1 <- end1
    neg1 <- strand(gr1.select) == "-"
    end5.1[which(neg1)] <- end1[which(neg1)]      
    end3.1[which(neg1)] <- start1[which(neg1)]      
    mid1 <- as.integer( (start1+end1)/2)
    
    start2=start(gr2.select)
    end2 = end(gr2.select)
    end5.2 <- start2
    end3.2 <- end2
    neg2 <- strand(gr2.select) == "-"
    end5.2[which(neg2)] <- end2[which(neg2)]      
    end3.2[which(neg2)] <- start2[which(neg2)]      
    mid2 <- as.integer( (start2+end2)/2)
    
    pos1 <- switch(r1.pos,
                   mid=mid1, 
                   start=start1,
                   end = end1, 
                   end5 = end5.1,
                   end3 = end3.1)
    pos2 <- switch(r2.pos,
                   mid=mid2, 
                   start=start2,
                   end = end2, 
                   end5 = end5.2,
                   end3 = end3.2)             
    tmp.dist <- as.integer(pos1 - pos2[l])
    if(strand==2){
      tmp.dist[which(neg2[l])] <- - tmp.dist[which(neg2[l])]
    }      
    if(strand==1){
      tmp.dist[which(neg1)] <- - tmp.dist[which(neg1)]
    }
    dist[which(select1)] <- tmp.dist
  }
  if(!return.nearest){
    dist
  }
  else{
    gr <- gr2[map.id,]
    values(gr)$dist <- dist
    gr$map.id = map.id
    gr
  }
}


# Compute differentially expressed genes
DESeq.genes.pw <- function(dat,cl, dds.file="dds.rda",mc.cores=4){
  require("DESeq2")
  require("parallel")
  if(is.null(dds.file) || !file.exists(dds.file)){
    df= data.frame(cl=cl)
    dds = DESeqDataSetFromMatrix(dat, colData=df, design=~ cl)
    dds=DESeq(dds)
    if(!is.null(dds.file)){
      save(dds, file=dds.file)
    }
  }
  tmp.cl = sort(unique(cl))
  pairs = data.frame(X=as.character(rep(tmp.cl,length(tmp.cl))), Y=as.character(rep(tmp.cl, rep(length(tmp.cl), length(tmp.cl)))), stringsAsFactors=F)
  pairs = pairs[pairs[,1]<pairs[,2],]
  de.df=sapply(1:nrow(pairs), function(i){
    print(pairs[i,])
    x=pairs[i,1]
    y=pairs[i,2]
    res=results(dds, contrast = c("cl", x,y))
    colnames(res)[2]="lfc"
    res
  },simplify=F)
  names(de.df) = paste(pairs[,1],pairs[,2],sep="_")
  return(de.df)
}

# Retrieve sequences for a GRanges object from BioStrings:
getSequence <- function(gr, genome=Mmusculus,use.strand=T)
{
  all.seq <- DNAStringSet(rep("", length(gr)))
  for(chr in (sort(unique(seqnames(gr))))){
    select <- as.vector(seqnames(gr)==chr)
    gr.select <- gr[select]
    chr.seq <- genome[[chr]]                            
    seq <- DNAStringSet(Views(chr.seq, start=start(gr.select), end=end(gr.select)))
    if(use.strand){
      neg <- as.vector(strand(gr.select) == "-" )
      seq[neg] <- reverseComplement(DNAStringSet(seq[neg]))
    }
    all.seq[select] <- seq      
  }
  all.seq
}


values_to_colors <- function (x, minval = 0, maxval = NULL, colorset = c("darkblue","dodgerblue", "gray80", "orangered", "red")) {
  heat_colors <- colorRampPalette(colorset)(1001)
  if (is.null(maxval)) {
    maxval <- max(x)
  }
  heat_positions <- unlist(round((x - minval)/(maxval - minval) * 
                                   1000 + 1, 0))
  colors <- heat_colors[heat_positions]
  colors
}