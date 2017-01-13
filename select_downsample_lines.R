library(dplyr)
options(stringsAsFactors=F)

args <- commandArgs(TRUE)

frags <- as.numeric(args[2]) * 1e6

print(args)

bamstats <- read.delim(args[1],sep=" ",header=F)

start_pairs <- as.numeric(bamstats[17,1])

rows_to_sample <- data.frame(id1=seq(1,start_pairs*2,2),
                             id2=seq(2,start_pairs*2+1,2)) %>%
  sample_n(as.numeric(frags)) %>%
  unlist() %>%
  as.data.frame()

names(rows_to_sample) <- "id"

rows_to_sample <- rows_to_sample %>% arrange(id)

write.table(rows_to_sample,args[3],quote=F,col.names=F,row.names=F)
