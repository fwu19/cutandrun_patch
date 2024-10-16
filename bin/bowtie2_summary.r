#!/usr/bin/env Rscript

# This script combines metrics from Cut&Run analysis

options(stringsAsFactors = F)
library(dplyr)

args <- as.vector(commandArgs(T))

## functions
bt2_metrics <-  function(bt2log){
  if(file.exists(bt2log)){
    if (file.size(bt2log) > 0){
    bt2 <- as.data.frame(matrix(gsub(' .*', '', trimws(scan(bt2log, what = 'character', sep = "\n"))), nrow = 1))
    colnames(bt2) <- c('input.reads','input.pairs','unmapped','uniquely.mapped','multimapper','overall.alignment.rate')
    return(cbind(data.frame(file = basename(bt2log)), bt2))
    }else{
      return(NULL)
    }
  }
}
 
## call functions
x <- do.call(rbind,lapply(args, bt2_metrics)) %>% 
  mutate(
    bt2_total_aligned_yeast = as.numeric(uniquely.mapped) + as.numeric(multimapper),
    yeast_alignment_rate = bt2_total_aligned_yeast / as.numeric(input.pairs)
    ) %>% 
  dplyr::select(-c(input.reads, overall.alignment.rate))

write.table(x, 'qc/read_metrics.yeast.csv',sep = ',',quote = F,row.names = F)
