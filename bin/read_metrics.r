#!/usr/bin/env Rscript

options(stringsAsFactors = F)
library(dplyr)

args <- commandArgs(T)
if (length(args) < 1){
  stop("Please provide the following positional arguments, <path/to/sample_sheet.csv> <path/to/meta_table_ctrl.csv> <path/to/fragment_length>")
}

## sample sheet ####
ss <- read.csv(args[1])
if ('id' %in% colnames(ss)){ss <- dplyr::select(ss, !id)}

## metadata ####
meta <- read.csv(args[2])

## run the following only once
if (! 'target_alignment_rate' %in% colnames(meta)){ 
  meta <- meta %>% 
  mutate(
    target_alignment_rate = bt2_total_aligned_target/bt2_total_reads_target
    )
  }
meta <- meta %>% 
  left_join(
    ss %>% dplyr::select(group, replicate, peak_type, sample_id, control_id, sample_group, sample_replicate, target), 
    by = c('group' = 'group', 'replicate' = 'replicate'))

write.table(meta, 'read_metrics.csv', sep = ',', quote = F, row.names = F)


## fragment lengths ####
flist <- list.files('./', pattern = 'fragmentLen', full.names = T, recursive = T)
bind_rows(lapply(
  flist,
  function(fname){
    if(file.size(fname) > 0 ){
      read.delim(fname, header = F, col.names = c('length','count'), colClasses = 'numeric') %>% 
        mutate(
          weight = count/sum(count), 
          sample_id = gsub('.fragmentLen.txt','',basename(fname)))
    }
  })) %>% 
  left_join(
    meta %>% select(sample_id, sample_group, sample_replicate, target),
    by = 'sample_id'
  ) %>% 
  saveRDS('fragment_length.rds')
