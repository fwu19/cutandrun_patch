#!/usr/bin/env Rscript

options(stringsAsFactors = F)
options(scipen = 99)
library(dplyr)
library(GenomicRanges)

## sample sheet ####
meta <- read.csv('read_metrics.csv') 

targets <- sort(setdiff(unique(meta$target), 'IgG')) # targets to make QC plots

peak.list <- grep('narrowPeak$|broadPeak$|bed$', list.files('./'), value = T)
rip.list <- grep('readsInPeak.csv$', list.files('./'), value = T)


## read original peaks ####
peaks <- lapply(
  peak.list,
  function(fname){
    if(file.size(fname) > 0){
      GenomicRanges::makeGRangesFromDataFrame(
        read.delim(fname, header = F)[1:3], seqnames.field = 'V1', start.field = 'V2', end.field = 'V3', starts.in.df.are.0based = T
      )
    }
  }); names(peaks) <- peak.list
saveRDS(peaks, 'original_peaks.rds')


## summarize original peaks ####
npeaks <- bind_rows(
  mapply(function(pk, fname){
    data.frame(file = fname, npeak = length(pk))
  }, peaks, names(peaks), SIMPLIFY = F
  )
) %>% 
  mutate(sample_id = gsub('.toIgG.*|.noControl.*|.IgG.norm.*|.top01.*','',basename(file))) %>% 
  left_join(
    meta %>% select(sample_id, sample_group, target, sample_replicate, bt2_total_aligned_target), 
    by = c("sample_id" = "sample_id")
  ) %>% 
  mutate(
    caller = case_when(
      grepl('stringent', basename(file)) ~ 'SEACR', 
      grepl('broad', basename(file)) ~ 'MACS2broad',
      grepl('narrow', basename(file)) ~ 'MACS2narrow',
      TRUE ~ 'others'
    ),
    toIgG = ifelse(grepl('toIgG|\\.IgG.norm', file),"IgG_controlled", "Target_only"),
    filtered = ifelse(grepl('filtered', file),'filtered','unfiltered')
  ) %>% 
  left_join(
    bind_rows(lapply(
      rip.list, 
      function(fname){
        if(file.size(fname) > 0){
          read.csv(fname, header = F, col.names = c('file', 'reads_in_peak'))
        }else{
          return(NULL)
        }
      }
    )) ,
    by = c('file' = 'file')
  ) %>% 
  mutate(frip = reads_in_peak/bt2_total_aligned_target)
write.table(npeaks, 'original_peak_metrics.csv', sep = ',', quote = F, row.names = F)

## compute peak widths ####
wpeaks <- bind_rows(mapply(
  function(pk, file){
    if(length(pk)>0){
      tab <- table(width(pk))
      data.frame(
        file = file,
        length = as.integer(names(tab)),
        count = as.vector(tab)
      ) %>% 
        mutate(weight = count/ sum(count))
      
    }
  }, peaks, names(peaks), SIMPLIFY = F)
) %>% 
  left_join(
    npeaks %>% select(file, sample_id, sample_group, target, sample_replicate, caller),
    by = c('file' = 'file')
  )

## save results ####
saveRDS(wpeaks, 'original_peak_widths.rds')

