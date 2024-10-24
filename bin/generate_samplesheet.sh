#!/usr/bin/env bash

module load fhR/4.1.2-foss-2021b; export R_LIBS=/fh/fast/_SR/Genomics/user/fwu/R/x86_64-pc-linux-gnu-library/4.1 

echo "id,fastq_1,fastq_2" >fq.csv
# for fqdir in $@; do
# 	paste -d "," <(ls $fqdir/*_R1_* 2>/dev/null | sed "s/.*\///g" | sed "s/_S[0-9]\+_.*//" ) <(ls $fqdir/*_R1_* 2>/dev/null ) <(ls $fqdir/*_R2_* 2>/dev/null ) 
# done | egrep "^PC.*K562" >>fq.csv

for fqdir in $@; do
	paste -d "," <(find $fqdir/ | egrep "Unaligned/Project" | egrep "/PC.*K562.*_R1_" 2>/dev/null | sort | sed "s/.*\///g" | sed "s/_S[0-9]\+_.*//" ) <(find $fqdir/ | egrep "Unaligned/Project" | egrep "/PC.*K562.*_R1_" 2>/dev/null | sort ) <(find $fqdir/ | egrep "Unaligned/Project" | egrep "/PC.*K562.*_R2_" 2>/dev/null | sort ) 
done >>fq.csv

R --vanilla "--args fq.csv $(date +%y%m%d)" <<code
options(stringsAsFactor=F)
options(scipen=99)
library(dplyr)
args <- commandArgs(T)
ss <- read.csv(args[1]) %>% 
  mutate(
    flowcell=basename(gsub('.Unaligned.*', '', fastq_1)),
    target=case_when(
      grepl('K27me3|K27M3', id) ~ 'H3K27me3',
      grepl('Pol2Ser5', id) ~ 'Pol2S5',
      grepl('CTCF', id) ~ 'CTCF',
      grepl('Myc', id) ~ 'cMyc',
      TRUE ~ 'TBD'
    ),
    peak_type = case_when(
      target %in% 'H3K27me3' ~ 'broad',
      target %in% c('cMyc', 'Pol2S5') ~ 'narrow',
      TRUE ~ 'both'
    )
    
  ) %>% 
  group_by(flowcell, target) %>% 
  mutate(
    n = n(),
    i = 1:n,
    sample_replicate=sapply(strsplit(flowcell, split='_'), function(v){paste0(v[1],v[length(v)])})
  ) %>% 
  mutate(
    sample_replicate=ifelse(n > 1, paste0(sample_replicate, letters[i]), sample_replicate),
    group=paste(target, sample_replicate, sep = '_'),
    replicate = 1,
    control = '',
    control_replicate = '',
    sample_group='PE50',
    Sequencing.Type='NextSeq P2',
    Read.Length='50 / 50',
    analysis_dir = args[2]
    ) %>% 
  mutate(
    is_control = ifelse(control == '', 1, 0),
    sample_id = paste(group, replicate, sep = '_R'),
    control_id = ifelse(control == '', '', paste(control, control_replicate, sep = '_R'))
  ) %>% 
  relocate(
    group, replicate, fastq_1, fastq_2, control, control_replicate, peak_type, sample_group, target, sample_replicate, flowcell, Sequencing.Type, Read.Length, i, n, analysis_dir, is_control, sample_id, control_id
  )
ss %>% write.table(paste('samples', args[2], 'csv', sep = '.'), sep = ',', quote = F, row.names = F)



code