#!/usr/bin/env Rscript

options(stringsAsFactors = F)
library(dplyr)
library(patchwork)
library(ggplot2)

args <- as.vector(commandArgs(T))

## read data ####
stopifnot(file.exists('read_metrics.csv'))
meta <- read.csv('read_metrics.csv')

dat <- list()

if (file.exists('fragment_length.rds')){
  dat$frag_lens <- readRDS('fragment_length.rds')
}

if (file.exists('original_peak_metrics.csv')){
  dat$npeaks <- read.csv('original_peak_metrics.csv')
}

if (file.exists('original_peak_widths.rds')){
  dat$wpeaks <- readRDS('original_peak_widths.rds')
}

if (file.exists('replicated_peak_metrics.csv')){
  dat$nreps <- read.csv('replicated_peak_metrics.csv')
}

if (file.exists('consensus_peak_metrics.csv')){
  dat$nconps <- read.csv('consensus_peak_metrics.csv')
}

if (file.exists('consensus_peaks.rds')){
  dat$conps <- readRDS('consensus_peaks.rds')
}

## initiate variables ####
funcs <- list()
figs <- list()
targets <- sort(setdiff(unique(meta$target), 'IgG')) # targets to make QC plots

## Number of original peaks ####
qc <- 'count_peaks' 
funcs[[qc]] <- function(df, tgt, var.x = 'npeak', var.y = 'sample_group', var.color = 'toIgG', add_facet = facet_grid(~caller, scales = 'free'), var.shape = 'filtered', xlab = 'Total Peaks (in thousand)', ylab = '', scale_color = scale_color_manual('', values = c("IgG_controlled"="indianred", "Target_only"="steelblue")), shape = '', plot.title = 'Peaks from each sample'){
  require(ggplot2)
  
  df <- as.data.frame(df) %>% 
    mutate(sample_replicate = ifelse(is.na(sample_replicate), 1, sample_replicate))
  df %>% 
    mutate(x = df[,var.x]/1e3, y = df[,var.y], color = df[,var.color], shape = df[,var.shape]) %>% 
    subset(target %in% tgt) %>%
    ggplot(mapping = aes(x = x, y = y, color = color, shape = shape))+
    geom_point()+
    add_facet+
    scale_color+
    scale_shape_manual(values = c(filtered=19, unfiltered=1))+
    labs(
      x = xlab, y = ylab, 
      shape = shape, title = plot.title)+
    theme_bw()+
    theme(
      text = element_text(size = 8),
      axis.text.x = element_text(angle = 90, hjust = 1),
      legend.position = 'top'
    )
}

df <- dat$npeaks %>% 
  mutate(id = paste(sample_group, sample_replicate, sep = '_'))

figs[[qc]] <- lapply(targets, funcs[[qc]], df = df, var.y = 'id')
names(figs[[qc]]) <- targets


## Width of filtered original peaks ####
qc <- 'peak_width'

funcs[[qc]] <- function(df, tgt, var.x = 'length', var.y = 'id', var.color = 'caller', xlab = 'Peak Width (bp)', ylab = '', scale_color = scale_color_manual('', values = c(SEACR='darkblue', MACS2narrow='red', MACS2broad='brown')), plot.title = 'Peak width'){
  require(ggplot2)
  
  ## determine data range
  length.max <- as.data.frame(df) %>% 
    subset(target %in% tgt) %>% 
    group_by(sample_group, caller) %>% 
    reframe(
      qt = quantile(length, 0.75)
    ) %>% 
    reframe(
      max = max(qt)
    ) %>% 
    as.matrix() %>% 
    as.vector()
  
  ## make plots
  as.data.frame(df) %>% 
    mutate(x = .[,var.x], y = .[,var.y], color = .[,var.color]) %>% 
    subset(target %in% tgt) %>%
    ggplot(mapping = aes(x = x, y = y, color = color, weight = weight))+
    geom_violin(bw=5, trim = T, draw_quantiles = 0.5, orientation = 'y')+
    scale_color+
    coord_cartesian(xlim = c(0, length.max))+
    labs(x=xlab, y=ylab, title = plot.title)+
    theme_bw()+
    theme(
      text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = 'top'
    )
}

df <- dat$wpeaks %>% 
  mutate(id = paste(sample_group, sample_replicate, sep = '_'))


figs[[qc]] <- lapply(targets, funcs[[qc]], df = df)
names(figs[[qc]]) <- targets


## Fraction of reads in filtered peaks ####
qc <- 'peak_frip'

funcs[[qc]] <- function(df, tgt, var.x = 'frip', var.y = 'id', var.color = 'caller', add_facet = NULL, xlab = 'Fraction of reads in peak', ylab = '', scale_color = scale_color_manual('', values = c(SEACR='darkblue', MACS2narrow='red', MACS2broad='brown')), plot.title = 'Fraction of reads in peak'){
  require(ggplot2)
  
  as.data.frame(df) %>% 
    mutate(x = .[,var.x], y = .[,var.y], color = .[,var.color]) %>% 
    subset(target %in% tgt) %>%
    ggplot(mapping = aes(x = x, y = y, color = color))+
    geom_point()+
    add_facet+
    scale_color+
    labs(
      x = xlab, y = ylab, 
      title = plot.title
    )+
    theme_bw()+
    theme(
      text = element_text(size = 8),
      axis.text.x = element_text(angle = 90, hjust = 1),
      plot.title = element_text(size = 10),
      legend.position = 'top'
    )
  
  
}

df <- dat$npeaks %>% 
  subset(filtered %in% 'filtered') %>% 
  mutate(id = paste(sample_group, sample_replicate, sep = '_'))

figs[[qc]] <- lapply(targets, funcs[[qc]], df = df)
names(figs[[qc]]) <- targets



## compare consensus peaks across sample groups ####
qc <- 'rep2conp'

funcs[[qc]] <- function(df, plot.title){
  df <- df %>% 
    mutate(
      shared = ifelse(grepl(',', sample.groups), 'shared', 'unique')
    ) %>% 
    rowwise() %>% 
    reframe(
      peak.id = peak.id,
      shared = shared,
      sample.group = unlist(strsplit(sample.groups, split = ','))
    )
  
  df %>% 
    ggplot(aes(y = sample.group, fill = factor(shared, levels = c('unique', 'shared'))))+
    geom_bar()+
    labs(y = '', x = 'peak count', fill = '', title = plot.title)+
    theme_bw(base_size = 8)
  
}

figs[[qc]] <- lapply(
  paste(rep(targets, each = 3), c('MACS2narrow', 'MACS2broad', 'SEACR'), sep = ':'),
  function(i){
    if (i %in% names(dat$conps)){
      funcs[[qc]](dat$conps[[i]], i)
    }else{
      plot_spacer()
    }
  }
)

## Compare MACS2 and SEACR by consensus peaks ####
qc <- 'seacr2macs'

funcs[[qc]] <- function(peak.list){
  require(GenomicRanges)
  
  tgt <- gsub(':.*', '', names(peak.list)[1])
  
  if(sum(sapply(peak.list, nrow)>0) < 2){return(NULL)}
  
  peak.list <- peak.list[sapply(peak.list, nrow)>0]
  peak.list <- lapply(
    peak.list, function(x){
      if(nrow(x)>0){
        makeGRangesFromDataFrame(x, ignore.strand = T, seqnames.field = 'seqnames', start.field = 'start', end.field = 'end', starts.in.df.are.0based = T)
      }else{
        NULL
      }
    }
  ); names(peak.list) <- gsub('.*:', '', names(peak.list))
  
  conp <- GRanges()
  for (i in 1:length(peak.list)){
    conp <- c(conp, peak.list[[i]])
  }
  conp <- reduce(conp)
  
  overlap2conp <- as.data.frame(sapply(
    peak.list, function(pk){countOverlaps(conp, pk)>0}
  ))
  
  venn::venn(
    overlap2conp,
    ggplot = T,
    box = F
  )+
    labs(title = tgt)+
    theme(
      text = element_text(size = 10),
      plot.margin = unit(rep(0.5, 4), 'line'),
      plot.title = element_text(size = 12)
    )
  
}

figs[[qc]] <- lapply(
  split(dat$conps, gsub(':.*', '', names(dat$conps))),
  funcs[[qc]]
)

## Save results ####
saveRDS(funcs, 'funcs.rds')
saveRDS(figs, 'figs.rds')
saveRDS(dat, 'data.rds')
