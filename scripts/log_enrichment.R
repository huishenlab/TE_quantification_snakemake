#!/usr/bin/env Rscripts

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))

options(scipen=999)

log_enrichment <- function(tes_cpm, sample){
  
  cells <- colnames(tes_cpm)[1:(ncol(tes_cpm) - 4)]
  for (i in 1:(ncol(tes_cpm) - 4)){
    cell <- colnames(tes_cpm)[i]
    te_cpm <- tes_cpm[,c(cell, "Geneid", "repName", "repClass", "repFamily")]
    te_cpm_grt_1_cpm <- te_cpm[te_cpm[,cell] > 1, ]
    
    te_cpm_grouped <- te_cpm %>% group_by(repFamily,repClass)
    te_cpm_grouped_family_counts <- te_cpm_grouped %>% dplyr::count(repFamily) 
    te_cpm_grouped_family_counts <- as.data.frame(te_cpm_grouped_family_counts)
    
    te_cpm_grt_1_grouped <- te_cpm_grt_1_cpm %>% group_by(repFamily,repClass)
    te_cpm_grt_1_grouped_family_counts <- te_cpm_grt_1_grouped %>% dplyr::count(repFamily) 
    te_cpm_grt_1_grouped_family_counts <- as.data.frame(te_cpm_grt_1_grouped_family_counts)
    
    te_cpm_lst_1_grouped_family_counts <- te_cpm_grouped_family_counts[!te_cpm_grouped_family_counts$repFamily %in% te_cpm_grt_1_grouped_family_counts$repFamily,]
    te_cpm_lst_1_grouped_family_counts$n <- rep(0,nrow(te_cpm_lst_1_grouped_family_counts))
    
    te_cpm_grt_lst_1_grouped_family_counts <- rbind(te_cpm_grt_1_grouped_family_counts, te_cpm_lst_1_grouped_family_counts)
    
    te_cpm_grt_lst_1_grouped_family_counts$enrich_numerator <- (te_cpm_grt_lst_1_grouped_family_counts$n)/sum(te_cpm_grt_lst_1_grouped_family_counts$n)
    te_cpm_grouped_family_counts$enrich_denominator <- (te_cpm_grouped_family_counts$n)/sum(te_cpm_grouped_family_counts$n)
    
    te_enrichment <- inner_join(te_cpm_grt_lst_1_grouped_family_counts, te_cpm_grouped_family_counts, by="repFamily")
    te_enrichment$enrichment <- (te_enrichment$enrich_numerator)/(te_enrichment$enrich_denominator)
    te_enrichment$log_enrichment <- log2(te_enrichment$enrichment)
    te_enrichment$cells <- rep(cell, nrow(te_enrichment))
    te_enrichment$samples <- rep(sample, nrow(te_enrichment))
    
    if (i == 1){
      combined_Te_enrichment <- te_enrichment
    }
    else{
      combined_Te_enrichment <- rbindlist(list(combined_Te_enrichment, te_enrichment))
    }
  }
  return(combined_Te_enrichment)
}

## load input variables to the create_matrix_file function from snakemake script

tes_cpm_path <- snakemake@input[['tes_cpm']]

output_log_enrich <- snakemake@output[['output_log_enrich']]

sample_name <- snakemake@params[['sample_name']]
output_dir <- snakemake@params[['output_dir']]
qc <- snakemake@params[['qc']]

## calculate log enrichment 

message("** Loading TE count matrix file **")

tes_cpm <- fread(tes_cpm_path, header = TRUE) %>% as.data.frame()
enrich_out <- log_enrichment(tes_cpm, sample_name) %>% as.data.frame()
fwrite(enrich_out %>% as.data.frame(),
         file= output_log_enrich,
         quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)

message("UPDATE: TE log enrichment calculation finished and saved as a .tsv file\n\n")
