#!/usr/bin/env Rscripts

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(optparse))

options(scipen=999)

create_matrix_file <- function(sample_names, samples_fc, single_cell, intronic_intergenic_tes, output_dir, qc, threads){

  cells <- sample_names
  fc_files <- samples_fc
  single_cell <- single_cell

  ## Create a combined matrix file from individual feature counts files

  message("** Combining feature counts to create a matrix file **")
  
  for (i in 1:length(cells)){
    cell <- cells[i]
    
    if (i == 1){
      feat_counts_file <- fread(samples_fc[i], select = c(1,7), 
                                nThread = threads)
      names(feat_counts_file)[ncol(feat_counts_file)] <- cell
      count_matrix <- feat_counts_file
    }
    else {
      feat_counts_file <- fread(samples_fc[i], select = c(7), 
                                nThread = threads)
      names(feat_counts_file)[ncol(feat_counts_file)] <- cell
      count_matrix[, cell] <- feat_counts_file[, ..cell]
    }
  }
  
  Geneids <- count_matrix$Geneid
  
  colnames(count_matrix) <- c("Geneid", cells)
  fwrite(count_matrix %>% as.data.frame(),
         file=paste0(output_dir, "/count_matrix.tsv"),
         quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)
  
  count_matrix <- count_matrix[,2:ncol(count_matrix)]

  message("UPDATE: Raw count matrix file generated\n\n")

  ## Compute cpm and extract TE's 

  message("** Computing CPM, extracting TE's and saving raw counts and cpm matrices for TE's **")
  
  intergenic_intronic_tes_df <- fread(intronic_intergenic_tes)  
  colnames(intergenic_intronic_tes_df) <- gsub(colnames(intergenic_intronic_tes_df)[length(colnames(intergenic_intronic_tes_df)) - 1], "Geneid", colnames(intergenic_intronic_tes_df))
  
  count_matrix$Geneid <- Geneids
  count_matrix_TE <- inner_join(count_matrix, intergenic_intronic_tes_df[,c(6,7,8,9)], by = "Geneid")
  
  count_matrix <- count_matrix[, 1:(ncol(count_matrix)-1)]
  count_matrix_cpm <- sweep(count_matrix, 2, colSums(count_matrix)/1000000, `/`)
  count_matrix_cpm$Geneid <- Geneids
  
  count_matrix_cpm_TE <- inner_join(count_matrix_cpm, intergenic_intronic_tes_df[,c(6,7,8,9)], by = "Geneid")
  
  fwrite(count_matrix_TE %>% as.data.frame(),
         file=paste0(output_dir, "/count_matrix_TE.tsv"),
         quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)
  
  fwrite(count_matrix_cpm_TE %>% as.data.frame(),
         file=paste0(output_dir, "/count_matrix_cpm_TE.tsv"),
         quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)

  message("UPDATE: CPM computation and TE extraction finished\n\n")

## perform quick quality control if the samples are single-cell

  if (single_cell == "TRUE"){
    if (qc == "TRUE"){

      message("** Performing quick QC using Scuttle and saving raw TE counts and CPM for filtered cells  **")

      qcstats <- scuttle::perCellQCMetrics(count_matrix) 
      qcfilter <- scuttle::quickPerCellQC(qcstats)
      count_matrix_filtered <- count_matrix[,!qcfilter$discard, with=FALSE]
      
      count_matrix_filtered$Geneid <- Geneids
      count_matrix_filtered_TE <- inner_join(count_matrix_filtered, intergenic_intronic_tes_df[,c(6,7,8,9)], by = "Geneid")
      
      count_matrix_filtered <- count_matrix_filtered[, 1:(ncol(count_matrix_filtered)-1)]
      count_matrix_filtered_cpm <- sweep(count_matrix_filtered, 2, colSums(count_matrix_filtered)/1000000, `/`)
      count_matrix_filtered_cpm$Geneid <- Geneids
      
      count_matrix_filtered_cpm_TE <- inner_join(count_matrix_filtered_cpm, intergenic_intronic_tes_df[,c(6,7,8,9)], by = "Geneid")
    }
    
    fwrite(count_matrix_filtered_TE %>% as.data.frame(),
           file=paste0(output_dir, "/count_matrix_filtered_TE.tsv"),
           quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)
    
    fwrite(count_matrix_filtered_cpm_TE %>% as.data.frame(),
           file=paste0(output_dir, "/count_matrix_filtered_cpm_TE.tsv"),
           quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)

    message("UPDATE: QC and filtering finished\n\n")
  }
}

## load input variables to the create_matrix_file function from snakemake script

sample_names <- snakemake@params[["sample_names"]]
samples_fc <- snakemake@input[["samples_fc"]]
intronic_intergenic_tes <- snakemake@input[["intronic_intergenic_tes"]]

single_cell <- snakemake@params[["single_cell"]]
output_dir <- snakemake@params[["output_dir"]]
qc <- snakemake@params[["qc"]]
threads <-snakemake@params[["threads"]]

## generate counts matrices

create_matrix_file(sample_names, samples_fc, single_cell, intronic_intergenic_tes, output_dir, qc, threads)
