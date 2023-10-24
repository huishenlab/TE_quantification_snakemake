# TE quantification snakemake workflow

## Overview

This is a pipeline for the quantification of Transposable Elements in the single-cell STORM-seq samples. It can also be used with other technologies like SMART-seq, SMART-seq2 etc. 

## Dependencies

### Command Line tools

1) featureCounts

### R libraries

1) data.table
2) dplyr
3) stringr
4) optparse
5) ggplot2
6) scales
7) ggh4x
8) Scuttle

## Steps in the pipeline

1) Aligned bam files are used as inputs to featureCounts with parameters `-F SAF -O -B -p --fracOverlap 0.1 -M -s 0 --fraction` specified in the config.yaml file. These can also be altered as per the requirement.
2) feature counts files for all the cells are then used to generate a combined raw count, counts per million, raw count for only intergenic and intronic TEs and counts per millions for only intergenic and intronic TEs matrices.
3) If filtering requirement is set to be True in the `config.yaml` file then a quick Scuttle filtering is performed to filter out low quality cells and generate filtered count matrices.
4) Counts per millions for only intergenic and intronic TEs matrix is used for the log enrichment calculation of TE's
6) Enrichment score is calculated as per the folloing formula.

$$ enrichment\ score = {{{Number\ of\ TE\ subfamilies\ >\ 1cpm} \over {Number\ of\ TEs\ >\ 1cpm}} \over {{Number\ of\ TE\ subfamilies} \over {Number\ of\ TEs}}} $$

7) If the log enrichment heatmap plot requirement is set to be True in the `config.yaml` file then a heatmap plot, generated using ggplot2, is also saved as a pdf file.
   
## Using the workflow

1) Cone the repo (https://github.com/AyushSemwal/TE_quantification_snakemake) using:
   *  HTTPS: `https://github.com/AyushSemwal/TE_quantification_snakemake.git`
   *  SSH: `git@github.com:AyushSemwal/TE_quantification_snakemake.git`
2) Unzip `hg38_pc_te_chrM.saf.tar.gz` and `intergenic_intronic_tes.txt.tar.gz` in the config folder.
3) Modify the `config.yaml` file in the config folder to specify `aligned_bam_dir` (aligned bam files directory), `output_dir` (directory where all output files and sub directories will be stored) and other parameters as per the requirements. Even though I have assigned intuitive names to the parameters, I have also added comments in front of them.
4) Populate the `samples.tsv` file such that first column is contains cell names (or sample names in case of bulk samples) and second columns contains the bam file names. Do not add a header to this file.
5) If slurm is available then submit the job by running `sbatch bin/workflow_sbatch.sh` from the parent directory else you can run `snakemake --use-conda --cores {num_cores}`.
