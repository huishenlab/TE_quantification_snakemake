# Import modules
import os

from snakemake.utils import min_version

# Set minimum Snakemake version
min_version('7.0')

##########################################
########### Load configuration ###########
##########################################

configfile: "config/config.yaml"

#######################################################
## following function specifies the output directory ##
#######################################################

def set_output_directory():
    if config['output_dir'] == '':
        return os.getcwd()
    else:
        return config['output_dir']

output_dir = set_output_directory()

#############################################
########### assign some variables ###########
#############################################

featureCounts_threads = config["resource_specs"]["featureCounts_threads"]
fread_threads = config["resource_specs"]["fread_threads"]

aligned_bam_dir = config["aligned_bam_dir"]

samples_tsv = config["samples"]["samples_tsv"]
single_cell = config["samples"]["single_cell"]

aligned_bam_dir = config["aligned_bam_dir"]
annotation_saf = config["featurecounts"]["annotation_saf"]

intronic_intergenic_tes = config['create_matrix']["intronic_intergenic_tes"]

#########################################################
## create symlinks of the imput bams in the output dir ##
#########################################################

SAMPLES = []

def create_symlinks():
    os.system("mkdir -p " + output_dir + "/" + "aligned_files")
    samples_file = open(samples_tsv,'r')
    for line in samples_file:
        fields = line.rstrip().split("\t")
        sample = fields[0]
        SAMPLES.append(sample)
        os.system("ln -s " + aligned_bam_dir + "/" + fields[1] + " " + output_dir + "/" + "aligned_files/" + sample + ".bam")

create_symlinks()

#########################################################################################################
## depending on qc requirement the following function specifies the count matrix files to be generated ##
#########################################################################################################

def count_matrix():
    if config['samples']['qc'] == "TRUE":
        return ["count_matrix_filtered_TE.tsv", "count_matrix_filtered_cpm_TE.tsv", "count_matrix.tsv", "count_matrix_TE.tsv", "count_matrix_cpm_TE.tsv"]
    else:
        return ["count_matrix.tsv", "count_matrix_TE.tsv", "count_matrix_cpm_TE.tsv"]

count_matrix_files = count_matrix()

##################################################################################################################
## depending on qc requirement the following function specifies the te's cpm file to be used for log enrichment ##
##################################################################################################################

def tes_cpm():
    if config['samples']['qc'] == "TRUE":
        return "count_matrix_filtered_cpm_TE.tsv"
    else:
        return "count_matrix_cpm_TE.tsv"
    
tes_cpm_file = tes_cpm()

######################################################################################################
## depending on qc requirement the following function specifies the name of the log enrichment file ##
######################################################################################################

def log_enrich():
    if config['samples']['qc'] == "TRUE":
        return "filtered_log_enrichment.tsv"
    else:
        return "log_enrichment.tsv"
    
log_enrich_file = log_enrich()

######################################################################################################
## depending on qc requirement the following function specifies the name of the log enrichment file ##
######################################################################################################

def log_enrich_plot():
    return "log_enrichment_plot.pdf"
    
log_enrich_plot_pdf = log_enrich_plot()

###############################################
## all the output files from different rules ##
###############################################

rule all:
    input:
        expand('{output_dir}/{log_enrich_plot_pdf}', log_enrich_plot_pdf=log_enrich_plot_pdf, output_dir=output_dir) if config['log_enrichment']['plot'] else [],
        expand('{output_dir}/{log_enrich_file}', log_enrich_file=log_enrich_file, output_dir=output_dir) if config['log_enrichment']['calculate'] else [],
        expand('{output_dir}/{count_matrix}', count_matrix=count_matrix_files, output_dir=output_dir),
        expand('{output_dir}/featureCounts/{sample}_fc.txt', sample=SAMPLES, output_dir=output_dir)

######## include rules ##########

include: "rules/featurecounts.smk"
include: "rules/create_matrix_file.smk"
include: "rules/log_enrichment.smk"
include: "rules/log_enrichment_plot.smk"

