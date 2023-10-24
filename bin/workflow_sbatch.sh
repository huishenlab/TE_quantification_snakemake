#!/bin/bash
#SBATCH --export=NONE
#SBATCH -J SNAKEMASTER
#SBATCH -o TE_quant_workflow.o
#SBATCH -e TE_quant_workflow.e
#SBATCH --ntasks 1
#SBATCH --time 240:00:00
#SBATCH --mem=8G

cd $SLURM_SUBMIT_DIR

mkdir -p logs
mkdir -p logs/workflows

CONFIG_FILE=config/config.yaml

# Add snakemake to PATH here
if [[ `which snakemake 2>&1 /dev/null` ]]; then
    snakemake_module="bbc2/snakemake/snakemake-7.25.0"

    module load $snakemake_module
fi

# save DAG job file with time stamp
TIME=$(date "+%Y-%m-%d_%H.%M.%S")
snakemake --configfile ${CONFIG_FILE} --dry-run > logs/workflows/workflow_${TIME}.txt
snakemake --configfile ${CONFIG_FILE} --dag | dot -Tpng > logs/workflows/workflow_${TIME}.png

# Default to using conda, if using environment modules, then replace --use-conda with --use-envmodules
# Note, this requires downloading mamba (conda install -n base -c conda-forge mamba)

snakemake --unlock

snakemake \
    --printshellcmds \
    --latency-wait 20 \
    --keep-going \
    --use-conda \
    --jobs 20 \
    --configfile ${CONFIG_FILE} \
    --cluster "mkdir -p logs/{rule}; sbatch \
        --export=ALL \
        --ntasks {threads} \
        --mem={resources.mem_gb}G \
        -t {resources.time} \
        -o logs/{rule}/{rule}-%j.log"