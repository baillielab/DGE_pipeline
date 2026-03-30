#!/bin/bash 

# === SLURM SETTINGS ===

#SBATCH --job-name=DGE-pipeline
#SBATCH --time=96:00:00
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

# =======================


# === USER INPUT ========

# Set base paths - This is the filepath to your analysis directory:
ROOT=

# =======================


# === FILE PATHS ========
CONT=/mnt/odap-beegfs/software/singularity-images

CONT_DGE=${CONT}/dge-pipeline-container.sif
SNAKEFILE=${ROOT}/Snakefile
CONFIG=${ROOT}/config.yaml

# =======================

# === RUN PIPELINE ======
singularity exec --bind ${ROOT},${CONT} ${CONT_DGE} \
conda run -n DGE_pipeline snakemake -s ${SNAKEFILE} --configfile ${CONFIG} --cores 64 -p --use-conda --use-singularity --rerun-incomplete

# =======================

