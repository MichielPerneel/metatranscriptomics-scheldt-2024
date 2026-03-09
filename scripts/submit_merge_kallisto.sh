#!/bin/bash
#PBS -N merge_kallisto
#PBS -o data/logs/merge_kallisto_output.log
#PBS -e data/logs/merge_kallisto_error.log
#PBS -l walltime=02:00:00
#PBS -l nodes=1:ppn=16
#PBS -l mem=128gb

# Change directory to where the job was submitted from
cd $PBS_O_WORKDIR

# Activate the snakemake conda environment, which contains pandas
source activate snakemake_7

# Run the script
python scripts/merge_kallisto.py -i data/quantification -o data/quantification/