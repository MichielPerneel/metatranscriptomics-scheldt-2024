#!/bin/bash

##################
# Inspiration: https://github.com/SchlossLab/snakemake_cluster_tutorial
####  PBS settings

####  Job name
#PBS -N submit_snakemake

####  Resources

#PBS -l nodes=1:ppn=1,mem=16000mb
#PBS -l walltime=72:00:00

####  Account and return

#PBS -M name@ugent.be
#PBS -m ae

#PBS -o /?/

#### Load Snakemake conda environment
source activate snakemake_7

##################

##  Change to the directory from which you submit the job, if running
##  from within a job
if [ -d "$PBS_O_WORKDIR" ] ; then
    cd $PBS_O_WORKDIR
fi

# Initiating snakemake and running workflow in cluster mode
snakemake --profile hpc_config/ --rerun-incomplete --use-conda --conda-prefix envs/ --conda-frontend conda --rerun-triggers mtime