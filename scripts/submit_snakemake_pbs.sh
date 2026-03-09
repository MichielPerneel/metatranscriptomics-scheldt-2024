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

#PBS -M nadine.rijsdijk@ugent.be
#PBS -m ae

#PBS -o /scratch/gent/vo/000/gvo00081/metatranscriptomics-workflow-michiel-nadine/submit_snakemake.$PBS_JOBID.out
#PBS -e /scratch/gent/vo/000/gvo00081/metatranscriptomics-workflow-michiel-nadine/submit_snakemake.$PBS_JOBID.err

#### Load Snakemake conda environment
source activate snakemake_7

##################

##  Change to the directory from which you submit the job, if running
##  from within a job
if [ -d "$PBS_O_WORKDIR" ] ; then
    cd $PBS_O_WORKDIR
fi

mkdir -p /kyukon/scratch/gent/vo/000/gvo00081/{tmp,conda_pkgs,cache,pip_cache}/$USER

export TMPDIR=/kyukon/scratch/gent/vo/000/gvo00081/tmp/$USER
export TEMP=$TMPDIR
export TMP=$TMPDIR
export CONDA_PKGS_DIRS=/kyukon/scratch/gent/vo/000/gvo00081/conda_pkgs/$USER
export XDG_CACHE_HOME=/kyukon/scratch/gent/vo/000/gvo00081/cache/$USER
export PIP_CACHE_DIR=/kyukon/scratch/gent/vo/000/gvo00081/pip_cache/$USER

# Initiating snakemake and running workflow in cluster mode
snakemake --profile hpc_config/ --rerun-incomplete --use-conda --conda-prefix envs/ --conda-frontend conda --rerun-triggers mtime --directory .