#!/bin/bash
#
#SBATCH --array=0-0
#SBATCH --cpus-per-task=6
#SBATCH --job-name=SequentialDesignWithtLimit
#SBATCH --output=slurm_%a.out
/imaging/local/software/R/3.5.3shlib/bin/Rscript --vanilla slurm_run.R
