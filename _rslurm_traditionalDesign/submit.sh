#!/bin/bash
#
#SBATCH --array=0-0
#SBATCH --cpus-per-task=6
#SBATCH --job-name=traditionalDesign
#SBATCH --output=slurm_%a.out
/usr/lib64/R/bin/Rscript --vanilla slurm_run.R
