#!/bin/bash
#SBATCH -p giocomo
#SBATCH --job-name=distance_tuning
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=attialex@stanford.edu
#SBATCH --mem-per-cpu=12G
#SBATCH --time=8:00:00
#SBATCH -o /scratch/users/attialex/reg.%N.%j.out # STDOUT
#SBATCH -e /scratch/users/attialex/reg.%N.%j.err # STDERR

module load matlab
echo "$(date): job $SLURM_JOBID starting on $SLURM_NODELIST"
cd $HOME
matlab -nodisplay -nosplash -r "run ./AlexA_Library/NeuroPixel/dark_distance/distance_tuning_loop.m,exit"