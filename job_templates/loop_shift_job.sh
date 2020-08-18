#!/bin/bash
#SBATCH -p giocomo
#SBATCH --job-name=find_shifts
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=attialex@stanford.edu
#SBATCH --mem-per-cpu=5G
#SBATCH --time=12:00:00
#SBATCH -o /scratch/users/attialex/reg.%N.%j.out # STDOUT
#SBATCH -e /scratch/users/attialex/reg.%N.%j.err # STDERR

module load matlab
echo "$(date): job $SLURM_JOBID starting on $SLURM_NODELIST"

matlab -nodisplay -nosplash -r "run $HOME/startup.m,$HOME/AlexA_Library/NeuroPixel/speed_sorting/loop_findshifts_v2.m,exit"