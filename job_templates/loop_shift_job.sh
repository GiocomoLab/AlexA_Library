#!/bin/bash
#SBATCH -p giocomo
#SBATCH --job-name=find_shifts
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=attialex@stanford.edu
#SBATCH --mem-per-cpu=4G
#SBATCH --time=4:00:00
#SBATCH -o /scratch/users/attialex/reg.%N.%j.out # STDOUT
#SBATCH -e /scratch/users/attialex/reg.%N.%j.err # STDERR

module load matlab
echo "$(date): job $SLURM_JOBID starting on $SLURM_NODELIST"
cd $HOME
matlab -nodisplay -nosplash -r "run ./AlexA_Library/NeuroPixel/speed_sorting/loop_findshifts_gaincontrast.m,exit"