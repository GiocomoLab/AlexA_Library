#!/bin/bash
#SBATCH -p owners
#SBATCH --job-name=xcorr_connect
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=attialex@stanford.edu
#SBATCH --mem-per-cpu=80G
#SBATCH --time=12:00:00
#SBATCH -o /scratch/users/attialex/xcorr.%N.%j.out # STDOUT
#SBATCH -e /scratch/users/attialex/xcorr.%N.%j.err # STDERR

module load matlab
echo "$(date): job $SLURM_JOBID starting on $SLURM_NODELIST"

matlab -nodisplay -nosplash -r "run $HOME/AlexA_Library/correlation/connectivity_dir.m, exit"