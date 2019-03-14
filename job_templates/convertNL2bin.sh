#!/bin/bash
#SBATCH -p owners
#SBATCH --job-name=nl2bin
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=attialex@stanford.edu
#SBATCH --mem-per-cpu=50G
#SBATCH --time=2:00:00
#SBATCH -o /scratch/users/attialex/nl2bin.%N.%j.out # STDOUT
#SBATCH -e /scratch/users/attialex/nl2bin.%N.%j.err # STDERR

module load matlab
echo "$(date): job $SLURM_JOBID starting on $SLURM_NODELIST"

matlab -nodisplay -nosplash -r "run $HOME/AlexA_Library/Neuropixel/convertNL2binScript.m, exit"