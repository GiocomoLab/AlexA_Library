#!/bin/bash
#SBATCH --job-name=cp_from_gdrive
#SBATCH -p giocomo
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=attialex@stanford.edu
#SBATCH --mem-per-cpu=40G
#SBATCH -o /scratch/users/attialex/slurm.%N.%j.out # STDOUT
#SBATCH -e /scratch/users/attialex/slurm.%N.%j.err # STDERR

module load matlab
echo "$(date): job $SLURM_JOBID starting on $SLURM_NODELIST"

matlab -nodisplay -nosplash -r "run $HOME/AlexA_Library/suite2p_helpers/create_h5_files.m, exit"