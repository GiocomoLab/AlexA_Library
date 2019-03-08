#!/bin/bash
#SBATCH --job-name=cp_from_gdrive
#SBATCH -p giocomo
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=attialex@stanford.edu
#SBATCH --mem-per-cpu=4G
#SBATCH -o /scratch/users/attialex/slurm.%N.%j.out # STDOUT
#SBATCH -e /scratch/users/attialex/slurm.%N.%j.err # STDERR

ml python/3.6.1

python3 $HOME/AlexA_Library/suite2p_helpers/run_single_file.py