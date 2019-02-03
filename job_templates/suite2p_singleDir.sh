#!/bin/bash
#SBATCH --job-name=suite2p
#SBATCH -p giocomo
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=attialex@stanford.edu
#SBATCH -o /scratch/users/attialex/slurm.%N.%j.out # STDOUT
#SBATCH -e /scratch/users/attialex/slurm.%N.%j.err # STDERR
#SBATCH --time=120:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=20G

echo "Hello World suite2p"
ml python/3.6.1

python3 $HOME/AlexA_Library/suite2p_helpers/run_single_file.py
