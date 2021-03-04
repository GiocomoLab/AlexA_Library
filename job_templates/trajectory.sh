#!/bin/bash
#SBATCH --job-name=trajectory_predict
#SBATCH -p giocomo
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=attialex@stanford.edu
#SBATCH --mem-per-cpu=12G
#SBATCH -o /scratch/users/attialex/slurm.%N.%j.out # STDOUT
#SBATCH -e /scratch/users/attialex/slurm.%N.%j.err # STDERR

ml py-tensorflow/2.1.0_py36

python3 trajectory_prediction.py /oak/stanford/groups/giocomo/attialex /oak/stanford/groups/giocomo/attialex/logs/trajectories/trajectories_1000.npy