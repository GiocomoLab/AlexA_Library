#!/bin/bash
#SBATCH --job-name=cp_from_gdrive
#SBATCH -p giocomo
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=attialex@stanford.edu
#SBATCH --mem 100 # memory pool for all cores

module load system rclone
echo "$(date): job $SLURM_JOBID starting on $SLURM_NODELIST"
rclone copy rclone copy VRGDrive:attialex /oak/stanford/groups/giocomo/attialex/DATA -P
