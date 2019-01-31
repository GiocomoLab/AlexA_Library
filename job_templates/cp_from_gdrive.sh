#!/bin/bash
#SBATCH --job-name=cp_from_gdrive
#SBATCH -p giocomo
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=attialex@stanford.edu
#SBATCH --begin=now
#SBATCH --dependency=singleton
#SBATCH --signal=B:SIGUSR1@90
# catch the SIGUSR1 signal
_requeue() {
    echo "$(date): job $SLURM_JOBID received SIGUSR1, re-queueing"
    scontrol requeue $SLURM_JOBID
}
trap '_requeue' SIGUSR1
module load system rclone
echo "$(date): job $SLURM_JOBID starting on $SLURM_NODELIST"
rclone copy rclone copy VRGDrive:attialex /oak/stanford/groups/giocomo/attialex/DATA 
wait
sleep 1m
date
sleep 1m