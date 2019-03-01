#!/bin/bash
#SBATCH --job-name=cp_bin
#SBATCH -p owners
#SBATCH --time=04:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=attialex@stanford.edu
#SBATCH --mem-per-cpu=2G
#SBATCH -o /scratch/users/attialex/copy/slurm.%N.%j.out # STDOUT
#SBATCH -e /scratch/users/attialex/copy/slurm.%N.%j.err # STDERR

module load system rclone
echo "$(date): job $SLURM_JOBID starting on $SLURM_NODELIST"
rclone --include "data.bin" copy $SCRATCH/register_jobs/bin_files $OAK/attialex/REGISTRATION
rclone check $SCRATCH/register_jobs/bin_files $OAK/attialex/REGISTRATION --size-only --one-way