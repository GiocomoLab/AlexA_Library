#!/bin/bash
#SBATCH -p owners
#SBATCH --job-name=rigid_reg
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=attialex@stanford.edu
#SBATCH --mem-per-cpu=8G
#SBATCH --time=12:00:00
#SBATCH -o /scratch/users/attialex/reg.%N.%j.out # STDOUT
#SBATCH -e /scratch/users/attialex/reg.%N.%j.err # STDERR

module load matlab
echo "$(date): job $SLURM_JOBID starting on $SLURM_NODELIST"

matlab -nodisplay -nosplash -r "run $HOME/AlexA_Library/NeuroPixel/speed_sorting/loop_files_speed.m, exit"