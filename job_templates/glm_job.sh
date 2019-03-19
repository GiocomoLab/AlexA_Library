#!/bin/bash
#SBATCH -p giocomo
#SBATCH --job-name=glmModel
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=attialex@stanford.edu
#SBATCH --mem-per-cpu=4G
#SBATCH --time=7:00:00
#SBATCH -o /scratch/users/attialex/glm.%N.%j.out # STDOUT
#SBATCH -e /scratch/users/attialex/gml.%N.%j.err # STDERR

module load matlab
echo "$(date): job $SLURM_JOBID starting on $SLURM_NODELIST"

matlab -nodisplay -nosplash -r "run $HOME/AlexA_Library/spline-lnp-model-VR/glmModel_dir.m, exit"