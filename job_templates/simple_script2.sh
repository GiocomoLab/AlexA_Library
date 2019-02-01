#!/bin/bash#
#SBATCH -p giocomo
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 100 # memory pool for all cores
#SBATCH -t 0-2:00 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
for i in {1..100000}; do
echo $RANDOM >> SomeRandomNumbers.txt
donesort SomeRandomNumbers.txt
