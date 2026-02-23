#!/bin/bash
#SBATCH -p ndl
#SBATCH --gres=gpu:4
#SBATCH -N 1
#SBATCH --verbose
##SBATCH --exclusive
#SBATCH --time=00:02:00
##SBATCH --output=sortie.txt
##SBATCH --switches=3
##SBATCH -c 2

module purge
module load nvhpc/24.5

set -x

srun -o affichage_top.txt --gres=gpu:4 testbatched_top.x
srun -o affichage.txt --gres=gpu:4 testbatched.x
