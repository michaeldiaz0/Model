#!/bin/bash
# Job name:
#SBATCH --job-name=test
#
# Account:
#SBATCH --account=fc_tropical
#
# Partition:
#SBATCH --partition=savio
#
# Wall clock limit:
#SBATCH --time=00:05:00
#

#SBATCH --nodes=4

#SBATCH --ntasks-per-node=20

#SBATCH --cpus-per-task=1

module load gcc/4.8.5
module load netcdf
module load fftw
module load openmpi

mpirun -np 80 /global/home/users/michaeldiaz/Model/solve.exe