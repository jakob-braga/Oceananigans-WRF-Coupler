#!/bin/bash
#SBATCH --account=def-fpoulin
#SBATCH --mail-user=jnbraga@uwaterloo.ca
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --job-name="wrf_test"
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4096
mpirun ./wrf.exe
