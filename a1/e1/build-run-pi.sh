#!/bin/bash
#SBATCH --job-name=pi
#SBATCH --ntasks=1 --nodes=1 --cpus-per-task=1
#SBATCH --mem-per-cpu=7G
#SBATCH --time=10:00
#SBATCH --partition=cpsc424
#SBATCH --output=pi.out-%j

module load intel
pwd
cat $SLURMD_NODENAME
make clean
make pi
time pi
time pi
time pi
