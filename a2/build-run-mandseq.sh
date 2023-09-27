#!/bin/bash

#SBATCH --partition=day
#SBATCH --reservation=cpsc424
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=7G
#SBATCH --job-name=MandelbrotFinalTest
#SBATCH --output=%x-%j.out

# Load Required Modules

module purge
module load intel
module list

# Task 1

make clean
echo "make mandseq"
make mandseq

echo ""
echo ""
echo "Serial version"

time ./bin/mandseq