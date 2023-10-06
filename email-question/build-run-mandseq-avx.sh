#!/bin/bash

#SBATCH --partition=day
#SBATCH --reservation=cpsc424
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=7G
#SBATCH --job-name=Part1MandelbrotSeq
#SBATCH --output=%x-%j.out

# Load Required Modules

module load intel

# Task 1

make clean 
echo "make mandseq-avx"
make mandseq-avx

echo "Serial version with AVX"

time ./bin/mandseq-avx
time ./bin/mandseq-avx
time ./bin/mandseq-avx
