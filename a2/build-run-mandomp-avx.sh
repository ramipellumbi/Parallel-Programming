#!/bin/bash

#SBATCH --partition=day
#SBATCH --reservation=cpsc424
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=7G
#SBATCH --job-name=Part2a-1a-MandelbrotOmp
#SBATCH --output=out/part2-1a-avx/%x-%j.out

# Load Required Modules

module load intel

# Task 2 part 1a AVX

make clean-mandomp-avx
echo "make mandomp-avx"
make mandomp-avx

echo ""
echo ""
echo "OMP version AVX"
export OMP_NUM_THREADS=2
unset OMP_SCHEDULE

echo "Number of threads = " $OMP_NUM_THREADS
echo "OMP_SCHEDULE = " $OMP_SCHEDULE

time ./bin/mandomp-avx