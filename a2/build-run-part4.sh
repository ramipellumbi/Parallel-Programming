#!/bin/bash

#SBATCH --partition=day
#SBATCH --reservation=cpsc424
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=7G
#SBATCH --job-name=Part4-MandelbrotOMPParallel
#SBATCH --output=out/part4/%x-%j.out

# Load Required Modules

module load intel

echo "cleaning"
make clean-mandomp-ts-avx-parallel

echo "compiling"
make mandomp-ts-avx-parallel

export OMP_NUM_THREADS=24
export OMP_SCHEDULE="guided"

time ./bin/mandomp-ts-avx-parallel
time ./bin/mandomp-ts-avx-parallel
time ./bin/mandomp-ts-avx-parallel
