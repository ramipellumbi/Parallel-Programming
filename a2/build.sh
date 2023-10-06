#!/bin/bash

#SBATCH --partition=day
#SBATCH --reservation=cpsc424
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=7G
#SBATCH --job-name=BUILD-MANDELBROT
#SBATCH --output=%x-%j.out

module load intel

make clean
make mandomp-avx
make mandomp-ts-avx
make mandomp-ts
make mandomp
make mandseq-avx
make mandseq