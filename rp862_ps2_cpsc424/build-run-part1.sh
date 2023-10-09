#!/bin/bash

#SBATCH --partition=day
#SBATCH --reservation=cpsc424
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=7G
#SBATCH --job-name=Part1MandelbrotSeq
#SBATCH --output=out/part1/%x-%j.out

# Load Required Modules

module load intel

# Task 1

make clean-mandseq
make clean-mandseq-avx
echo "make mandseq"
make mandseq
make mandseq-avx

echo ""
echo ""
echo "Serial version"

for i in {1..3}
do 
    time ./bin/mandseq
done

echo ""
echo ""
echo "Serial version AVX"

for i in {1..3}
do 
    time ./bin/mandseq-avx
done