#!/bin/bash

#SBATCH --partition=day
#SBATCH --reservation=cpsc424
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=7G
#SBATCH --job-name=Part2a-1a-MandelbrotOmp
#SBATCH --output=out/part2-1a/%x-%j.out

# Load Required Modules

module load intel

# Task 2 part 1a

make clean-mandomp
echo "make mandomp"
make mandomp 

echo ""
echo ""
echo "OMP version"
export OMP_NUM_THREADS=2
unset OMP_SCHEDULE

echo "Number of threads = " $OMP_NUM_THREADS
echo "OMP_SCHEDULE = " $OMP_SCHEDULE

time ./bin/mandomp
time ./bin/mandomp
time ./bin/mandomp