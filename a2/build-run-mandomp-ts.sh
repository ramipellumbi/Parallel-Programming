#!/bin/bash

#SBATCH --partition=day
#SBATCH --reservation=cpsc424
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=7G
#SBATCH --job-name=Part2a-1b-MandelbrotOmpThreadSafe
#SBATCH --output=out/part2-1b/%x-%j.out

# Load Required Modules

module load intel

# Task 2 part 1b

make clean
echo "make mandomp-ts"
make mandomp-ts 

echo ""
echo ""
echo "OMP THREAD SAFE version"

for k in 1 2 4 12 24
do
    export OMP_NUM_THREADS=$k
    unset OMP_SCHEDULE

    echo ""
    echo ""
    echo "Number of threads = " $OMP_NUM_THREADS
    echo "OMP_SCHEDULE = " $OMP_SCHEDULE

    for i in {1..3}
    do
        time ./bin/mandomp-ts
    done
done
