#!/bin/bash

#SBATCH --partition=day
#SBATCH --reservation=cpsc424
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=7G
#SBATCH --job-name=Part2a-MandelbrotOMP
#SBATCH --output=out/part2/%x-%j.out

# Load Required Modules

module load intel

echo "cleaning"
make clean-mandomp
make clean-mandomp-avx
make clean-mandomp-ts
make clean-mandomp-ts-avx

echo "compiling"
make mandomp
make mandomp-avx
make mandomp-ts
make mandomp-ts-avx

# Part 2 - 1a
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

echo ""
echo ""
echo "OMP version"

time ./bin/mandomp-avx
time ./bin/mandomp-avx
time ./bin/mandomp-avx

# Part 2 - 1b

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
        echo "NON AVX"
        time ./bin/mandomp-ts
        echo ""
        echo "AVX"
        time ./bin/mandomp-ts-avx
    done
done

# Part 2 - 2
schedules=("static,1" "static,100" "dynamic" "dynamic,250" "guided")
for schedule in "${schedules[@]}"; 
do
    for k in 2 4 12 24
    do 
        export OMP_NUM_THREADS=$k
        export OMP_SCHEDULE=$schedule

        echo ""
        echo ""
        echo "Number of threads = " $OMP_NUM_THREADS
        echo "OMP_SCHEDULE = " $OMP_SCHEDULE

        for i in {1..3}
        do
            echo "NON AVX"
            time ./bin/mandomp-ts
            echo ""
            echo "AVX"
            time ./bin/mandomp-ts-avx
        done
    done 
done
