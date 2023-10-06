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

echo ""
echo ""

# Part 2 - 1a

# Part 2 - 1b

# Part 2 - 2
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
        time ./bin/mandomp-ts-avx
    done
done

# Part 2 - 3
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
            time ./bin/mandomp-ts-avx
        done
    done 
done
