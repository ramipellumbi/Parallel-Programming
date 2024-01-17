#!/bin/bash

#SBATCH --partition=day
#SBATCH --reservation=cpsc424
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=7G
#SBATCH --job-name=Part3-MandelbrotOMPTasks
#SBATCH --output=out/part3/%x-%j.out

# Load Required Modules

module load intel

echo "cleaning"
make clean-mandomp-tasks
make clean-mandomp-tasks-columns
make clean-mandomp-tasks-columns-shared
echo ""

echo "compiling"
make mandomp-tasks
make mandomp-tasks-columns
make mandomp-tasks-columns-shared

echo ""
echo ""
echo "MANDOMP TASKS SINGLE CELL"

for k in 1 2 4 12 24
do
    export OMP_NUM_THREADS=$k
    unset OMP_SCHEDULE

    echo ""
    echo ""
    echo "Number of threads = " $OMP_NUM_THREADS

    for i in {1..3}
    do
        echo ""
        echo "Per cell"
        time ./bin/mandomp-tasks
        echo ""
        echo "Per column"
        time ./bin/mandomp-tasks-columns
        echo ""
        echo "Shared Task Creation"
        time ./bin/mandomp-tasks-columns-shared
    done
done