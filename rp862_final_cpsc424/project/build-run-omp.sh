#!/bin/bash
#SBATCH --partition=day
#SBATCH --reservation=cpsc424
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=7G
#SBATCH --job-name=rp862-final-all-omp-mm
#SBATCH --output=out/%x-%j.out
#SBATCH --time=55:00

module purge
module load intel
module list
pwd

echo $SLURM_JOB_NODELIST
echo $SLURM_NTASKS_PER_NODE
make -f Makefile-omp clean
make -f Makefile-omp e-omp-ts

sizes="1024,1024,1024 8192,8192,8192"
export OMP_NUM_THREADS=24
for tuple in $sizes
do 
    IFS=',' read -ra ADDR <<< "$tuple"
    N=${ADDR[0]}
    P=${ADDR[1]}
    M=${ADDR[2]}

    echo "Running N=$N, P=$P, M=$M"
    for k in 1 2 3
    do 
        time ./bin/e-omp-ts $N $P $M
    done
done

