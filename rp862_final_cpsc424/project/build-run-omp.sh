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
make -f Makefile-omp c-omp

export OMP_NUM_THREADS=24

time ./bin/c-omp

