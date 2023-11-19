#!/bin/bash
#SBATCH --partition=day
#SBATCH --reservation=cpsc424
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=7G
#SBATCH --job-name=rp862-final-all-serial-mm
#SBATCH --output=out/%x-%j.out
#SBATCH --time=55:00

module purge
module load intel
module list
pwd

echo $SLURM_JOB_NODELIST
echo $SLURM_NTASKS_PER_NODE
make -f Makefile-serial clean
make -f Makefile-serial

time ./bin/a-serial
time ./bin/b-serial-blocking

