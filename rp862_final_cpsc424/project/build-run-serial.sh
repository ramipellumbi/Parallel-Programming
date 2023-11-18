#!/bin/bash
#SBATCH --partition=day
#SBATCH --reservation=cpsc424
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=7G
#SBATCH --job-name=rp862-final-all-serial-mm
#SBATCH --output=out/%x-%j.out
#SBATCH --time=30:00

module purge
module load iomkl
module list
pwd

echo $SLURM_JOB_NODELIST
echo $SLURM_NTASKS_PER_NODE
make -f Makefile-serial clean
make -f Makefile-serial

for k in 1 2 3
do 
    time ./bin/a-serial
    time ./bin/b-serial-avx
done
