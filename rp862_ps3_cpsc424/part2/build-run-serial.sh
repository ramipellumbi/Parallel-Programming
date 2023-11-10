#!/bin/bash
#SBATCH --partition=day
#SBATCH --reservation=cpsc424
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=7G
#SBATCH --job-name=Serial
#SBATCH --output=out/%x-%j.out
#SBATCH --time=30:00

module purge
module load iomkl
module list
pwd

echo $SLURM_JOB_NODELIST
echo $SLURM_NTASKS_PER_NODE
make clean
make serial
time ./bin/serial
time ./bin/serial
time ./bin/serial
