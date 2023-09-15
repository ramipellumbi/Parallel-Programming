#!/bin/bash
#SBATCH --job-name=triad
#SBATCH --ntasks=1 --nodes=1 --cpus-per-task=1
#SBATCH --mem-per-cpu=7G
#SBATCH --time=10:00
#SBATCH --partition=day
#SBATCH --reservation=cpsc424
#SBATCH --output=triad.out-%j

module load intel
pwd
echo $SLURMD_NODENAME
make clean
make triad

for k in {26..26}
do
  echo
  echo "Running with k = $k"
  time ./triad $k
done
