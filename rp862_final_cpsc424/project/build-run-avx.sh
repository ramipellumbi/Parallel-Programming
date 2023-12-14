#!/bin/bash
#SBATCH --partition=day
#SBATCH --reservation=cpsc424
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=7G
#SBATCH --job-name=rp862-final-all-avx-mm
#SBATCH --output=out/3-avx/%x-%j.out
#SBATCH --time=55:00

module purge
module load intel
module list
pwd

echo $SLURM_JOB_NODELIST
echo $SLURM_NTASKS_PER_NODE
make -f Makefile-avx clean
make -f Makefile-avx t9-omp-divisible-avx-blocking

# N,P,M sizes delimited by space
sizes="1200,1200,1200 3600,3600,3600 1440,1440,1440 4320,4320,4320 8400,8400,8400 9600,9600,9600"
export OMP_NUM_THREADS=24
for tuple in $sizes
do 
    IFS=',' read -ra ADDR <<< "$tuple"
    N=${ADDR[0]}
    P=${ADDR[1]}
    M=${ADDR[2]}

    echo "Running N=$N, P=$P, M=$M"
    for k in {1..3}
    do
        time ./bin/t9-omp-divisible-avx-blocking $N $P $M
    done
done
