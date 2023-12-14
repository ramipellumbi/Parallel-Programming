#!/bin/bash
#SBATCH --partition=day
#SBATCH --reservation=cpsc424
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=45G
#SBATCH --job-name=rp862-final-all-serial-mm
#SBATCH --output=out/1-serial/%x-%j.out
#SBATCH --time=55:00

module purge
module load intel
module list
pwd

echo $SLURM_JOB_NODELIST
echo $SLURM_NTASKS_PER_NODE
make -f Makefile-serial clean
make -f Makefile-serial

# N,P,M sizes delimited by space
ijk_sizes="1024,1024,1024, 2048,2048,2048"
sizes="1024,1024,1024 8192,8192,8192 1024,1024,8192 8192,8192,1024 8192,1024,8192"

for tuple in $sizes
do 
    IFS=',' read -ra ADDR <<< "$tuple"
    N=${ADDR[0]}
    P=${ADDR[1]}
    M=${ADDR[2]}

    echo "Running N=$N, P=$P, M=$M"
    for k in {1..3}
    do
        time ./bin/t2-serial-kij $N $P $M
        time ./bin/t3-serial-blocking $N $P $M
    done
done

for tuple in $ijk_sizes
do 
    IFS=',' read -ra ADDR <<< "$tuple"
    N=${ADDR[0]}
    P=${ADDR[1]}
    M=${ADDR[2]}

    echo "Running N=$N, P=$P, M=$M"
    time ./bin/t1-serial-ijk $N $P $M
    time ./bin/t1-serial-ijk $N $P $M
    time ./bin/t1-serial-ijk $N $P $M
done
