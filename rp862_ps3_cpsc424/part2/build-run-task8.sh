#!/bin/bash
#SBATCH --reservation=cpsc424
#SBATCH --partition=day

# set total number of MPI processes
#SBATCH --ntasks=7

# set number of MPI processes per node
# (number of nodes is calculated by Slurm)
#SBATCH --ntasks-per-node=2 --ntasks-per-socket=1

# set number of cpus per MPI process
#SBATCH --cpus-per-task=1

# set memory per cpu
#SBATCH --mem-per-cpu=7G
#SBATCH --job-name=T8-MPI_RUN
#SBATCH --output=out/%x-%j.out
#SBATCH --time=30:00

module purge
module load iomkl
module list
pwd

# echo some environment variables
echo "Nodes allocated:"
echo $SLURM_JOB_NODELIST
echo "Number of tasks per node:"
echo $SLURM_NTASKS_PER_NODE
echo "Number of tasks per socket:"
echo $SLURM_NTASKS_PER_SOCKET

make clean
make task8

echo " "
echo "Task 8"

# Run the application with the manager task alone on the first node and two tasks per node on the other nodes
for k in 1 2 3
do 
    # Run the MPI program with the custom hostfile
    time mpirun --report-bindings -np 7 ./bin/task8
done

echo "All Done!"
