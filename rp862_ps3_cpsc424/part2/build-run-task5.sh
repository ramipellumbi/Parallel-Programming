#!/bin/bash
#SBATCH --reservation=cpsc424
#SBATCH --partition=day

# set total number of MPI processes
#SBATCH --ntasks=8

# set number of MPI processes per node
# (number of nodes is calculated by Slurm)
#SBATCH --ntasks-per-node=4 --ntasks-per-socket=2

# set number of cpus per MPI process
#SBATCH --cpus-per-task=1
#SBATCH --distribution=cyclic:cyclic

# set memory per cpu
#SBATCH --mem-per-cpu=7G
#SBATCH --job-name=MPI_RUN
#SBATCH --output=out/%x-%j.out
#SBATCH --time=30:00

module purge
module load iomkl
module list
pwd

# echo some environment variables
echo $SLURM_JOB_NODELIST
echo $SLURM_NTASKS_PER_NODE
echo $SLURM_NTASKS_PER_SOCKET

make clean
make task5

echo " "
echo "Task 5"
for ntasks in 2 4 8
do 
    for k in 1 2 3 
    do 
        time mpirun --report-bindings -np ${ntasks} ./bin/task5 ${ntasks}
    done
done

echo "All Done!"
