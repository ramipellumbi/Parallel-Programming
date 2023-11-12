#!/bin/bash
#SBATCH --reservation=cpsc424
#SBATCH --partition=day
#SBATCH --job-name=T8-MPI_RUN
#SBATCH --output=out/%x-%j.out
#SBATCH --time=30:00

# Request 4 nodes in total
#SBATCH --ntasks-per-node=4 --ntasks-per-socket=2
# Request 7 tasks in total, which will include 1 manager and 6 workers
#SBATCH --ntasks=7
# Request 1 CPU per task
#SBATCH --cpus-per-task=1
# Request 7 GB of memory per CPU
#SBATCH --mem-per-cpu=7G
# Use exclusive mode to prevent sharing of nodes between jobs
#SBATCH --exclusive

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

# Expand the nodelist using scontrol and create a hostfile
scontrol show hostname $SLURM_JOB_NODELIST > hostfile

# Edit the hostfile to specify the number of slots:
# The manager node (first line in the hostfile) should have 1 slot.
# The other nodes should have 2 slots.
MANAGER_NODE=$(head -n 1 hostfile)
WORKER_NODES=$(tail -n +2 hostfile)

# Create a new hostfile for mpirun
echo "$MANAGER_NODE slots=1" > my_hostfile
for NODE in $WORKER_NODES; do
    echo "$NODE slots=2" >> my_hostfile
done

# Run the application with the manager task alone on the first node and two tasks per node on the other nodes
for k in 1 
do 
    # Run the MPI program with the custom hostfile
    time mpirun --report-bindings -np 7 ./bin/task8
    # time mpirun --report-bindings --hostfile my_hostfile -np 7 --map-by ppr:2:node ./bin/task8
done

echo "All Done!"
