#!/bin/bash
#SBATCH --partition=day
#SBATCH --reservation=cpsc424
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem-per-cpu=7G
#SBATCH --time=20:00
#SBATCH --job-name=MPI_Hello
#SBATCH --output=%x-%j.out

# Load necessary module files
module load iomkl
module list

# Print initial working directory
echo " "
echo " "
echo "Working Directory:"
pwd

echo " "
echo " "
echo "Making task1"
make clean
make task1

# Print the node list
echo " "
echo " "
echo "Node List:"
echo $SLURM_NODELIST
echo "ntasks-per-node = " $SLURM_NTASKS_PER_NODE

# Run the program 3 times
echo " "
echo " "
echo "Run 1"
time mpirun -n 4 task1
echo " "
echo " "
echo "Run 2"
time mpirun -n 4 task1
echo " "
echo " "
echo "Run 3"
time mpirun -n 4 task1
