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
#SBATCH --output=%x-%j.out
#SBATCH --time=30:00

module purge
module load iomkl
module list
pwd
# echo some environment variables
echo $SLURM_JOB_NODELIST
echo $SLURM_NTASKS_PER_NODE
# Do a clean build
make -f Makefile-mpi clean
# My MPI program is task5, task6, task7 or task8
make -f Makefile-mpi task5
#make -f Makefile-mpi task6
#make -f Makefile-mpi task7
#make -f Makefile-mpi task8
# The following mpirun command will pick up required info on nodes and cpus from Slurm. 
# You can use mpirun's -n option to reduce the number of MPI processes started on the cpus. (At most 1 MPI proc per Slurm task.)
# You can use mpirun options to control the layout of MPI processes---e.g., to spread processes out onto multiple nodes
# In this example, we've asked Slurm for 4 tasks (2 each on 2 nodes), but we've asked mpirun for two MPI procs, which will go onto 1 node.
# (If "-n 2" is omitted, you'll get 4 MPI procs---1 per Slurm task)
echo " "
echo "Task 5"
time mpirun --report-bindings ./task5
#echo " "
#echo "Task 6"
#time mpirun --report-bindings ./task6
#echo " "
#echo "Task 7"
#time mpirun --report-bindings ./task7
#echo " "
#echo "Task 8"
#time mpirun --report-bindings ./task8
echo "All Done!"
