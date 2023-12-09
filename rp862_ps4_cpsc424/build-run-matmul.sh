#!/bin/bash

#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=45G
#SBATCH --partition=gpu
#SBATCH --reservation=cpsc424gpu
#SBATCH -t 10:00
#SBATCH --gpus=1

echo "***Purging module files"
echo ""
module purge
echo ""
echo "***Loading CUDA module file"
echo ""
module load CUDA
echo ""
module list

echo ""
echo "***Running nvidia-smi"
echo ""
nvidia-smi
echo ""
echo ""

echo "***Running deviceQuery"
/vast/palmer/apps/avx.grace/software/CUDAcore/11.3.1/extras/demo_suite/deviceQuery
echo ""

echo "***Building matmul"
make clean
make matmul

# Now run the code. Note that if you turn on the error check using a
# cpu matmul code to check the answers, you will need more time for
# the job (possibly as much as 2 hours if you run all 4 test cases)
echo ""
echo "***Running matmul"
time ./matmul 1024 32 32
echo ""
time ./matmul 2048 32 32
echo ""
time ./matmul 4096 32 32
echo ""
time ./matmul 8192 32 32
echo ""
echo "***All Done."
