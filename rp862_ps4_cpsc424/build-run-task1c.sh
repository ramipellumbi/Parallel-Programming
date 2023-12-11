#!/bin/bash

#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=45G
#SBATCH --partition=gpu
#SBATCH --reservation=cpsc424gpu
#SBATCH -t 120:00
#SBATCH --gpus=1
#SBATCH --job-name=T1a
#SBATCH --output=out/task1/%x-%j.out


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

echo "***Building task 1"
make clean
make serial
make task1 

# PART B - EDIT FP to double first and comment out the above

sizes="30536,30536,30536"
for tuple in $sizes
do
    IFS=',' read -ra ADDR <<< "$tuple"
    n=${ADDR[0]}
    m=${ADDR[1]}
    p=${ADDR[2]}
    blockx=32
    blocky=32

    Grid_Dim_x=$((($m + $blockx - 1)/$blockx))
    Grid_Dim_y=$((($n + $blocky - 1)/$blocky))
    echo "Running block BLOCK_DIM_X=$blockx BLOCK_DIM_Y=$blocky"
    echo "With GRID_DIM_X=$Grid_Dim_x GRID_DIM_Y=$Grid_Dim_y"
    time ./bin/task1 $n $p $m $blockx $blocky $Grid_Dim_x $Grid_Dim_y
done