# Parallel Programming

## Introduction

This repository contains code done for assignments as part of Yale University's
CPSC 524: Parallel Programming course. The course was taught by the excellent
Andrew Sherman.

## Assignments

### Assignment 1: [Divide and Vector Triad Performance](1-CPU-Profiling/report.pdf)

This assignment investigated the characteristics of the processors on the Grace Linux Cluster. First, a program is written to approximate π by numerically integrating the function

`f(x) = 4.0 / (1.0 + x*x)`

from 0 to 1. The divide latency is estimated by assuming that floating point division cannot be pipelined and
assuming that the calculation is dominated by the cost of the divide. Next, a program is written to perform
the vector triad benchmark, an operation that evaluates the memory bandwidth and computational speed
for vector operations. The MFlops are plotted vs. log N , where N is the size of the arrays in the vector triad
benchmark. All benchmarks are run on one core of a compute node.

### Assignment 2: OpenMP - [Mandelbrot Set](2-OpenMP/docs/report.pdf)

This assignment focused on assessing the performance gains achievable through
the application of loop directives and tasks in parallelizing the computation of the Mandelbrot set,
detailed as follows:

1. Define a rectangular region in the complex plane bounded by lower left corner `-2.0 + 0.0i` and
   upper right corner `0.5 + 1.25i`.
2. Discretize the rectangular region into cells with side length of `0.001`.
3. For each cell in the discretized rectangle, generate a random point `c` in the cell.
4. Compute the Mandelbrot update `z <- z*z + c` for up to 25,000 iterations.
   - If `|z * z| > 2` is met at any iteration, the point `c` is marked outside the set.
   - Conversely, if all iterations complete without this condition being met, the point is
     considered inside the Mandelbrot set.
5. Compute the Mandelbrot area as

```
A = 2.0 * (#In / (#Out + #In)) * Area(rectangular region)
```

OpenMP was used to parallelize the computation of the Mandelbrot set. Further, I
implemented AVX-512 instructions to speed up the computation - beating the performance of
the professor's solution by ~3x.

Additionally, the assignment showed the shortcoming of naive random number generation in parallel.
The linked report details the implementation of the Leapfrog method to generate random numbers in parallel.

### Assignment 3: MPI - [Matrix Multiplication](3-MPI/docs/report.pdf)

This assignment investigates the implementation and performance of an MPI- based program tasked with
parallel computation of large double precision matrices. The program’s objective is to effectively utilize multiple processors to perform the multiplication, optimizing the use of system resources and minimizing computation time. The implementation required careful consideration of data distribution and process synchronization, ensuring accuracy while maximizing computational efficiency. Performance analysis focuses on the scalability of the solution across various matrix sizes and the number of processors used. The linked report outlines the design considerations, describes the implementation strategy, and presents an evaluation of the program’s computational performance.

### Assignment 4: CUDA - [Matrix Multiplication](4-Cuda-GPU-Programming/docs//report.pdf)

This assignment explores GPU programming using CUDA, centered around the task of matrix multiplication.
The primary challenge involves constructing a CUDA kernel capable of multiplying two random rectangular
matrices. This initial task serves as a gateway into the realms of parallel computing and efficient memory
management, foundational elements in leveraging GPU architecture for computational tasks.
The report addresses more advanced techniques, such as the utilization of shared memory and the strategic
optimization of thread computations, essential for enhancing the computational performance.

### Final Project: [Compute Bound CPU Matrix Multiplication](5-Compute-Bound-Matrix-Multiply)

#### Abstract

Inspired by a YouTube video [Adding Nested Loops Makes this Algorithm 120x FASTER](https://www.youtube.com/watch?v=QGYvbsHDPxo&ab_channel=DepthBuffer), this project aimed
to implement a highly optimized matrix multiplication algorithm on CPUs. The video demonstrates how
the performance of a naive matrix multiplication algorithm can be improved by simply reordering the loops
and implementing a blocking strategy. The video also highlights the importance of memory access patterns
in matrix multiplication.
Matrix multiplication is a cornerstone operation in numerous computational fields, ranging from scientific
computing to machine learning. At its core, the operation involves the element-wise multiplication and
summation of elements across two matrices to produce a third matrix. The theoretical simplicity of this
operation belies its computational complexity, particularly when dealing with large matrices.

Matrix multiplication scales with the size of the input matrices, often resulting in a significant computational challenge
for even modern day processors. This challenge is accentuated by the fact that CPUs, with their limited
number of cores and sequential processing capabilities, are often outperformed by GPUs in parallelizable
tasks like matrix multiplication. However, understanding and optimizing matrix multiplication on CPUs is
crucial, as CPUs are more universally accessible and are often the primary computing resource available in
many environments.

The difficulty in optimizing matrix multiplication on CPUs stems from several factors. First, the compu-
tational intensity: as the size of the matrices increases, the number of calculations grows cubically, leading
to a steep increase in the required computational resources. Second, memory access patterns play a critical
role: as the size of the matrices increase, the total memory accesses increase quadratically. Efficient matrix
multiplication algorithms must minimize cache misses and effectively utilize the CPU cache hierarchy. This
is challenging due to the non-contiguous memory access patterns inherent in naive matrix multiplication.

The current state of the art in matrix multiplication optimization are built on top of Basic Linear Algebra
Subprograms (BLAS). The magic of BLAS lies in its ability to significantly optimize these computationally
intensive operations. These routines are meticulously engineered to exploit the underlying architecture of
CPUs to their fullest, leveraging techniques such as loop unrolling, blocking for cache, and efficient use of
SIMD instructions. These optimizations allow BLAS to achieve performance levels that are often an order
of magnitude faster than naive implementations. This project started as ‘investigate this videos claims’ and
ended with the goal of matching the performance of Intel’s Math Kernel Library (MKL) implementation of
dgemm: a highly optimized matrix multiplication routine for double precision matrices.
