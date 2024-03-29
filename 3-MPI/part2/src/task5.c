#include <math.h>
#include <matmul.h>
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <utilities.h>

int main(int argc, char **argv)
{
    double start_time, end_time, exe_time;
    int sizes[4] = {1000, 2000, 4000, 8000};
    char files[4][75] = {"/gpfs/gibbs/project/cpsc424/shared/assignments/assignment3/data/C-1000.dat",
                         "/gpfs/gibbs/project/cpsc424/shared/assignments/assignment3/data/C-2000.dat",
                         "/gpfs/gibbs/project/cpsc424/shared/assignments/assignment3/data/C-4000.dat",
                         "/gpfs/gibbs/project/cpsc424/shared/assignments/assignment3/data/C-8000.dat"};

    MPI_Init(&argc, &argv);
    MPI_Status status;

    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size); // how many processes are there
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // which process am I

    int next_rank = (rank + 1) % size;        // process this rank will send to
    int prev_rank = (rank - 1 + size) % size; // process this rank will receive from

    if (rank == 0)
    {
        // Print a table heading
        printf("Matrix multiplication times:\n   N      TIME (secs)    F-norm of Error\n -----   -------------  -----------------\n");
    }

    for (int run = 0; run < 4; run++)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0)
        {
            srand(12345);
            start_time = MPI_Wtime();
        }
        // matrices to multiply and the result matrix - only the manager allocates these
        double *A = NULL;
        double *B = NULL;
        double *C = NULL;

        int N = sizes[run];
        long sizeNxN = N * N;

        int block_size = N / size;                      // number of rows in each block
        int size_BLOCKxN = block_size * N;              // total number of doubles in a block
        int size_BLOCKxBLOCK = block_size * block_size; // total number of doubles in the multiplication of blockA and blockB

        // all processes allocate memory of their blocks and secondary B buffer
        double *blockA = (double *)calloc(size_BLOCKxN, sizeof(double));
        double *blockB = (double *)calloc(size_BLOCKxN, sizeof(double));
        double *blockC = (double *)calloc(size_BLOCKxN, sizeof(double));
        // tempB is used for each process to receive data from prior process
        double *tempB = (double *)calloc(size_BLOCKxN, sizeof(double));

        // only the manager initializes A, B, C only
        if (rank == 0)
        {
            A = (double *)calloc(sizeNxN, sizeof(double));
            B = (double *)calloc(sizeNxN, sizeof(double));
            C = (double *)calloc(sizeNxN, sizeof(double));

            // Load A row by row
            for (int i = 0; i < sizeNxN; i++)
            {
                A[i] = ((double)rand() / (double)RAND_MAX);
            }
            // Load B column by column
            for (int i = 0; i < sizeNxN; i++)
            {
                B[i] = ((double)rand() / (double)RAND_MAX);
            }
        }

        // send blocks of A to all other processes in a communicator
        MPI_Scatter(A, size_BLOCKxN, MPI_DOUBLE, blockA, size_BLOCKxN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // send blocks of B to all other processes in a communicator
        MPI_Scatter(B, size_BLOCKxN, MPI_DOUBLE, blockB, size_BLOCKxN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // send blocks of C to all other processes in a communicator
        MPI_Scatter(C, size_BLOCKxN, MPI_DOUBLE, blockC, size_BLOCKxN, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Ring Passing B - MPI processes pass their block columns of B to the next higher-ranked process for processing
        for (int step = 0; step < size; step++)
        {
            // perform multiply of permanant block of A on this process with current block of B on this process
            // result is placed into correct place in blockC
            gemm(block_size, N, blockA, blockB, blockC, ((rank - step + size) % size) * block_size);

            // each process now hands its blockB to the next process in the ring
            if (rank % 2 == 0)
            {
                // to avoid deadlock, even ranks send first and then receive
                MPI_Send(blockB, size_BLOCKxN, MPI_DOUBLE, next_rank, 0, MPI_COMM_WORLD);
                MPI_Recv(tempB, size_BLOCKxN, MPI_DOUBLE, prev_rank, 0, MPI_COMM_WORLD, &status);
            }
            else
            {
                // odd ranks receive first and then send
                MPI_Recv(tempB, size_BLOCKxN, MPI_DOUBLE, prev_rank, 0, MPI_COMM_WORLD, &status);
                MPI_Send(blockB, size_BLOCKxN, MPI_DOUBLE, next_rank, 0, MPI_COMM_WORLD);
            }

            // ensure the received tempB is moved to blockB for the next iterations multiply
            double *swap = blockB;
            blockB = tempB;
            tempB = swap;

            MPI_Barrier(MPI_COMM_WORLD);
        }

        // Gather the p blockC into C from the processes
        MPI_Gather(blockC, size_BLOCKxN, MPI_DOUBLE, C, size_BLOCKxN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        // manager memory cleanup and printouts
        if (rank == 0)
        {
            end_time = MPI_Wtime();
            exe_time = end_time - start_time;
            // compare resulting C to actual
            double *Ctrue = (double *)calloc(sizeNxN, sizeof(double));
            FILE *fptr = fopen(files[run], "r");
            fread(Ctrue, sizeof(double), sizeNxN, fptr);
            fclose(fptr);

            double Fnorm = 0.;
            for (int i = 0; i < sizeNxN; i++)
            {
                Fnorm += (Ctrue[i] - C[i]) * (Ctrue[i] - C[i]);
            }
            Fnorm = sqrt(Fnorm);

            // print a table row
            printf("  %5d    %9.4f  %17.12f\n", N, exe_time, Fnorm);
            write_data_to_file("out/results.csv", "task5", N, size, exe_time, Fnorm);

            // manager frees memory of manager only arrays
            free(A);
            free(B);
            free(C);
            free(Ctrue);
        }

        // all processes free memory of their blocks and secondary buffers
        free(blockA);
        free(blockB);
        free(blockC);
        free(tempB);
    }

    MPI_Finalize();

    return 0;
}