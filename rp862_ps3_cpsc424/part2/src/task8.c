#include <math.h>
#include <matmul.h>
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <utilities.h>

void swap(double **arr_1, double **arr_2)
{
    double *temp = *arr_1;
    *arr_1 = *arr_2;
    *arr_2 = temp;
}

int main(int argc, char **argv)
{
    double start_time, end_time, exe_time;
    const char *file = "/gpfs/gibbs/project/cpsc424/shared/assignments/assignment3/data/C-7633.dat";
    int N = 7633;

    MPI_Init(&argc, &argv);
    MPI_Status status;

    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size); // how many processes are there
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // which process am I

    int next_rank = (rank + 1) % size;        // process this rank will send to
    int prev_rank = (rank - 1 + size) % size; // process this rank will receive from

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

    long sizeNxN = N * N;

    /**
     * If N % p == 0, then N / p evenly distributes the blocks
     *
     * When N % p != 0, then the remainder will be less than p and
     * the remaining rows are placed in the first N % p ranks
     *
     * This ensures that among the p blocks, all rows are accounted for
     */
    int base = N / size;
    int remainder = N % size;

    /**
     * Need to keep track of the displacements from block to block
     * since some will be (base + 1 x N) and some will be (base x N)
     */
    int *block_sizes = (int *)calloc(size, sizeof(int));
    int *sizes_BLOCKxN = (int *)calloc(size, sizeof(int));
    int *displacements = (int *)calloc(size, sizeof(int));

    int sum = 0;
    int max_block_size = base + 1;
    int max_size_BLOCKxN = (base + 1) * N;
    for (int i = 0; i < size; i++)
    {
        if (i < remainder)
        {
            block_sizes[i] = base + 1;
            sizes_BLOCKxN[i] = (base + 1) * N;
        }
        else
        {
            block_sizes[i] = base;
            sizes_BLOCKxN[i] = base * N;
        }
        displacements[i] = sum;
        sum += sizes_BLOCKxN[i];
    }

    // initialize A, B, C to receive the appropriate initial block sizes
    double *blockA = (double *)calloc(sizes_BLOCKxN[rank], sizeof(double));
    // every blockB is allocated to max possible size to make receiving from ranks more simple
    double *blockB = (double *)calloc(max_size_BLOCKxN, sizeof(double));
    double *blockC = (double *)calloc(sizes_BLOCKxN[rank], sizeof(double));
    // initialized to handle receiving the prior processes blockB - initialized to max possible size to make swapping simple
    double *tempB = (double *)calloc(max_size_BLOCKxN, sizeof(double));

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
    MPI_Scatterv(A, sizes_BLOCKxN, displacements, MPI_DOUBLE, blockA, sizes_BLOCKxN[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // send blocks of B to all other processes in a communicator
    MPI_Scatterv(B, sizes_BLOCKxN, displacements, MPI_DOUBLE, blockB, sizes_BLOCKxN[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // send blocks of C to all other processes in a communicator
    MPI_Scatterv(C, sizes_BLOCKxN, displacements, MPI_DOUBLE, blockC, sizes_BLOCKxN[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Ring Passing B - MPI processes pass their block columns of B to the next higher-ranked process for processing
    for (int step = 0; step < size; step++)
    {
        /**
         * Process i and process j may have different memory allocations for their blocks of B
         * Thus, blockB and tempB need to be reallocated within the current rank to
         * prepare for the next receive
         *
         * At each step we receive the block of B initially allocated to (prev_rank - step + size) % size
         */
        int rank_receiving_from = (prev_rank - step + size) % size; // at step i, this is i^th left neighbors left neighbor's assignment
        int rank_sending = (rank - step + size) % size;             // at step i, this is the i^th left neighbors initial assignment

        int receiving_BLOCKxN = sizes_BLOCKxN[rank_receiving_from]; // this is the receiving buffer size
        int sending_BLOCKxN = sizes_BLOCKxN[rank_sending];          // this is the buffer size this rank is sending

        int previous_block_size = block_sizes[rank_receiving_from]; // this is the block_size from the receiving rank
        int current_block = block_sizes[rank_sending];              // this is the block_size this rank is currently handling

        // The computed col_offset is where we offset to for placing the result in blockC to account for variable block lengths
        int cumulative_col_size = 0;
        for (int i = 0; i < rank_sending; ++i)
        {
            cumulative_col_size += block_sizes[i];
        }
        int col_offset = cumulative_col_size % N;

        // perform multiply of the static block of A on this process with current block of B on this process
        gemm_k(block_sizes[rank], N, current_block, blockA, blockB, blockC, col_offset);

        // each process now hands its blockB to the next process in the ring
        if (rank % 2 == 0)
        {
            // to avoid deadlock, even ranks send first and then receive
            MPI_Send(blockB, sending_BLOCKxN, MPI_DOUBLE, next_rank, 0, MPI_COMM_WORLD);
            MPI_Recv(tempB, receiving_BLOCKxN, MPI_DOUBLE, prev_rank, 0, MPI_COMM_WORLD, &status);
        }
        else
        {
            // odd ranks receive first and then send
            MPI_Recv(tempB, receiving_BLOCKxN, MPI_DOUBLE, prev_rank, 0, MPI_COMM_WORLD, &status);
            MPI_Send(blockB, sending_BLOCKxN, MPI_DOUBLE, next_rank, 0, MPI_COMM_WORLD);
        }

        double *swap = blockB;
        blockB = tempB;
        tempB = swap;

        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Gather the p blockC into C from the processes, accounting for the variable sizing
    MPI_Gatherv(blockC, sizes_BLOCKxN[rank], MPI_DOUBLE, C, sizes_BLOCKxN, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    // manager memory cleanup and printouts
    if (rank == 0)
    {
        end_time = MPI_Wtime();
        exe_time = end_time - start_time;
        // compare resulting C to actual
        double *Ctrue = (double *)calloc(sizeNxN, sizeof(double));
        FILE *fptr = fopen(file, "r");
        fread(Ctrue, sizeof(double), sizeNxN, fptr);
        fclose(fptr);

        // Print a table heading
        printf("Matrix multiplication times:\n   N      TIME (secs)    F-norm of Error\n -----   -------------  -----------------\n");

        double Fnorm = 0.;
        for (int i = 0; i < sizeNxN; i++)
        {
            Fnorm += (Ctrue[i] - C[i]) * (Ctrue[i] - C[i]);
        }
        Fnorm = sqrt(Fnorm);

        // print a table row
        printf("  %5d    %9.4f  %17.12f\n", N, exe_time, Fnorm);
        write_data_to_file("out/results.csv", "task8", N, size, exe_time, Fnorm);

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

    MPI_Finalize();

    return 0;
}
