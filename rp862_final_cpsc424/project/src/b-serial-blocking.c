#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <timing.h>
#include <utilities.h>

#define BLOCK_SIZE 32

/**
 * Returns the wall clock time elapsed in the naive matrix multiplication between A and B
 *
 * @param A double[M*N] (row wise storage) -> first N entries are row 1
 * @param B double[N*K] (column wise storage) -> first N entries are column 1
 * @param C double[M*k] (row wise storage) -> first N entries are row 1
 * @param M Number of rows of A
 * @param N Number of columns of A (rows of B)
 * @param K Number or columns of B
 */
double matrix_multiply_blocking(double *A, double *B, double *C, int M, int N, int K)
{
    double wc_start, wc_end;
    double cpu_start, cpu_end;
    int iA, jB, iC;

    timing(&wc_start, &cpu_start);
    // Calculate the bounds for the blocked multiplication
    int M_block_max = (M / BLOCK_SIZE) * BLOCK_SIZE;
    int K_block_max = (K / BLOCK_SIZE) * BLOCK_SIZE;
    int N_block_max = (N / BLOCK_SIZE) * BLOCK_SIZE;

    // Blocked multiplication
    for (int ii = 0; ii < M_block_max; ii += BLOCK_SIZE)
    {
        for (int jj = 0; jj < K_block_max; jj += BLOCK_SIZE)
        {
            for (int kk = 0; kk < N_block_max; kk += BLOCK_SIZE)
            {
                for (int i = ii; i < ii + BLOCK_SIZE; i++)
                {
                    iA = i * N;
                    for (int j = jj; j < jj + BLOCK_SIZE; j++)
                    {
                        jB = j * N;
                        iC = i * K + j;
                        for (int k = kk; k < kk + BLOCK_SIZE; k++)
                        {
                            C[iC] += A[iA + k] * B[jB + k];
                        }
                    }
                }
            }
        }
    }

    // A is M x N and we only processed a sub-matrix M1 x N1
    // B is N x K and we only processed a sub-matrix N1 x K1
    for (int i = 0; i < M_block_max; i++)
    {
        iA = i * N;
        for (int j = 0; j < K; j++)
        {
            jB = j * N;
            iC = i * K + j;
            for (int k = N_block_max; k < N; k++)
            {
                C[iC] += A[iA + k] * B[jB + k];
            }
        }
    }

    for (int i = M_block_max; i < M; i++)
    {
        iA = i * N;
        for (int j = 0; j < K; j++)
        {
            jB = j * N;
            iC = i * K + j;
            for (int k = 0; k < N; k++)
            {
                C[iC] += A[iA + k] * B[jB + k];
            }
        }
    }

    for (int i = 0; i < M_block_max; i++)
    {
        iA = i * N;
        for (int j = K_block_max; j < K; j++)
        {
            jB = j * N;
            iC = i * K + j;
            for (int k = 0; k < N_block_max; k++)
            {
                C[iC] += A[iA + k] * B[jB + k];
            }
        }
    }

    timing(&wc_end, &cpu_end);
    double elapsed_time = wc_end - wc_start;

    return elapsed_time;
}

int main(int argc, char **argv)
{
    double *A, *B, *C;
    size_t sizes[5] = {1000, 2000, 4000, 7633, 8000};
    char files[5][75] = {"/gpfs/gibbs/project/cpsc424/shared/assignments/assignment3/data/C-1000.dat",
                         "/gpfs/gibbs/project/cpsc424/shared/assignments/assignment3/data/C-2000.dat",
                         "/gpfs/gibbs/project/cpsc424/shared/assignments/assignment3/data/C-4000.dat",
                         "/gpfs/gibbs/project/cpsc424/shared/assignments/assignment3/data/C-7633.dat",
                         "/gpfs/gibbs/project/cpsc424/shared/assignments/assignment3/data/C-8000.dat"};

    // Print a table heading
    printf("Matrix multiplication times:\n   N      TIME (secs)    F-norm of Error\n -----   -------------  -----------------\n");

    // Now run the four test cases
    for (int run = 0; run < 5; run++)
    {
        srand(12345);
        size_t size = sizes[run];
        int size_A = size * size;
        int size_B = size * size;
        int size_C = size * size;

        A = (double *)calloc(size_A, sizeof(double));
        B = (double *)calloc(size_B, sizeof(double));
        C = (double *)calloc(size_C, sizeof(double));

        for (int i = 0; i < size_A; i++)
        {
            A[i] = ((double)rand() / (double)RAND_MAX);
        }
        for (int i = 0; i < size_B; i++)
        {
            B[i] = ((double)rand() / (double)RAND_MAX);
        }

        double wctime = matrix_multiply_blocking(A, B, C, size, size, size);

        free(A);
        free(B);

        double Fnorm = compute_fnorm(files[run], C, size_C);

        // Print a table row
        printf("  %5d    %9.4f  %17.12f\n", size, wctime, Fnorm);
        write_data_to_file("out/results.csv", "b-serial-blocking", size, BLOCK_SIZE, 1, wctime, Fnorm);

        free(C);
    }
}
