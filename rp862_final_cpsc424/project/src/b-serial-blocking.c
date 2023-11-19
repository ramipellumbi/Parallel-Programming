#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <timing.h>
#include <utilities.h>

static const int BLOCK_SIZE = 8;

/**
 * Returns the wall clock time elapsed in the naive matrix multiplication between A and B
 *
 * @param A double[M*N] (row wise storage)
 * @param B double[N*K] (column wise storage)
 * @param C double[L] (row wise storage) where L >= M*K.
 * @param M Number of rows of A
 * @param N Number of columns of A (rows of B)
 * @param K Number or columns of B
 */
double matrix_multiply_naive(double *A, double *B, double *C, int M, int N, int K)
{
    double wc_start, wc_end;
    double cpu_start, cpu_end;

    timing(&wc_start, &cpu_start);
    for (int ii = 0; ii < N; ii += BLOCK_SIZE)
    {
        for (int jj = 0; jj < N; jj += BLOCK_SIZE)
        {
            for (int kk = 0; kk < N; kk += BLOCK_SIZE)
            {
                for (int i = ii; i < ii + BLOCK_SIZE; i++)
                {
                    int iA = i * N;
                    for (int j = jj; j < jj + BLOCK_SIZE; j++)
                    {
                        int jB = j * N;
                        int iC = i * N + j;
                        for (int k = kk; k < kk + BLOCK_SIZE; k++)
                        {
                            C[iC] += A[iA + k] * B[jB + k];
                        }
                    }
                }
            }
        }
    }
    timing(&wc_end, &cpu_end);
    double elapsed_time = wc_end - wc_start;

    return elapsed_time;
}

typedef struct Size
{
    int M;
    int N;
    int K;
} Size;

int main(int argc, char **argv)
{
    double *A, *B, *C;
    Size sizes[4] = {{
                         M : 1000,
                         N : 1000,
                         K : 1000,
                     },
                     {
                         M : 2000,
                         N : 2000,
                         K : 2000,
                     },
                     {
                         M : 4000,
                         N : 4000,
                         K : 4000,
                     },
                     {
                         M : 8000,
                         N : 8000,
                         K : 8000,
                     }};
    char files[4][75] = {"/gpfs/gibbs/project/cpsc424/shared/assignments/assignment3/data/C-1000.dat",
                         "/gpfs/gibbs/project/cpsc424/shared/assignments/assignment3/data/C-2000.dat",
                         "/gpfs/gibbs/project/cpsc424/shared/assignments/assignment3/data/C-4000.dat",
                         "/gpfs/gibbs/project/cpsc424/shared/assignments/assignment3/data/C-8000.dat"};

    // Print a table heading
    printf("Matrix multiplication times:\n   N      TIME (secs)    F-norm of Error\n -----   -------------  -----------------\n");

    // Now run the four test cases
    for (int run = 0; run < 4; run++)
    {
        srand(12345);
        Size size = sizes[run];
        int M = size.M;
        int N = size.N;
        int K = size.K;

        int size_A = M * N;
        int size_B = N * K;
        int size_C = M * K;

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

        double wctime = matrix_multiply_naive(A, B, C, M, N, K);

        free(A);
        free(B);

        double Fnorm = compute_fnorm(files[run], C, size_C);

        // Print a table row
        printf("  %5d    %9.4f  %17.12f\n", N, wctime, Fnorm);
        write_data_to_file("out/results.csv", "b-serial-blocking", N, 1, wctime, Fnorm);

        free(C);
    }
}
