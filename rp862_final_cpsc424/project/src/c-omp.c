#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdalign.h>
#include <stdlib.h>
#include <stdio.h>
#include <utilities.h>
#include <timing.h>

static const int BLOCK_SIZE = 16;
// alignas(512) static double blockA[BLOCK_SIZE * BLOCK_SIZE];
// alignas(512) static double blockB[BLOCK_SIZE * BLOCK_SIZE];
// alignas(512) static double blockC[BLOCK_SIZE * BLOCK_SIZE];
// #pragma omp threadprivate(blockA, blockB, blockC)

/**
 * Returns the wall clock time elapsed in the matrix multiplication between A and B
 *
 * @param A double[M*N] (row wise storage)
 * @param B double[N*K] (column wise storage)
 * @param C double[L] (row wise storage) where L >= M*K.
 * @param M Number of rows of A
 * @param N Number of columns of A (rows of B)
 * @param K Number or columns of B
 */
double matrix_multiply_omp(double *A, double *B, double *C, int M, int N, int K)
{
    double wc_start, wc_end;
    double cpu_start, cpu_end;
    int blockNum = M / BLOCK_SIZE;

    timing(&wc_start, &cpu_start);
#pragma omp parallel default(none) shared(A, B, C, blockNum)
    {
#pragma omp for
        for (size_t bi = 0; bi < blockNum; bi++)
        {
            for (size_t bj = 0; bj < blockNum; bj++)
            {
                for (size_t bk = 0; bk < blockNum; bk++)
                {
                    for (int i = 0; i < BLOCK_SIZE; i++)
                    {
                        for (int j = 0; j < BLOCK_SIZE; j++)
                        {
                            int iC = bi * BLOCK_SIZE * blockNum * BLOCK_SIZE + i * blockNum * BLOCK_SIZE + bj * BLOCK_SIZE + j;
                            double dot_product = 0.;
                            C[iC] = 0;
#pragma omp simd reduction(+ : dot_product)
                            for (int k = 0; k < BLOCK_SIZE; k++)
                            {
                                size_t iA = bi * BLOCK_SIZE * blockNum * BLOCK_SIZE + i * blockNum * BLOCK_SIZE + bk * BLOCK_SIZE + k;
                                size_t jB = bj * BLOCK_SIZE * blockNum * BLOCK_SIZE + k * blockNum * BLOCK_SIZE + bj * BLOCK_SIZE + j;
                                dot_product += A[iA] * B[jB];
                            }
                            C[iC] += dot_product;
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
    printf("\n Matrix multiplication times:\n   N      TIME (secs)    F-norm of Error\n -----   -------------  -----------------\n");

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

        // Allocate size bytes of memory, aligned to the alignment specified in align, and return a pointer to the allocated memory.
        double *A = (double *)malloc(M * K * sizeof(double));
        double *B = (double *)malloc(K * N * sizeof(double));
        C = (double *)calloc(size_C, sizeof(double));

        for (int i = 0; i < size_A; i++)
        {
            A[i] = ((double)rand() / (double)RAND_MAX);
        }
        for (int i = 0; i < size_B; i++)
        {
            B[i] = ((double)rand() / (double)RAND_MAX);
        }

        double wctime = matrix_multiply_omp(A, B, C, M, N, K);

        free(A);
        free(B);

        double Fnorm = compute_fnorm(files[run], C, size_C);

        // Print a table row
        printf("  %5d    %9.4f  %17.12f\n", N, wctime, Fnorm);
        write_data_to_file("out/results.csv", "c-omp", N, 1, wctime, Fnorm);

        free(C);
    }
}
