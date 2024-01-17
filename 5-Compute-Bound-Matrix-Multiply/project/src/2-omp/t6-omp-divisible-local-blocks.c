#include <mkl.h>
#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdalign.h>
#include <stdlib.h>
#include <stdio.h>
#include <utilities.h>
#include <timing.h>

#define BLOCK_SIZE 32
alignas(64) static double blockA[BLOCK_SIZE][BLOCK_SIZE];
alignas(64) static double blockB[BLOCK_SIZE][BLOCK_SIZE];
alignas(64) static double blockC[BLOCK_SIZE][BLOCK_SIZE];

// ensure there is no false sharing between threads by giving each thread its own local matrix copies
#pragma omp threadprivate(blockA, blockB, blockC)

/**
 * Returns the wall clock time elapsed in the matrix multiplication between A and B using tiled matrix multiplication,
 * OMP, and local block copies
 *
 * @param A double[N*P] (row wise storage) -> first P entries are row 1
 * @param B double[P*M] (row wise storage) -> first M entries are row 1
 * @param C double[N*M] (row wise storage) -> first M entries are row 1
 * @param N Number of rows of A
 * @param P Number of columns of A (rows of B)
 * @param M Number or columns of B
 */
double matrix_multiply_blocking(double *A, double *B, double *C, int N, int P, int M)
{
    double wc_start, wc_end;
    double cpu_start, cpu_end;
    size_t iA, jB, iC;

    timing(&wc_start, &cpu_start);
#pragma omp parallel default(none) shared(N, P, M, A, B, C) private(iA, jB, iC)
    {
// Compute A1 * B1
#pragma omp for schedule(runtime)
        for (size_t ii = 0; ii < N; ii += BLOCK_SIZE)
        {
            for (size_t jj = 0; jj < M; jj += BLOCK_SIZE)
            {
                // clear blockC
                for (size_t i = 0; i < BLOCK_SIZE; i++)
                {
                    for (size_t j = 0; j < BLOCK_SIZE; j++)
                    {
                        blockC[i][j] = 0.;
                    }
                }

                for (size_t kk = 0; kk < P; kk += BLOCK_SIZE)
                {
                    // Copy blocks of A and B into blockA and blockB
                    for (size_t i = 0; i < BLOCK_SIZE; i++)
                    {
                        iA = (ii + i) * P + kk;
                        jB = (kk + i) * M + jj;
                        for (size_t j = 0; j < BLOCK_SIZE; j++)
                        {
                            blockA[i][j] = A[iA + j];
                            blockB[i][j] = B[jB + j];
                        }
                    }

                    // Perform block multiplication
                    for (size_t k = 0; k < BLOCK_SIZE; k++)
                    {
                        for (size_t i = 0; i < BLOCK_SIZE; i++)
                        {
                            double r = blockA[i][k];
                            for (size_t j = 0; j < BLOCK_SIZE; j++)
                            {
                                blockC[i][j] += r * blockB[k][j];
                            }
                        }
                    }
                }

                // copy local block C back into C
                for (size_t i = 0; i < BLOCK_SIZE; i++)
                {
                    iC = (ii + i) * M + jj;
                    for (size_t j = 0; j < BLOCK_SIZE; j++)
                    {
                        C[iC + j] -= blockC[i][j];
                    }
                }
            }
        }
    }

    timing(&wc_end, &cpu_end);
    double elapsed_time = wc_end - wc_start;

    return elapsed_time;
}

int main(int argc, char **argv)
{
    if (argc != 4)
    {
        printf("Usage: a-serial <N> <P> <M>\n");
        exit(-1);
    }
    int num_threads = get_environment_value("OMP_NUM_THREADS");
    if (num_threads == -1)
    {
        fprintf(stderr, "OMP_NUM_THREADS not set");
        num_threads = 1;
    }

    // ------------------- data load -----------------
    srand(12345);
    int N = atoi(argv[1]);
    int P = atoi(argv[2]);
    int M = atoi(argv[3]);

    if (N % BLOCK_SIZE != 0 || P % BLOCK_SIZE != 0 || M % BLOCK_SIZE != 0)
    {
        fprintf(stderr, "Matrix dimensions must be a multiple of %d", BLOCK_SIZE);
        exit(-1);
    }

    double *A, *B, *C;
    double wc_start, wc_end;
    double cpu_start, cpu_end;

    load_random_matrices(&A, &B, N, P, M);
    C = (double *)calloc(N * M, sizeof(double));

    // ------------------- perform BLAS first -----------------
    timing(&wc_start, &cpu_start);
    double alpha = 1.0, beta = 0.0;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                N, M, P, alpha, A, P, B, M, beta, C, M);
    timing(&wc_end, &cpu_end);
    double wctime_blas = wc_end - wc_start;

    // ------------------- run custom multiply -----------------
    double wctime = matrix_multiply_blocking(A, B, C, N, P, M);

    // ------------------- check custom multiply creates correct results -----------------
    double error = compute_relative_error(A, B, C, N, P, M);

    // Print a table row
    printf("\n(%d, %d, %d) %9.4f  %f\n", N, P, M, wctime, error);
    write_data_to_file("out/results-omp.csv", "omp-local", N, P, M, BLOCK_SIZE, num_threads, wctime, wctime_blas, error);

    free(A);
    free(B);
    free(C);
}
