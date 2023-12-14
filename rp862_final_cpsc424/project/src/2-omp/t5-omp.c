#include <mkl.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <timing.h>
#include <utilities.h>

#define BLOCK_SIZE 16

/**
 * Returns the wall clock time elapsed in the matrix multiplication between A and B using tiled matrix multiplication with OMP
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
    /**
     * Can think of A * B as multiplication of block matrices
     *
     *     A1, A2
     * A = A3, A4
     *
     * B = B1, B2
     *     B3, B4
     *
     * Yielding
     *
     * C = A1*B1 + A2*B3, A1*B2 + A2*B4
     *     A3*B1 + A4*BB, A3*B2 + A4*B4
     *
     * The matrices A and B are divided as above so that A1 and B1
     * make use of efficient cache utilization by multiplying appropriate
     * sub blocks. The remaining computation is done in cleanup loops.
     */
    double wc_start, wc_end;
    double cpu_start, cpu_end;
    size_t iA, jB, iC;

    timing(&wc_start, &cpu_start);
    // transpose B
    double *B_TR = (double *)malloc(P * M * sizeof(double));
    for (size_t i = 0; i < P; i++)
    {
        for (size_t j = 0; j < M; j++)
        {
            B_TR[j * P + i] = B[i * M + j];
        }
    }

    // Calculate the bounds for the blocked multiplication
    size_t N_block_max = (N / BLOCK_SIZE) * BLOCK_SIZE;
    size_t P_block_max = (P / BLOCK_SIZE) * BLOCK_SIZE;
    size_t M_block_max = (M / BLOCK_SIZE) * BLOCK_SIZE;

#pragma omp parallel default(none) shared(N_block_max, P_block_max, M_block_max, A, B_TR, C, M, P, N) private(iA, jB, iC)
    {
// Compute A1 * B1
#pragma omp for schedule(runtime)
        for (size_t ii = 0; ii < N_block_max; ii += BLOCK_SIZE)
        {
            for (size_t jj = 0; jj < M_block_max; jj += BLOCK_SIZE)
            {
                for (size_t kk = 0; kk < P_block_max; kk += BLOCK_SIZE)
                {
                    for (size_t i = ii; i < ii + BLOCK_SIZE; i++)
                    {
                        iA = i * P;
                        for (size_t j = jj; j < jj + BLOCK_SIZE; j++)
                        {
                            iC = i * M + j;
                            jB = j * P;
                            double sum = 0.;
#pragma omp simd reduction(+ : sum)
                            for (size_t k = kk; k < kk + BLOCK_SIZE; k++)
                            {
                                sum += A[iA + k] * B_TR[jB + k];
                            }
                            C[iC] -= sum;
                        }
                    }
                }
            }
        }

// This for loop does A2*B3 and A2*B4 and places into appropriate block of C
#pragma omp for schedule(runtime)
        for (size_t i = 0; i < N_block_max; i++)
        {
            size_t iA = i * P;
            for (size_t j = 0; j < M; j++)
            {
                jB = j * P;
                iC = i * M + j;
                double sum = 0.;
#pragma omp simd reduction(+ : sum)
                for (size_t k = P_block_max; k < P; k++)
                {
                    sum += A[iA + k] * B_TR[jB + k];
                }
                C[iC] -= sum;
            }
        }

// This for loop does [A3 A4] * B
#pragma omp for schedule(runtime)
        for (size_t i = N_block_max; i < N; i++)
        {
            iA = i * P;
            for (size_t j = 0; j < M; j++)
            {
                jB = j * P;
                iC = i * M + j;
                double sum = 0.;
#pragma omp simd reduction(+ : sum)
                for (size_t k = 0; k < P; k++)
                {
                    sum += A[iA + k] * B_TR[jB + k];
                }
                C[iC] -= sum;
            }
        }

// This for loop does A1*B2
#pragma omp for schedule(runtime)
        for (size_t i = 0; i < N_block_max; i++)
        {
            size_t iA = i * P;
            for (size_t j = M_block_max; j < M; j++)
            {
                jB = j * P;
                iC = i * M + j;
                double sum = 0.;
#pragma omp simd reduction(+ : sum)
                for (size_t k = 0; k < P_block_max; k++)
                {
                    sum += A[iA + k] * B_TR[jB + k];
                }
                C[iC] -= sum;
            }
        }
    }

    timing(&wc_end, &cpu_end);
    double elapsed_time = wc_end - wc_start;

    free(B_TR);

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
    write_data_to_file("out/results-omp.csv", "omp", N, P, M, BLOCK_SIZE, num_threads, wctime, wctime_blas, error);

    free(A);
    free(B);
    free(C);
}