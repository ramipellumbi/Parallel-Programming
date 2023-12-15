#include <immintrin.h>
#include <malloc.h>
#include <math.h>
#include <mkl.h>
#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdalign.h>
#include <stdlib.h>
#include <stdio.h>
#include <utilities.h>
#include <timing.h>

#define KC 240       // KC * COL_BLOCK doubles fit in L1, 32K
#define MC 240       // KC * MC doubles fit in L2, 1024k
#define ROW_BLOCK 10 // operate on 10 rows at a time
#define COL_BLOCK 16 // operate on 16 columns at a time
#define ALIGNMENT 64 // align memory addresses on ALIGNMENT boundary

/**
 * Return next multiple of size from num if num is not a multiple of size
 *
 * Returns num if num is a multiple of size
 */
size_t get_padded_value(size_t num, size_t size)
{
    return ((num + size - 1) / size) * size;
}

double *pad_matrix(const double *original, int original_rows, int original_cols, int padded_rows, int padded_cols)
{
    double *padded = (double *)_mm_malloc(padded_rows * padded_cols * sizeof(double), ALIGNMENT);
    memset(padded, 0., padded_cols * padded_rows * sizeof(double));

    for (int i = 0; i < original_rows; ++i)
    {
        memcpy(&padded[i * padded_cols], &original[i * original_cols], original_cols * sizeof(double));
    }

    return padded;
}

/**
 * Returns the wall clock time elapsed in the serial matrix multiplication between A and B
 * using blocking + AVX instructions
 *
 * Inspired by
 *
 * https://www.cs.utexas.edu/users/flame/pubs/GotoTOMS_final.pdf
 *
 * and mainly
 *
 * https://www.cs.utexas.edu/users/flame/pubs/blis3_ipdps14.pdf
 * https://cs.stanford.edu/people/shadjis/blas.html
 *
 * @param A double[N*P] (row wise storage) -> first P entries are row 1
 * @param B double[P*M] (row wise storage) -> first M entries are row 1
 * @param C double[N*M] (row wise storage) -> first M entries are row 1
 * @param N Number of rows of A
 * @param P Number of columns of A (rows of B)
 * @param M Number or columns of B
 */
double matrix_multiply_blocking(const double *A, const double *B, double *C, size_t N, size_t P, size_t M)
{
    double wc_start, wc_end;
    double cpu_start, cpu_end;
    size_t iA, jB, iC;

    // pad the matrices to handle all matrix sizes
    size_t PADDED_M = get_padded_value(M, COL_BLOCK);
    size_t PADDED_P = get_padded_value(P, KC);
    size_t PADDED_N = get_padded_value(N, MC);

    // expect NC to be a multiple of COL_BLOCK
    printf("\nOld: %d %d %d", N, P, M);
    printf("\nNew: %d %d %d", PADDED_N, PADDED_P, PADDED_M);

    double *padded_A = pad_matrix(A, N, P, PADDED_N, PADDED_P);
    double *padded_B = pad_matrix(B, P, M, PADDED_P, PADDED_M);
    double *padded_C = pad_matrix(C, N, M, PADDED_N, PADDED_M);

    size_t NC = PADDED_M;

    timing(&wc_start, &cpu_start);
#pragma omp parallel default(none) shared(NC, PADDED_M, PADDED_P, PADDED_N, padded_A, padded_B, padded_C)
    {

        // jc = 0,...,M-1 in steps of NC
        for (size_t jc = 0; jc < PADDED_M; jc += NC)
        {
            // pc = 0,...,P-1 in steps of KC
            for (size_t pc = 0; pc < PADDED_P; pc += KC)
            {
                // ic = 0,...,N-1 in steps of MC
                for (size_t ic = 0; ic < PADDED_N; ic += MC)
                {
#pragma omp for
                    // jr = 0,...,NC-1 in steps of COL_BLOCK
                    for (size_t jr = 0; jr < NC; jr += COL_BLOCK)
                    {
                        // ir = 0,...,MC-1 in steps of ROW_BLOCK
                        for (size_t ir = 0; ir < MC; ir += ROW_BLOCK)
                        {
                            // 512 bit wide registers (operate on 8 doubles at once)
                            __m512d mB0, mB1, mA0, mA1;

                            // set 0 to all 8 doubles in the registers
                            __m512d result0_0 = _mm512_set1_pd(0);
                            __m512d result1_0 = _mm512_set1_pd(0);
                            __m512d result2_0 = _mm512_set1_pd(0);
                            __m512d result3_0 = _mm512_set1_pd(0);
                            __m512d result4_0 = _mm512_set1_pd(0);
                            __m512d result5_0 = _mm512_set1_pd(0);
                            __m512d result6_0 = _mm512_set1_pd(0);
                            __m512d result7_0 = _mm512_set1_pd(0);
                            __m512d result8_0 = _mm512_set1_pd(0);
                            __m512d result9_0 = _mm512_set1_pd(0);

                            __m512d result0_1 = _mm512_set1_pd(0);
                            __m512d result1_1 = _mm512_set1_pd(0);
                            __m512d result2_1 = _mm512_set1_pd(0);
                            __m512d result3_1 = _mm512_set1_pd(0);
                            __m512d result4_1 = _mm512_set1_pd(0);
                            __m512d result5_1 = _mm512_set1_pd(0);
                            __m512d result6_1 = _mm512_set1_pd(0);
                            __m512d result7_1 = _mm512_set1_pd(0);
                            __m512d result8_1 = _mm512_set1_pd(0);
                            __m512d result9_1 = _mm512_set1_pd(0);

                            for (size_t k = 0; k < KC; k++)
                            {
                                // load 16 consecutive doubles from k'th row of B (8 to mB0 and 8 to mB1)
                                mB0 = _mm512_load_pd(&padded_B[PADDED_M * (k + pc) + jc + jr]);
                                mB1 = _mm512_load_pd(&padded_B[PADDED_M * (k + pc) + jc + jr + 8]);

                                // Load a single value for the k'th col of A
                                // Note: the addresses below must be aligned on a 64-byte boundary
                                mA0 = _mm512_set1_pd(padded_A[k + pc + (ic + ir) * PADDED_P]);
                                mA1 = _mm512_set1_pd(padded_A[k + pc + (ic + ir + 1) * PADDED_P]);

                                result0_0 = _mm512_fmadd_pd(mB0, mA0, result0_0);
                                result0_1 = _mm512_fmadd_pd(mB1, mA0, result0_1);
                                result1_0 = _mm512_fmadd_pd(mB0, mA1, result1_0);
                                result1_1 = _mm512_fmadd_pd(mB1, mA1, result1_1);

                                mA0 = _mm512_set1_pd(padded_A[k + pc + (ic + ir + 2) * PADDED_P]);
                                mA1 = _mm512_set1_pd(padded_A[k + pc + (ic + ir + 3) * PADDED_P]);

                                result2_0 = _mm512_fmadd_pd(mB0, mA0, result2_0);
                                result2_1 = _mm512_fmadd_pd(mB1, mA0, result2_1);
                                result3_0 = _mm512_fmadd_pd(mB0, mA1, result3_0);
                                result3_1 = _mm512_fmadd_pd(mB1, mA1, result3_1);

                                mA0 = _mm512_set1_pd(padded_A[k + pc + (ic + ir + 4) * PADDED_P]);
                                mA1 = _mm512_set1_pd(padded_A[k + pc + (ic + ir + 5) * PADDED_P]);

                                result4_0 = _mm512_fmadd_pd(mB0, mA0, result4_0);
                                result4_1 = _mm512_fmadd_pd(mB1, mA0, result4_1);
                                result5_0 = _mm512_fmadd_pd(mB0, mA1, result5_0);
                                result5_1 = _mm512_fmadd_pd(mB1, mA1, result5_1);

                                mA0 = _mm512_set1_pd(padded_A[k + pc + (ic + ir + 6) * PADDED_P]);
                                mA1 = _mm512_set1_pd(padded_A[k + pc + (ic + ir + 7) * PADDED_P]);

                                result6_0 = _mm512_fmadd_pd(mB0, mA0, result6_0);
                                result6_1 = _mm512_fmadd_pd(mB1, mA0, result6_1);
                                result7_0 = _mm512_fmadd_pd(mB0, mA1, result7_0);
                                result7_1 = _mm512_fmadd_pd(mB1, mA1, result7_1);

                                mA0 = _mm512_set1_pd(padded_A[k + pc + (ic + ir + 8) * PADDED_P]);
                                mA1 = _mm512_set1_pd(padded_A[k + pc + (ic + ir + 9) * PADDED_P]);

                                result8_0 = _mm512_fmadd_pd(mB0, mA0, result8_0);
                                result8_1 = _mm512_fmadd_pd(mB1, mA0, result8_1);
                                result9_0 = _mm512_fmadd_pd(mB0, mA1, result9_0);
                                result9_1 = _mm512_fmadd_pd(mB1, mA1, result9_1);
                            }

                            *((__m512d *)(&padded_C[(ic + ir + 0) * PADDED_M + jc + jr])) += result0_0;
                            *((__m512d *)(&padded_C[(ic + ir + 0) * PADDED_M + jc + jr + 8])) += result0_1;
                            *((__m512d *)(&padded_C[(ic + ir + 1) * PADDED_M + jc + jr])) += result1_0;
                            *((__m512d *)(&padded_C[(ic + ir + 1) * PADDED_M + jc + jr + 8])) += result1_1;
                            *((__m512d *)(&padded_C[(ic + ir + 2) * PADDED_M + jc + jr])) += result2_0;
                            *((__m512d *)(&padded_C[(ic + ir + 2) * PADDED_M + jc + jr + 8])) += result2_1;
                            *((__m512d *)(&padded_C[(ic + ir + 3) * PADDED_M + jc + jr])) += result3_0;
                            *((__m512d *)(&padded_C[(ic + ir + 3) * PADDED_M + jc + jr + 8])) += result3_1;
                            *((__m512d *)(&padded_C[(ic + ir + 4) * PADDED_M + jc + jr])) += result4_0;
                            *((__m512d *)(&padded_C[(ic + ir + 4) * PADDED_M + jc + jr + 8])) += result4_1;
                            *((__m512d *)(&padded_C[(ic + ir + 5) * PADDED_M + jc + jr])) += result5_0;
                            *((__m512d *)(&padded_C[(ic + ir + 5) * PADDED_M + jc + jr + 8])) += result5_1;
                            *((__m512d *)(&padded_C[(ic + ir + 6) * PADDED_M + jc + jr])) += result6_0;
                            *((__m512d *)(&padded_C[(ic + ir + 6) * PADDED_M + jc + jr + 8])) += result6_1;
                            *((__m512d *)(&padded_C[(ic + ir + 7) * PADDED_M + jc + jr])) += result7_0;
                            *((__m512d *)(&padded_C[(ic + ir + 7) * PADDED_M + jc + jr + 8])) += result7_1;
                            *((__m512d *)(&padded_C[(ic + ir + 8) * PADDED_M + jc + jr])) += result8_0;
                            *((__m512d *)(&padded_C[(ic + ir + 8) * PADDED_M + jc + jr + 8])) += result8_1;
                            *((__m512d *)(&padded_C[(ic + ir + 9) * PADDED_M + jc + jr])) += result9_0;
                            *((__m512d *)(&padded_C[(ic + ir + 9) * PADDED_M + jc + jr + 8])) += result9_1;
                        }
                    }
                }
            }
        }
    }

    free(padded_A);
    free(padded_B);

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < M; j++)
        {
            C[i * M + j] = padded_C[i * PADDED_M + j];
        }
    }
    free(padded_C);

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
        exit(-1);
    }

    // ------------------- data load -----------------
    srand(12345);
    int N = atoi(argv[1]);
    int P = atoi(argv[2]);
    int M = atoi(argv[3]);

    double *A, *B, *C, *BLAS_C;
    double wc_start, wc_end;
    double cpu_start, cpu_end;

    A = (double *)_mm_malloc(N * P * sizeof(double), ALIGNMENT);
    B = (double *)_mm_malloc(P * M * sizeof(double), ALIGNMENT);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < P; j++)
        {
            int idx = i * P + j;
            A[idx] = ((double)rand() / (double)RAND_MAX);
        }
    }

    // load B row-wise
    for (int i = 0; i < P; i++)
    {
        for (int j = 0; j < M; j++)
        {
            int idx = i * M + j;
            B[idx] = ((double)rand() / (double)RAND_MAX);
        }
    }

    BLAS_C = (double *)_mm_malloc(N * M * sizeof(double), ALIGNMENT);
    C = (double *)_mm_malloc(N * M * sizeof(double), ALIGNMENT);
    memset(BLAS_C, 0, N * M * sizeof(double));
    memset(C, 0, N * M * sizeof(double));

    // ------------------- perform BLAS first -----------------
    timing(&wc_start, &cpu_start);
    double alpha = 1.0, beta = 0.0;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                N, M, P, alpha, A, P, B, M, beta, BLAS_C, M);
    timing(&wc_end, &cpu_end);
    double wctime_blas = wc_end - wc_start;

    // ------------------- run custom multiply -----------------
    double wctime = matrix_multiply_blocking(A, B, C, N, P, M);

    // ------------------- check custom multiply creates correct results -----------------
    double error = compute_relative_error_between(A, B, C, BLAS_C, N, P, M);

    // Print a table row
    printf("\n(%d, %d, %d) %9.4f  %f\n", N, P, M, wctime, error);
    write_data_to_file_avx("out/results-avx.csv", "avx-omp-non-divisible", N, P, M, KC, MC, ROW_BLOCK, COL_BLOCK, num_threads, ALIGNMENT, wctime, wctime_blas, error);

    _mm_free(A);
    _mm_free(B);
    _mm_free(C);
    _mm_free(BLAS_C);
}
