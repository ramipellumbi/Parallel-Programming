#include <immintrin.h>
#include <malloc.h>
#include <math.h>
#include <mkl.h>
#include <stdio.h>
#include <string.h>
#include <stdalign.h>
#include <stdlib.h>
#include <stdio.h>
#include <utilities.h>
#include <timing.h>

#define KC 240       // KC * COL_BLOCK doubles fit in L1, 32K
#define MC 240       // KC * MC doubles fit in L2, 1024k
#define ROW_BLOCK 10 // operate on 6 rows at a time
#define COL_BLOCK 16 // operate on 16 columns at a time
#define ALIGNMENT 64 // align memory addresses on ALIGNMENT boundary

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

    size_t NC = M;

    timing(&wc_start, &cpu_start);
#pragma omp parallel
    {

        // jc = 0,...,M-1 in steps of NC
        for (size_t jc = 0; jc < M; jc += NC)
        {
            // pc = 0,...,P-1 in steps of KC
            for (size_t pc = 0; pc < P; pc += KC)
            {
                // ic = 0,...,N-1 in steps of MC
                for (size_t ic = 0; ic < N; ic += MC)
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
                                mB0 = _mm512_load_pd(&B[M * (k + pc) + jc + jr]);
                                mB1 = _mm512_load_pd(&B[M * (k + pc) + jc + jr + 8]);

                                // Load a single value for the k'th col of A
                                // Note: the addresses below must be aligned on a 64-byte boundary
                                mA0 = _mm512_set1_pd(A[k + pc + (ic + ir) * P]);
                                mA1 = _mm512_set1_pd(A[k + pc + (ic + ir + 1) * P]);

                                result0_0 = _mm512_fmadd_pd(mB0, mA0, result0_0);
                                result0_1 = _mm512_fmadd_pd(mB1, mA0, result0_1);
                                result1_0 = _mm512_fmadd_pd(mB0, mA1, result1_0);
                                result1_1 = _mm512_fmadd_pd(mB1, mA1, result1_1);

                                mA0 = _mm512_set1_pd(A[k + pc + (ic + ir + 2) * P]);
                                mA1 = _mm512_set1_pd(A[k + pc + (ic + ir + 3) * P]);

                                result2_0 = _mm512_fmadd_pd(mB0, mA0, result2_0);
                                result2_1 = _mm512_fmadd_pd(mB1, mA0, result2_1);
                                result3_0 = _mm512_fmadd_pd(mB0, mA1, result3_0);
                                result3_1 = _mm512_fmadd_pd(mB1, mA1, result3_1);

                                mA0 = _mm512_set1_pd(A[k + pc + (ic + ir + 4) * P]);
                                mA1 = _mm512_set1_pd(A[k + pc + (ic + ir + 5) * P]);

                                result4_0 = _mm512_fmadd_pd(mB0, mA0, result4_0);
                                result4_1 = _mm512_fmadd_pd(mB1, mA0, result4_1);
                                result5_0 = _mm512_fmadd_pd(mB0, mA1, result5_0);
                                result5_1 = _mm512_fmadd_pd(mB1, mA1, result5_1);

                                mA0 = _mm512_set1_pd(A[k + pc + (ic + ir + 6) * P]);
                                mA1 = _mm512_set1_pd(A[k + pc + (ic + ir + 7) * P]);

                                result6_0 = _mm512_fmadd_pd(mB0, mA0, result6_0);
                                result6_1 = _mm512_fmadd_pd(mB1, mA0, result6_1);
                                result7_0 = _mm512_fmadd_pd(mB0, mA1, result7_0);
                                result7_1 = _mm512_fmadd_pd(mB1, mA1, result7_1);

                                mA0 = _mm512_set1_pd(A[k + pc + (ic + ir + 8) * P]);
                                mA1 = _mm512_set1_pd(A[k + pc + (ic + ir + 9) * P]);

                                result8_0 = _mm512_fmadd_pd(mB0, mA0, result8_0);
                                result8_1 = _mm512_fmadd_pd(mB1, mA0, result8_1);
                                result9_0 = _mm512_fmadd_pd(mB0, mA1, result9_0);
                                result9_1 = _mm512_fmadd_pd(mB1, mA1, result9_1);
                            }

                            *((__m512d *)(&C[(ic + ir + 0) * M + jc + jr])) += result0_0;
                            *((__m512d *)(&C[(ic + ir + 0) * M + jc + jr + 8])) += result0_1;
                            *((__m512d *)(&C[(ic + ir + 1) * M + jc + jr])) += result1_0;
                            *((__m512d *)(&C[(ic + ir + 1) * M + jc + jr + 8])) += result1_1;
                            *((__m512d *)(&C[(ic + ir + 2) * M + jc + jr])) += result2_0;
                            *((__m512d *)(&C[(ic + ir + 2) * M + jc + jr + 8])) += result2_1;
                            *((__m512d *)(&C[(ic + ir + 3) * M + jc + jr])) += result3_0;
                            *((__m512d *)(&C[(ic + ir + 3) * M + jc + jr + 8])) += result3_1;
                            *((__m512d *)(&C[(ic + ir + 4) * M + jc + jr])) += result4_0;
                            *((__m512d *)(&C[(ic + ir + 4) * M + jc + jr + 8])) += result4_1;
                            *((__m512d *)(&C[(ic + ir + 5) * M + jc + jr])) += result5_0;
                            *((__m512d *)(&C[(ic + ir + 5) * M + jc + jr + 8])) += result5_1;
                            *((__m512d *)(&C[(ic + ir + 6) * M + jc + jr])) += result6_0;
                            *((__m512d *)(&C[(ic + ir + 6) * M + jc + jr + 8])) += result6_1;
                            *((__m512d *)(&C[(ic + ir + 7) * M + jc + jr])) += result7_0;
                            *((__m512d *)(&C[(ic + ir + 7) * M + jc + jr + 8])) += result7_1;
                            *((__m512d *)(&C[(ic + ir + 8) * M + jc + jr])) += result8_0;
                            *((__m512d *)(&C[(ic + ir + 8) * M + jc + jr + 8])) += result8_1;
                            *((__m512d *)(&C[(ic + ir + 9) * M + jc + jr])) += result9_0;
                            *((__m512d *)(&C[(ic + ir + 9) * M + jc + jr + 8])) += result9_1;
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

    // data validation (for correctness in this program)
    if (M % COL_BLOCK != 0)
    {
        fprintf(stderr, "Matrix dimension M must be a multiple of %d", COL_BLOCK);
        exit(-1);
    }

    if (N % MC != 0 || MC % ROW_BLOCK != 0)
    {
        fprintf(stderr, "Matrix dimension N must be a multiple of %d", ROW_BLOCK);
        exit(-1);
    }

    if (P % KC != 0)
    {
        fprintf(stderr, "Matrix dimension P must be a multiple of %d", KC);
        exit(-1);
    }

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
    write_data_to_file_avx("out/results-avx.csv", "avx-omp-divisible", N, P, M, KC, MC, ROW_BLOCK, COL_BLOCK, num_threads, ALIGNMENT, wctime, wctime_blas, error);

    _mm_free(A);
    _mm_free(B);
    _mm_free(C);
    _mm_free(BLAS_C);
}
