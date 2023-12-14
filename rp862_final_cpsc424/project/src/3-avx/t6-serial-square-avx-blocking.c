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

double compute_relative_error_2(double *A, double *B, double *C, double *C2, size_t N, size_t P, size_t M)
{
    double error, suma, sumb, sumc, ai, bi, ci;
    suma = 0.;
    sumb = 0;
    sumc = 0;
    for (size_t i = 0; i < N * P; i++)
    {
        ai = A[i];
        suma += ai * ai;
    }
    for (size_t i = 0; i < P * M; i++)
    {
        bi = B[i];
        sumb += bi * bi;
    }
    for (size_t i = 0; i < N * M; i++)
    {
        ci = C[i] - C2[i];
        sumc += ci * ci;
    }
    suma = sqrt(suma);
    sumb = sqrt(sumb);
    sumc = sqrt(sumc);
    error = sumc / (suma * sumb);

    return error;
}

/**
 * Returns the wall clock time elapsed in the naive matrix multiplication between A and B
 *
 * @param A double[N*P] (row wise storage) -> first P entries are row 1
 * @param B double[P*M] (row wise storage) -> first M entries are row 1
 * @param C double[N*M] (row wise storage) -> first M entries are row 1
 * @param N Number of rows of A
 * @param P Number of columns of A (rows of B)
 * @param M Number or columns of B
 */
double matrix_multiply_blocking(const double *__restrict__ A, const double *__restrict__ B, double *C, int N, int P, int M)
{
    double wc_start, wc_end;
    double cpu_start, cpu_end;
    int iA, jB, iC;

    // ikj loop (same cache analysis as kij)
    const int nc = M;
    const int kc = 240; // kc * 16 doubles fit in L1, 32K
    const int mc = 120; // kc * mc doubles fit in L2, 1024K
    // note: kc * nc fits in L3, which is 35.75 M BUT SHARED AMONGST ALL CORES

    // consider blocks of 6x32 (32 because 512 bit registers fit 4 doubles)
    const int nr = 16;
    const int mr = 6;

    timing(&wc_start, &cpu_start);
    // jc = 0,...,M-1 in steps of nc
    for (size_t jc = 0; jc < M; jc += nc)
    {
        // pc = 0,...,P-1 in steps of kc
        for (size_t pc = 0; pc < P; pc += kc)
        {
            // ic = 0,...,N-1 in steps of mc
            for (size_t ic = 0; ic < N; ic += mc)
            {
                // jr = 0,...,nc-1 in steps of cols
                for (size_t jr = 0; jr < nc; jr += nr)
                {
                    // ir = 0,...,mc-1 in steps of rows
                    for (size_t ir = 0; ir < mc; ir += mr)
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

                        __m512d result0_1 = _mm512_set1_pd(0);
                        __m512d result1_1 = _mm512_set1_pd(0);
                        __m512d result2_1 = _mm512_set1_pd(0);
                        __m512d result3_1 = _mm512_set1_pd(0);
                        __m512d result4_1 = _mm512_set1_pd(0);
                        __m512d result5_1 = _mm512_set1_pd(0);

                        for (size_t k = 0; k < kc; k++)
                        {
                            // load 16 consecutive doubles from k'th row of B (8 to mB0 and 8 to mB1)
                            mB0 = _mm512_load_pd(&B[M * (k + pc) + jc + jr + 8 * 0]);
                            mB1 = _mm512_load_pd(&B[M * (k + pc) + jc + jr + 8 * 1]);

                            mA0 = _mm512_set1_pd(A[k + pc + (ic + ir + 0) * P]);
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
                        }

                        *((__m512d *)(&C[(ic + ir + 0) * M + jc + jr + 0 * 8])) += result0_0;
                        *((__m512d *)(&C[(ic + ir + 0) * M + jc + jr + 1 * 8])) += result0_1;
                        *((__m512d *)(&C[(ic + ir + 1) * M + jc + jr + 0 * 8])) += result1_0;
                        *((__m512d *)(&C[(ic + ir + 1) * M + jc + jr + 1 * 8])) += result1_1;
                        *((__m512d *)(&C[(ic + ir + 2) * M + jc + jr + 0 * 8])) += result2_0;
                        *((__m512d *)(&C[(ic + ir + 2) * M + jc + jr + 1 * 8])) += result2_1;
                        *((__m512d *)(&C[(ic + ir + 3) * M + jc + jr + 0 * 8])) += result3_0;
                        *((__m512d *)(&C[(ic + ir + 3) * M + jc + jr + 1 * 8])) += result3_1;
                        *((__m512d *)(&C[(ic + ir + 4) * M + jc + jr + 0 * 8])) += result4_0;
                        *((__m512d *)(&C[(ic + ir + 4) * M + jc + jr + 1 * 8])) += result4_1;
                        *((__m512d *)(&C[(ic + ir + 5) * M + jc + jr + 0 * 8])) += result5_0;
                        *((__m512d *)(&C[(ic + ir + 5) * M + jc + jr + 1 * 8])) += result5_1;
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

    // ------------------- data load -----------------
    srand(12345);
    int N = atoi(argv[1]);
    int P = atoi(argv[2]);
    int M = atoi(argv[3]);

    double *A, *B, *C, *BLAS_C;
    double wc_start, wc_end;
    double cpu_start, cpu_end;

    A = (double *)_mm_malloc(N * P * sizeof(double), 512);
    B = (double *)_mm_malloc(P * M * sizeof(double), 512);
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

    BLAS_C = (double *)_mm_malloc(N * M * sizeof(double), 512);
    C = (double *)_mm_malloc(N * M * sizeof(double), 512);
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
    double error = compute_relative_error_2(A, B, C, BLAS_C, N, P, M);

    // Print a table row
    printf("\n(%d, %d, %d) %9.4f  %f\n", N, P, M, wctime, error);
    write_data_to_file("out/results-omp.csv", "f-avx-blocking", N, P, M, 6, 1, wctime, wctime_blas, error);

    _mm_free(A);
    _mm_free(B);
    _mm_free(C);
    _mm_free(BLAS_C);
}
