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

#define MIN(a, b) (((a) < (b)) ? (a) : (b))

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
    const int kc = 240; // Cache block size for K dimension
    const int mc = 480; // Cache block size for M dimension
    const int nc = M;   // Cache block size for N dimension

    timing(&wc_start, &cpu_start);
    for (int jc = 0; jc < M; jc += nc)
    {
        int nc_eff = (jc + nc > M) ? M - jc : nc;

        for (int pc = 0; pc < P; pc += kc)
        {
            int kc_eff = (pc + kc > P) ? P - pc : kc;

            for (int ic = 0; ic < N; ic += mc)
            {
                int mc_eff = (ic + mc > N) ? N - ic : mc;

#pragma omp parallel for
                for (int jr = 0; jr < nc_eff; ++jr)
                {
                    for (int ir = 0; ir < mc_eff; ++ir)
                    {
                        double sum = C[(ic + ir) * M + jc + jr];
                        for (int k = 0; k < kc_eff; ++k)
                        {
                            sum += A[(ic + ir) * P + pc + k] * B[(pc + k) * M + jc + jr];
                        }
                        C[(ic + ir) * M + jc + jr] = sum;
                    }
                }
            }
        }
    }
    timing(&wc_end, &cpu_end);
    return wc_end - wc_start;
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
    write_data_to_file("out/results-avx.csv", "h-avx-blocking", N, P, M, 6, 1, wctime, wctime_blas, error);

    _mm_free(A);
    _mm_free(B);
    _mm_free(C);
    _mm_free(BLAS_C);
}
