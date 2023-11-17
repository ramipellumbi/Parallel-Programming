#include <omp.h>
#include <timing.h>

static const int PACKING_SIZE = 8;

double matrix_multiply_naive(double *A, double *B, double *C, int M, int N, int K, int offset)
{
    double wc_start, wc_end;
    double cpu_start, cpu_end;

    timing(&wc_start, &cpu_start);
    for (int i = 0; i < N; i++)
    {
        int iA = i * N;
        for (int j = 0; j < N; j++)
        {
            int jB = j * N;
            int iC = i * N + j + offset;

            C[iC] = 0;
            for (int k = 0; k < N; k++)
            {
                C[iC] += A[iA + k] * B[jB + k];
            }
        }
    }
    timing(&wc_end, &cpu_end);
    double elapsed_time = wc_end - wc_start;

    return elapsed_time;
}

double matrix_multiply_avx(double *A, double *B, double *C, int M, int N, int K, int offset)
{
    double wc_start, wc_end;
    double cpu_start, cpu_end;
    int AVX_ITERS = N / PACKING_SIZE * PACKING_SIZE;

    timing(&wc_start, &cpu_start);
    for (int i = 0; i < N; i++)
    {
        int iA = i * N;
        for (int j = 0; j < N; j++)
        {
            int jB = j * N;
            int iC = i * N + j + offset;

            for (int k = 0; k < AVX_ITERS; k += PACKING_SIZE)
            {
                __m512d A_for_multiply = __mm512_loadu_pd(A[iA + k]);
                __m512d B_for_multiply = __mm512_loadu_pd(B[jB + k]);
                __m512d C_partial_result = __mm512_fmadd_pd(
                    A_for_multiply,
                    B_for_multiply,
                    __m512_set1_pd(0.0));

                // place the result into iC:iC+PACKING_SIZE
            }

            C[iC] = 0;
            // cleanup loop
            for (int k = AVX_ITERS; k < N; k++)
            {
                C[iC] += A[iA + k] * B[jB + k];
            }
        }
    }

    timing(&wc_end, &cpu_end);
    double elapsed_time = wc_end - wc_start;

    return elapsed_time;
}