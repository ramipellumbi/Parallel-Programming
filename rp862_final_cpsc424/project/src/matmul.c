#include <immintrin.h>
#include <nmmintrin.h>
#include <timing.h>

static const int PACKING_SIZE = 8;

double matrix_multiply_naive(double *A, double *B, double *C, int M, int N, int K, int offset)
{
    double wc_start, wc_end;
    double cpu_start, cpu_end;

    timing(&wc_start, &cpu_start);
    for (int i = 0; i < M; i++)
    {
        int iA = i * N;
        for (int j = 0; j < K; j++)
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
    for (int i = 0; i < M; i++)
    {
        int iA = i * N;
        for (int j = 0; j < K; j++)
        {
            int jB = j * N;
            int iC = i * N + j + offset;
            C[iC] = 0;

            __m512d C_accumulated = _mm512_setzero_pd();

            for (int k = 0; k < AVX_ITERS; k += PACKING_SIZE)
            {
                __m512d A_for_multiply = _mm512_load_pd(&A[iA + k]);
                __m512d B_for_multiply = _mm512_load_pd(&B[jB + k]);
                C_accumulated = _mm512_fmadd_pd(A_for_multiply, B_for_multiply, C_accumulated);
            }

            C[iC] += _mm512_reduce_add_pd(C_accumulated);

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