#include <matmul.h>
#include <immintrin.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <utilities.h>

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

        // Allocate size bytes of memory, aligned to the alignment specified in align, and return a pointer to the allocated memory.
        double *A = (double *)_mm_malloc(M * K * sizeof(double), 64);
        double *B = (double *)_mm_malloc(K * N * sizeof(double), 64);
        C = (double *)calloc(size_C, sizeof(double));

        fill_matrix(A, size_A);
        fill_matrix(B, size_B);

        double wctime = matrix_multiply_avx(A, B, C, M, N, K, 0);

        _mm_free(A);
        _mm_free(B);

        double Fnorm = compute_fnorm(files[run], C, size_C);

        // Print a table row
        printf("  %5d    %9.4f  %17.12f\n", N, wctime, Fnorm);
        write_data_to_file("out/results.csv", "b-serial-avx", N, 1, wctime, Fnorm);

        free(C);
    }
}
