#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <timing.h>
#include <utilities.h>

/**
 * Returns the wall clock time elapsed in the naive matrix multiplication between A and B
 *
 * @param A double[M*N] (row wise storage) -> first N entries are row 1
 * @param B double[N*K] (column wise storage) -> first N entries are column 1
 * @param C double[M*k] (row wise storage) -> first N entries are row 1
 * @param M Number of rows of A
 * @param N Number of columns of A (rows of B)
 * @param K Number or columns of B
 */
double matrix_multiply_naive(double *A, double *B, double *C, int M, int N, int K)
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
            int iC = i * N + j;

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

int main(int argc, char **argv)
{
    double *A, *B, *C;
    size_t sizes[5] = {1000, 2000, 4000, 7633, 8000};
    char files[5][75] = {"/gpfs/gibbs/project/cpsc424/shared/assignments/assignment3/data/C-1000.dat",
                         "/gpfs/gibbs/project/cpsc424/shared/assignments/assignment3/data/C-2000.dat",
                         "/gpfs/gibbs/project/cpsc424/shared/assignments/assignment3/data/C-4000.dat",
                         "/gpfs/gibbs/project/cpsc424/shared/assignments/assignment3/data/C-7633.dat",
                         "/gpfs/gibbs/project/cpsc424/shared/assignments/assignment3/data/C-8000.dat"};

    // Print a table heading
    printf("Matrix multiplication times:\n   N      TIME (secs)    F-norm of Error\n -----   -------------  -----------------\n");

    // Now run the four test cases
    for (int run = 0; run < 5; run++)
    {
        srand(12345);
        size_t size = sizes[run];
        int size_A = size * size;
        int size_B = size * size;
        int size_C = size * size;

        A = (double *)calloc(size_A, sizeof(double));
        B = (double *)calloc(size_B, sizeof(double));
        C = (double *)calloc(size_C, sizeof(double));

        for (int i = 0; i < size_A; i++)
        {
            A[i] = ((double)rand() / (double)RAND_MAX);
        }
        for (int i = 0; i < size_B; i++)
        {
            B[i] = ((double)rand() / (double)RAND_MAX);
        }

        double wctime = matrix_multiply_naive(A, B, C, size, size, size);

        free(A);
        free(B);

        double Fnorm = compute_fnorm(files[run], C, size_C);

        // Print a table row
        printf("  %5d    %9.4f  %17.12f\n", size, wctime, Fnorm);
        write_data_to_file("out/results-serial.csv", "a-serial", size, 1, 1, wctime, Fnorm);

        free(C);
    }
}
