#include <mkl.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <timing.h>
#include <utilities.h>

/**
 * Returns the wall clock time elapsed in the matrix multiplication between A and B using kij access pattern
 *
 * @param A double[N*P] (row wise storage) -> first P entries are row 1
 * @param B double[P*M] (row wise storage) -> first M entries are row 1
 * @param C double[N*M] (row wise storage) -> first M entries are row 1
 * @param N Number of rows of A
 * @param P Number of columns of A (rows of B)
 * @param M Number or columns of B
 *
 * This function takes in the correct C as computed by BLAS and detracts from it
 */
double matrix_multiply_kij(double *A, double *B, double *C, int N, int P, int M)
{
    double wc_start, wc_end;
    double cpu_start, cpu_end;

    timing(&wc_start, &cpu_start);
    for (int k = 0; k < P; k++)
    {
        for (int i = 0; i < N; i++)
        {
            // get row i column k of A for reuse
            int iA = i * P + k;
            double r = A[iA];
            for (int j = 0; j < M; j++)
            {
                // row k column j of B
                int iB = k * M + j;
                // row i column j of c
                int iC = i * M + j;
                C[iC] -= r * B[iB];
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
    double wctime = matrix_multiply_kij(A, B, C, N, P, M);

    // ------------------- check custom multiply creates correct results -----------------
    double error = compute_relative_error(A, B, C, N, P, M);

    // Print a table row
    printf("\n(%d, %d, %d) %9.4f  %f\n", N, P, M, wctime, error);
    write_data_to_file("out/results-serial.csv", "b-serial-kij", N, P, M, 1, 1, wctime, wctime_blas, error);

    free(A);
    free(B);
    free(C);
}
