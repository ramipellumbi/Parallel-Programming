#define FP float
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

void write_data_to_file(const char *filename,   // name of csv
                        const char *program,    // 'task1', 'task2', etc.
                        const char *multiplier, // 'cpu' or 'gpu'
                        int n,
                        int p,
                        int m,
                        int block_x,
                        int block_y,
                        int grid_x,
                        int grid_y,
                        float exe_time)
{
    // check if the file exists
    FILE *check_file = fopen(filename, "r");
    bool does_file_exist = false;
    if (check_file != NULL)
    {
        fclose(check_file);
        does_file_exist = true;
    }

    // open the file with means to append
    FILE *fp = fopen(filename, "a");
    if (fp == NULL)
    {
        fprintf(stderr, "Could not open or create file: %s\n", filename);
        return;
    }

    // if the file does not exist, add the header row
    if (!does_file_exist)
    {
        fprintf(fp, "program,multiplier,precision,n,p,m,block_x,block_y,grid_x,grid_y,exe_time\n");
    }

    // add data to row
    fprintf(fp, "\"%s\",\"%s\",\"%s\",%d,%d,%d,%d,%d,%d,%d,%f\n", program, multiplier, TOSTRING(FP), n, p, m, block_x, block_y, grid_x, grid_y, exe_time);
    fclose(fp);
}

/**
 * @param A is n x p - STORED ROW-WISE
 * @param B is p x m - STORED ROW-WISE
 * @param C is n x m - STORED ROW-WISE
 */
void cpu_matrixmult_rectangular_kij(FP *A, FP *B, FP *C, int n, int p, int m)
{
    FP r;
    // the dot product between a row of A and column of B is between p numbers
    for (size_t k = 0; k < p; k++)
    {
        // there are n rows in A
        for (size_t i = 0; i < n; i++)
        {
            size_t ia = i * p + k; // row i column k of A
            r = A[ia];

            // there are m columns in b
            for (size_t j = 0; j < m; j++)
            {
                size_t ib = k * m + j; // row k column j of B
                size_t ic = i * m + j; // row i column j of C

                C[ic] += r * B[ib];
            }
        }
    }
}

/**
 * Print matrix M (stored row wise)
 */
void print_matrix(FP *M, const char *name, size_t nrow, size_t ncol)
{
    printf("Matrix %s\n", name);
    printf("[\n");
    for (size_t i = 0; i < nrow; i++)
    {
        printf("[");
        for (size_t j = 0; j < ncol; j++)
        {
            printf("%f,", M[i * ncol + j]);
        }
        printf("],\n");
    }
    printf("]\n");
}

int main(int argc, char **argv)
{

    FP *A, *B, *C;
    int n, p, m;

    cudaEvent_t start, stop; // using cuda events to measure time
    float elapsed_time_ms;

    if (argc != 4)
    {
        printf("Usage: task1 <matrix dim n> <matrix dim p> <matrix dim m>\n");
        exit(-1);
    }

    n = atoi(argv[1]);
    p = atoi(argv[2]);
    m = atoi(argv[3]);
    printf("%d,%d,%d", n, p, m);

    size_t size_A = n * p * sizeof(FP);
    size_t size_B = p * m * sizeof(FP);
    size_t size_C = n * m * sizeof(FP);

    A = (FP *)malloc(size_A);
    B = (FP *)malloc(size_B);
    C = (FP *)malloc(size_C);

    srand(12345);
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < p; j++)
        {
            A[i * p + j] = (FP)rand() / (FP)RAND_MAX;
        }
    }

    for (size_t i = 0; i < p; i++)
    {
        for (size_t j = 0; j < m; j++)
        {
            B[i * m + j] = (FP)rand() / (FP)RAND_MAX;
        }
    }

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < m; j++)
        {
            C[i * m + j] = 0.0;
        }
    }

    // ------------- COMPUTATION DONE ON HOST CPU ----------------------------
    cudaEventCreate(&start); // instrument code to measure start time
    cudaEventCreate(&stop);

    cudaEventRecord(start, 0);
    cpu_matrixmult_rectangular_kij(A, B, C, n, p, m); // do calculation on host
    cudaEventRecord(stop, 0);                         // instrument code to measue end time
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsed_time_ms, start, stop);

    write_data_to_file("out/task1.csv", "task1a", "cpu", n, p, m, 0, 0, 0, 0, elapsed_time_ms);

    // -------------- clean up ---------------------------------------
    free(A);
    free(B);
    free(C);

    return 0;
}