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
__global__ void gpu_matrixmult_rectangular(FP *A, FP *B, FP *C, int n, int p, int m)
{

    int col = threadIdx.x + blockDim.x * blockIdx.x;
    int row = threadIdx.y + blockDim.y * blockIdx.y;

    int indexb = col;
    int index = row * m + col;

    if (col < m && row < n)
    {
        int start_index_a = row * p;
        int max_index_a = row * p + p;
        for (int indexa = start_index_a; indexa < max_index_a; indexa++)
        {
            C[index] += A[indexa] * B[indexb];
            indexb += m;
        }
    }
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

    FP *A, *B, *C, *HOST_C;    // matrices on host
    FP *dev_A, *dev_B, *dev_C; // matrices on device
    int n, p, m;               // matrix dimensions

    int gpucount = 0;  // Count of available GPUs
    int gpunum = 0;    // Device number to use
    int Grid_Dim_x = 1;  // Grid dimension, x and y
    int Grid_Dim_y = 1;
    int Block_Dim_x = 1; // Block dimension, x and y
    int Block_Dim_y = 1;

    cudaEvent_t start, stop; // using cuda events to measure time
    float elapsed_time_ms;   // which is applicable for asynchronous code also
    cudaError_t errorcode;

    // --------------------SET PARAMETERS AND DATA -----------------------
    errorcode = cudaGetDeviceCount(&gpucount);
    if (errorcode == cudaErrorNoDevice)
    {
        printf("No GPUs are visible\n");
        exit(-1);
    }
    else
    {
        printf("Device count = %d\n", gpucount);
    }

    if (argc != 8)
    {
        printf("Usage: task1 <matrix dim n> <matrix dim p> <matrix dim m> <block dim x> <block dim y> <grid dim x> <grid dim y>\n");
        exit(-1);
    }

    n = atoi(argv[1]);
    p = atoi(argv[2]);
    m = atoi(argv[3]);

    Block_Dim_x = atoi(argv[4]); // Rectangular block
    Block_Dim_y = atoi(argv[5]);
    if (Block_Dim_x * Block_Dim_y > 1024)
    {
        printf("Error, too many threads in block\n");
        exit(-1);
    }

    Grid_Dim_x = atoi(argv[6]);
    Grid_Dim_y = atoi(argv[7]);
    if (Grid_Dim_x * Block_Dim_x < m)
    {
        printf("%d,%d", Grid_Dim_x, Block_Dim_x);
        printf("Error, number of threads in x dimensions less than number of array elements\n");
        exit(-1);
    }
    if (Grid_Dim_y * Block_Dim_y < n)
    {
        printf("Error, number of threads in y dimensions less than number of array elements\n");
        exit(-1);
    }

    cudaSetDevice(gpunum);
    printf("Using device %d\n", gpunum);

    printf("Matrix Dimension = %d\n", n);
    printf("Block_Dim_x = %d, Block_Dim_y = %d, Grid_Dim_x = %d, Grid_Dim_y = %d\n", Block_Dim_x, Block_Dim_y, Grid_Dim_x, Grid_Dim_y);

    dim3 Grid(Grid_Dim_x, Grid_Dim_y);    // Grid structure
    dim3 Block(Block_Dim_x, Block_Dim_y); // Block structure

    size_t size_A = n * p * sizeof(FP);
    size_t size_B = p * m * sizeof(FP);
    size_t size_C = n * m * sizeof(FP);

    A = (FP *)malloc(size_A);      // dynamically allocated memory for arrays on host
    B = (FP *)malloc(size_B);      // dynamically allocated memory for arrays on host
    C = (FP *)malloc(size_C);      // results from GPU
    HOST_C = (FP *)malloc(size_C); // results from rectangular kij algorithm

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
        for (size_t j = 0; j < n; j++)
        {
            B[i * n + j] = (FP)rand() / (FP)RAND_MAX;
        }
    }

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < m; j++)
        {
            C[i * m + j] = 0.;
            HOST_C[i * m + j] = 0;
        }
    }

    // ------------- COMPUTATION DONE ON GPU ----------------------------
    cudaMalloc((void **)&dev_A, size_A); // allocate memory on device
    cudaMalloc((void **)&dev_B, size_B);
    cudaMalloc((void **)&dev_C, size_C);
    cudaMemcpy(dev_A, A, size_A, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_B, B, size_B, cudaMemcpyHostToDevice);

    cudaEventCreate(&start); // instrument code to measure start time
    cudaEventCreate(&stop);

    cudaEventRecord(start, 0);

    gpu_matrixmult_rectangular<<<Grid, Block>>>(dev_A, dev_B, dev_C, n, p, m);

    cudaEventRecord(stop, 0); // instrument code to measure end time
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsed_time_ms, start, stop);

    cudaMemcpy(C, dev_C, size_C, cudaMemcpyDeviceToHost);

    printf("Time to calculate results on GPU: %f ms.\n", elapsed_time_ms); // exec. time
    write_data_to_file("out/task1.csv", "task1a", "gpu", n, p, m, Block_Dim_x, Block_Dim_y, Grid_Dim_x, Grid_Dim_y,elapsed_time_ms);

    // ------------- COMPUTATION DONE ON HOST CPU ----------------------------
    // cudaEventRecord(start, 0); // use same timing

    // cpu_matrixmult_rectangular_kij(A, B, HOST_C, n, p, m); // do calculation on host

    // cudaEventRecord(stop, 0); // instrument code to measue end time
    // cudaEventSynchronize(stop);
    // cudaEventElapsedTime(&elapsed_time_ms, start, stop);

    // printf("Time to calculate results on CPU: %f ms.\n", elapsed_time_ms); // exec. time
    // write_data_to_file("out/task1.csv", "task1a", "cpu", n, p, m, Block_Dim_x, Block_Dim_y, Grid_Dim_x, Grid_Dim_y,elapsed_time_ms);

    // ------------------- check device creates correct results -----------------
    // double error, suma, sumb, sumc, ai, bi, ci;
    // suma = 0.;
    // sumb = 0;
    // sumc = 0;
    // for (size_t i = 0; i < n * p; i++)
    // {
    //     ai = (double)A[i];
    //     suma += ai * ai;
    // }
    // for (size_t i = 0; i < p * m; i++)
    // {
    //     bi = (double)B[i];
    //     sumb += bi * bi;
    // }
    // for (size_t i = 0; i < n * m; i++)
    // {
    //     ci = (double)(C[i] - HOST_C[i]);
    //     sumc += ci * ci;
    // }
    // suma = sqrt(suma);
    // sumb = sqrt(sumb);
    // sumc = sqrt(sumc);
    // error = sumc / (suma * sumb);
    // printf("Approximate relative error between GPU and CPU: %e\n", error);

    // -------------- clean up ---------------------------------------
    free(A);
    free(B);
    free(C);
    free(HOST_C);
    cudaFree(dev_A);
    cudaFree(dev_B);
    cudaFree(dev_C);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    return 0;
}