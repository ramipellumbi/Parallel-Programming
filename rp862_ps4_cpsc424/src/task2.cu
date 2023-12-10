#define FP double
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
 *
 * This function was modified and extensively based on Figure 4.16
 * and Figure 4.20 along with section 4.7 in Chapter 4 of Kirk and Hwu's
 * `Programming Massively Parallel Processors`
 */
__global__ void gpu_matrixmult_rectangular_shared(FP *A, FP *B, FP *C, int n, int p, int m, int TILE_WIDTH)
{
    // the textbook uses two double arrays - I interpreted not copying
    // the data structure as continuing to use single arrays
    extern __shared__ FP tiles[];
    FP *Ads = &tiles[0];
    FP *Bds = &tiles[TILE_WIDTH * TILE_WIDTH];

    int bx = blockIdx.x;
    int by = blockIdx.y;
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    int row = by * TILE_WIDTH + ty;
    int col = bx * TILE_WIDTH + tx;

    int tile_idx = ty * TILE_WIDTH + tx;

    FP Cvalue = 0;
    // Loop over the A and B tiles required to compute the C element
    for (size_t ph = 0; ph < ceil((double)p / (double)TILE_WIDTH); ph++)
    {
        // collaborative loading of A and B tiles into shared memory
        int col_bound_A = ph * TILE_WIDTH + tx;
        int row_bound_B = ph * TILE_WIDTH + ty;

        int indexa = row * p + col_bound_A;
        int indexb = row_bound_B * m + col;
        if (row < n && col_bound_A < p)
        {
            Ads[tile_idx] = A[indexa];
        }
        if (col < m && row_bound_B < p)
        {
            Bds[tile_idx] = B[indexb];
        }
        __syncthreads();

        for (size_t k = 0; k < TILE_WIDTH; k++)
        {
            // row ty column k of A with row k column tx of B
            Cvalue += Ads[ty * TILE_WIDTH + k] * Bds[k * TILE_WIDTH + tx];
        }
        __syncthreads();
    }
    if (col < m && row < n)
    {
        C[row * m + col] = Cvalue;
    }
}

// void cpu_matrixmult_rectangular_kij(FP *A, FP *B, FP *C, int n, int p, int m)
// {
//     FP r;
//     // the dot product between a row of A and column of B is between p numbers
//     for (size_t k = 0; k < p; k++)
//     {
//         // there are n rows in A
//         for (size_t i = 0; i < n; i++)
//         {
//             size_t ia = i * p + k; // row i column k of A
//             r = A[ia];

//             // there are m columns in b
//             for (size_t j = 0; j < m; j++)
//             {
//                 size_t ib = k * m + j; // row k column j of B
//                 size_t ic = i * m + j; // row i column j of C

//                 C[ic] += r * B[ib];
//             }
//         }
//     }
// }

int main(int argc, char **argv)
{

    FP *A, *B, *C; // matrices on host
    FP *HOST_C;    
    FP *dev_A, *dev_B, *dev_C; // matrices on device
    int n, p, m;               // matrix dimensions

    int gpucount = 0;   // Count of available GPUs
    int gpunum = 0;     // Device number to use
    int Grid_Dim_x = 1; // Grid dimension, x and y
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
    int TILE_WIDTH = Block_Dim_x;
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

    A = (FP *)malloc(size_A); // dynamically allocated memory for arrays on host
    B = (FP *)malloc(size_B); // dynamically allocated memory for arrays on host
    C = (FP *)malloc(size_C); // results from GPU
    // HOST_C = (FP *)malloc(size_C);

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
            C[i * m + j] = 0.;
            // HOST_C[i * m + j] = 0.;
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

    size_t TW = 2 * TILE_WIDTH * TILE_WIDTH * sizeof(FP);
    gpu_matrixmult_rectangular_shared<<<Grid, Block, TW>>>(dev_A, dev_B, dev_C, n, p, m, TILE_WIDTH);

    cudaEventRecord(stop, 0); // instrument code to measure end time
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsed_time_ms, start, stop);

    cudaMemcpy(C, dev_C, size_C, cudaMemcpyDeviceToHost);

    printf("Time to calculate results on GPU: %f ms.\n", elapsed_time_ms); // exec. time
    write_data_to_file("out/task2.csv", "task2", "gpu", n, p, m, Block_Dim_x, Block_Dim_y, Grid_Dim_x, Grid_Dim_y, elapsed_time_ms);

    // // ------------- COMPUTATION DONE ON HOST CPU ----------------------------
    // // DEBUGGING USE ONLY (AND FOR LIMITED NUMBERS OF TIMING RUNS)

    // cudaEventRecord(start, 0); // use same timing
    // // cudaEventSynchronize(start); // not needed

    // cpu_matrixmult_rectangular_kij(A, B, HOST_C, n, p, m); // do calculation on host (NOTE: This computes the diff with GPU result.)

    // cudaEventRecord(stop, 0); // instrument code to measue end time
    // cudaEventSynchronize(stop);
    // cudaEventElapsedTime(&elapsed_time_ms, start, stop);

    // printf("Time to calculate results on CPU: %f ms.\n", elapsed_time_ms); // exec. time

    // // ------------------- check device creates correct results -----------------

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
    //     ci = (double)C[i] - (double)HOST_C[i];
    //     sumc += ci * ci;
    // }
    // suma = sqrt(suma);
    // sumb = sqrt(sumb);
    // sumc = sqrt(sumc);
    // error = sumc / (suma * sumb);
    // printf("Approximate relative error between GPU and CPU: %e\n", error);
    // END OF OPTIONAL SECTION THAT CAN BE OMITTED
    // -------------- clean up ---------------------------------------

    // -------------- clean up ---------------------------------------
    free(A);
    free(B);
    free(C);
    // free(HOST_C);
    cudaFree(dev_A);
    cudaFree(dev_B);
    cudaFree(dev_C);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    return 0;
}