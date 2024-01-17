#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

int get_environment_value(const char *env_name)
{
    char *value_str = getenv(env_name);
    if (value_str == NULL)
    {
        return -1;
    }
    int value = atoi(value_str);

    return value;
}

void write_data_to_file(const char *filename,
                        const char *program,
                        int N,
                        int P,
                        int M,
                        int block_size,
                        int np,
                        double exe_time,
                        double blas_exe_time,
                        double f_norm)
{
    int num_cores = get_environment_value("SLURM_CPUS_PER_TASK");
    if (num_cores == -1)
    {
        fprintf(stderr, "\nSLURM_CPUS_PER_TASK not set");
    }

    int num_tasks_per_node = get_environment_value("SLURM_NTASKS_PER_NODE");
    if (num_tasks_per_node == -1)
    {
        fprintf(stderr, "\nSLURM_NTASKS_PER_NODE not set");
    }

    int num_tasks_per_socket = get_environment_value("SLURM_NTASKS_PER_SOCKET");
    if (num_tasks_per_socket == -1)
    {
        fprintf(stderr, "\nSLURM_NTASKS_PER_SOCKET not set");
    }

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
        fprintf(fp, "program,num_cores,np,num_tasks_per_node,num_tasks_per_socket,N,P,M,block_size,exe_time,blas_exe_time,f_norm\n");
    }

    // add data to row
    fprintf(fp, "\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%f,%f,%.12e\n", program, num_cores, np, num_tasks_per_node, num_tasks_per_socket, N, P, M, block_size, exe_time, blas_exe_time, f_norm);
    fclose(fp);
}

void write_data_to_file_avx(const char *filename,
                            const char *program,
                            int N,
                            int P,
                            int M,
                            int KC,
                            int MC,
                            int ROW_BLOCK,
                            int COL_BLOCK,
                            int np,
                            int ALIGNMENT,
                            double exe_time,
                            double blas_exe_time,
                            double f_norm)
{
    int num_cores = get_environment_value("SLURM_CPUS_PER_TASK");
    if (num_cores == -1)
    {
        fprintf(stderr, "\nSLURM_CPUS_PER_TASK not set");
    }

    int num_tasks_per_node = get_environment_value("SLURM_NTASKS_PER_NODE");
    if (num_tasks_per_node == -1)
    {
        fprintf(stderr, "\nSLURM_NTASKS_PER_NODE not set");
    }

    int num_tasks_per_socket = get_environment_value("SLURM_NTASKS_PER_SOCKET");
    if (num_tasks_per_socket == -1)
    {
        fprintf(stderr, "\nSLURM_NTASKS_PER_SOCKET not set");
    }

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
        fprintf(fp, "program,num_cores,np,num_tasks_per_node,num_tasks_per_socket,N,P,M,KC,MC,ROW_BLOCK,COL_BLOCK,ALIGNMENT,exe_time,blas_exe_time,f_norm\n");
    }

    // add data to row
    fprintf(fp, "\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f,%f,%.12e\n", program, num_cores, np, num_tasks_per_node, num_tasks_per_socket, N, P, M, KC, MC, ROW_BLOCK, COL_BLOCK, ALIGNMENT, exe_time, blas_exe_time, f_norm);
    fclose(fp);
}

/**
 * Load random matrices A and B
 * Both are stored in a single array row-wise
 *
 * A is n x p
 * B is p x m
 */
void load_random_matrices(double **A, double **B, size_t n, size_t p, size_t m)
{
    *A = (double *)malloc(n * p * sizeof(double));
    *B = (double *)malloc(p * m * sizeof(double));

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < p; j++)
        {
            int idx = i * p + j;
            (*A)[idx] = ((double)rand() / (double)RAND_MAX);
        }
    }

    // load B row-wise
    for (int i = 0; i < p; i++)
    {
        for (int j = 0; j < m; j++)
        {
            int idx = i * m + j;
            (*B)[idx] = ((double)rand() / (double)RAND_MAX);
        }
    }
}

double compute_relative_error(double *A, double *B, double *C, size_t N, size_t P, size_t M)
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
        ci = C[i];
        sumc += ci * ci;
    }
    suma = sqrt(suma);
    sumb = sqrt(sumb);
    sumc = sqrt(sumc);
    error = sumc / (suma * sumb);

    return error;
}

double compute_relative_error_between(double *A, double *B, double *C, double *C2, size_t N, size_t P, size_t M)
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
