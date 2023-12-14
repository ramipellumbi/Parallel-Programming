#ifndef UTILITIES_H
#define UTILITIES_H

int get_environment_value(const char *env_name);

/**
 * write runs to csv for tabulation
 */
void write_data_to_file(const char *filename,
                        const char *program,
                        int N,
                        int P,
                        int M,
                        int block_size,
                        int np,
                        double exe_time,
                        double blas_exe_time,
                        double f_norm);

void load_random_matrices(double **A, double **B, size_t n, size_t p, size_t m);

double compute_relative_error(double *A, double *B, double *C, size_t n, size_t p, size_t m);

double compute_relative_error_between(double *A, double *B, double *C, double *C2, size_t n, size_t p, size_t m)

#endif