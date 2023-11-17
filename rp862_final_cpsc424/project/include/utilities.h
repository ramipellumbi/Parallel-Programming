
/**
 * Fill double[size] with random values
 *
 * @param M double[size] array
 * @param size of the M array
 */
void fill_matrix(double *M, int size);

/**
 * Compute the Frobenius norm betweeen the known result matrix stored at filename
 * and the computed matrix
 */
double compute_fnorm(const char *filename, double *computed, int size_computed);

/**
 * write runs to csv for tabulation
 */
void write_data_to_file(const char *filename,
                        const char *program,
                        int N,
                        int np,
                        double exe_time,
                        double f_norm);