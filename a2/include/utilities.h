#ifndef UTILITIES_H
#define UTILITIES_H

double compute_mandelbrot_area_estimate(int n_i, int total_it);

/**
 * Get the integer value of an environment variable if possible
 */
int get_environment_value(const char *env_name);

/**
 * write runs to csv for tabulation
 */
void write_data_to_file(const char *filename,
                        const char *program,
                        unsigned seed,
                        float wc_time,
                        float area);

#endif