#include "drand.h"

#ifndef UTILITIES_H
#define UTILITIES_H

static inline double get_random_double_in_bounds(double min, double max)
{
    return min + drand() * (max - min);
}

static inline double get_random_double_in_bounds_ts(double min, double max)
{
    return min + drand_ts() * (max - min);
}

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
                        int num_cores,
                        int num_threads,
                        unsigned seed,
                        float wc_time,
                        float area);

#endif