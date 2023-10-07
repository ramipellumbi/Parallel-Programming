#include <immintrin.h>
#include <nmmintrin.h>

#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include "drand-ts.h"
#include "mandelbrot.h"
#include "timing.h"
#include "utilities.h"

static const size_t NUM_X_ITERATIONS = 2500;
static const size_t NUM_Y_ITERATIONS = 1250;
static const size_t MAX_ITERATIONS = 25000;
static const double CELL_SIDE_LENGTH = 0.001;
static const int PACKING_SIZE = 8;

int main(int argc, char *argv[])
{
    // number of iterations for y that is a multiple of PACKING_SIZE
    int NUM_Y_PS = NUM_Y_ITERATIONS / PACKING_SIZE * PACKING_SIZE;

    // initialize timing measures
    double start_wc_time = 0.0, end_wc_time = 0.0;
    double start_cpu_time = 0.0, end_cpu_time = 0.0;

    int number_of_cells_inside_mandelbrot_set = 0, total_iterations = 0;

    // We have 8 packed doubles in a 512 bit register
    // This stores the corner offsets - 0, 0.001, 0.002, ..., 0.007
    // so that at every x value we can assess 8 y values at a time
    __m512d pxs_deltas512 = _mm512_mul_pd(
        _mm512_set_pd(0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0),
        _mm512_set1_pd(CELL_SIDE_LENGTH));

    // set seed
    unsigned seed = 144545;
    dsrand_ts(seed);

    // arrays to store the 16 needed random numbers (8 real, 8 imaginary) each iteration
    double *random_x;
    double *random_y;

    // start the timing
    timing(&start_wc_time, &start_cpu_time);

#pragma omp parallel shared(number_of_cells_inside_mandelbrot_set, total_iterations, pxs_deltas512, NUM_Y_PS) private(random_x, random_y) default(none)
    {
        random_x = (double *)malloc(PACKING_SIZE * sizeof(double));
        random_y = (double *)malloc(PACKING_SIZE * sizeof(double));

// for each x value
#pragma omp single
        {
            for (size_t n = 0; n < NUM_X_ITERATIONS; n++)
            {
                double current_bottom_left_x = -2.0 + CELL_SIDE_LENGTH * n;
                double max_x = current_bottom_left_x + CELL_SIDE_LENGTH;

#pragma omp task firstprivate(current_bottom_left_x, max_x, random_x, random_y)
                {
                    // let's compute 8 y-values simultaneously
                    for (size_t m = 0; m < NUM_Y_PS; m += PACKING_SIZE)
                    {
                        for (int i = 0; i < PACKING_SIZE; ++i)
                        {
                            random_x[i] = drand_ts();
                            random_y[i] = drand_ts();
                        }
                        // grab 8 random numbers for the 8 needed random x coordinates
                        __m512d random_numbers_x = _mm512_loadu_pd(random_x);
                        // grab 8 random numbers for the 8 needed random y coordinates
                        __m512d random_numbers_y = _mm512_loadu_pd(random_y);

                        // get the 8 bottom left corners for this iteration of the loop
                        __m512d bottom_left_y_values = _mm512_add_pd(
                            _mm512_set1_pd(0.0),
                            _mm512_add_pd(_mm512_set1_pd(m * CELL_SIDE_LENGTH), pxs_deltas512));
                        // get the 8 top left corners for this iteration of the loop
                        __m512d top_left_y_values = _mm512_add_pd(
                            _mm512_set1_pd(CELL_SIDE_LENGTH),
                            bottom_left_y_values);

                        // compute the random x coordinates for the 8 cells
                        __m512d x_values = _mm512_fmadd_pd(
                            random_numbers_x,
                            _mm512_set1_pd(max_x - current_bottom_left_x),
                            _mm512_set1_pd(current_bottom_left_x));

                        // compute the random y coordinates for the 8 cells
                        __m512d y_values = _mm512_fmadd_pd(
                            random_numbers_y,
                            _mm512_sub_pd(top_left_y_values, bottom_left_y_values),
                            bottom_left_y_values);

                        // These are the 8 z values
                        __m512d z_re = _mm512_set1_pd(0.0);
                        __m512d z_im = _mm512_set1_pd(0.0);

                        __mmask8 diverged_indices = mandelbrot_iteration_avx(
                            x_values,
                            y_values,
                            MAX_ITERATIONS);

                        // the 1's in this mask are the iterations that did NOT diverge
                        __mmask8 indices_in_set = ~diverged_indices;
                        int count = _popcnt32((unsigned int)indices_in_set);
#pragma omp atomic
                        number_of_cells_inside_mandelbrot_set += count;

#pragma omp atomic
                        total_iterations += PACKING_SIZE;
                    }

                    // cleanup by performing the naive serial implementation
                    for (size_t m = NUM_Y_PS; m < NUM_Y_ITERATIONS; m++)
                    {
                        double current_bottom_left_y = 0.0 + CELL_SIDE_LENGTH * m;
                        double max_y = current_bottom_left_y + CELL_SIDE_LENGTH;

                        double c_re = current_bottom_left_x + drand_ts() * (max_x - current_bottom_left_x);
                        double c_im = current_bottom_left_y + drand_ts() * (max_y - current_bottom_left_y);

                        int counter = mandelbrot_iteration(c_re, c_im, MAX_ITERATIONS);

#pragma omp atomic
                        number_of_cells_inside_mandelbrot_set += counter;
#pragma omp atomic
                        total_iterations++;
                    }
                }
            }
        }

        // free the arrays for holding x and y
        free(random_x);
        free(random_y);
    }

    // end the timing
    timing(&end_wc_time, &end_cpu_time);

    double elapsed_wc_time = end_wc_time - start_wc_time;
    double elapsed_cpu_time = end_cpu_time - start_cpu_time;
    double area = compute_mandelbrot_area_estimate(number_of_cells_inside_mandelbrot_set,
                                                   total_iterations);

    write_data_to_file("out/tasks.csv",
                       "per-column-task",
                       seed,
                       elapsed_wc_time,
                       area);
    printf("\narea estimate = %f\n", area);
    printf("elapsed wall clock time = %f\n", elapsed_wc_time);
    printf("elapsed cpu time = %f\n", elapsed_cpu_time);

    return 0;
}