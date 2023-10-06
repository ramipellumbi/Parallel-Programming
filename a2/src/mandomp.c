#include <omp.h>
#include <stdio.h>

#include "drand.h"
#include "mandelbrot.h"
#include "timing.h"
#include "utilities.h"

static const size_t NUM_X_ITERATIONS = 2500;
static const size_t NUM_Y_ITERATIONS = 1250;
static const size_t MAX_ITERATIONS = 25000;
static const double CELL_SIDE_LENGTH = 0.001;

int main(int argc, char *argv[])
{
    // set seed
    unsigned seed = 144545;
    dsrand(seed);

    // initialize timing measures
    double start_wc_time = 0.0, end_wc_time = 0.0;
    double start_cpu_time = 0.0, end_cpu_time = 0.0;

    int number_of_cells_inside_mandelbrot_set = 0, total_iterations = 0;
    int number_of_cells_inside_mandelbrot_set_th, total_iterations_th;

    // start the timing
    timing(&start_wc_time, &start_cpu_time);
#pragma omp parallel shared(number_of_cells_inside_mandelbrot_set, total_iterations, seed) private(number_of_cells_inside_mandelbrot_set_th, total_iterations_th) default(none)
    {
        number_of_cells_inside_mandelbrot_set_th = 0;
        total_iterations_th = 0;

#pragma omp for
        for (size_t i = 0; i < NUM_X_ITERATIONS; ++i)
        {
            double current_bottom_left_x = -2.0 + i * CELL_SIDE_LENGTH;
            double cell_max_x = current_bottom_left_x + CELL_SIDE_LENGTH;

            for (size_t j = 0; j < NUM_Y_ITERATIONS; ++j)
            {
                double current_bottom_left_y = 0.0 + j * CELL_SIDE_LENGTH;
                double cell_max_y = current_bottom_left_y + CELL_SIDE_LENGTH;
                double random_x = current_bottom_left_x + (cell_max_x - current_bottom_left_x) * drand();
                double random_y = current_bottom_left_y + (cell_max_y - current_bottom_left_y) * drand();

                int increment = mandelbrot_iteration(random_x, random_y, MAX_ITERATIONS);

                number_of_cells_inside_mandelbrot_set_th += increment;
                total_iterations_th += 1;
            }
        }

// add each threads results to the desired totals
#pragma omp atomic
        number_of_cells_inside_mandelbrot_set += number_of_cells_inside_mandelbrot_set_th;

#pragma omp atomic
        total_iterations += total_iterations_th;
    }
    // end the timing
    timing(&end_wc_time, &end_cpu_time);

    double elapsed_wc_time = end_wc_time - start_wc_time;
    double elapsed_cpu_time = end_cpu_time - start_cpu_time;

    double area = compute_mandelbrot_area_estimate(number_of_cells_inside_mandelbrot_set,
                                                   total_iterations);

    write_data_to_file("out/omp.csv",
                       "omp-not-ts",
                       seed,
                       elapsed_wc_time,
                       area);

    printf("\narea estimate = %f\n", area);
    printf("elapsed wall clock time = %f\n", elapsed_wc_time);
    printf("elapsed cpu time = %f\n", elapsed_cpu_time);

    return 0;
}