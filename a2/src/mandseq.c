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

    // start the timing
    timing(&start_wc_time, &start_cpu_time);
    for (size_t n = 0; n < NUM_X_ITERATIONS; n++)
    {
        double current_bottom_left_x = -2.0 + CELL_SIDE_LENGTH * n;
        double max_x = current_bottom_left_x + CELL_SIDE_LENGTH;
        for (size_t m = 0; m < NUM_Y_ITERATIONS; m++)
        {
            double current_bottom_left_y = 0.0 + CELL_SIDE_LENGTH * m;
            double max_y = current_bottom_left_y + CELL_SIDE_LENGTH;
            double random_x = current_bottom_left_x + (max_x - current_bottom_left_x) * drand();
            double random_y = current_bottom_left_y + (max_y - current_bottom_left_y) * drand();

            int increment = mandelbrot_iteration(random_x, random_y, MAX_ITERATIONS);

            number_of_cells_inside_mandelbrot_set += increment;
            total_iterations += 1;
        }
    }
    // end the timing
    timing(&end_wc_time, &end_cpu_time);

    double elapsed_wc_time = end_wc_time - start_wc_time;
    double elapsed_cpu_time = end_cpu_time - start_cpu_time;

    double area = compute_mandelbrot_area_estimate(number_of_cells_inside_mandelbrot_set,
                                                   total_iterations);

    write_data_to_file("out/serial.csv",
                       "serial",
                       seed,
                       elapsed_wc_time,
                       area);
    printf("\narea estimate = %f\n", area);
    printf("elapsed wall clock time = %f\n", elapsed_wc_time);
    printf("elapsed cpu time = %f\n", elapsed_cpu_time);

    return 0;
}