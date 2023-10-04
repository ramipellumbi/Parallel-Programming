#include <stdio.h>

#include "mandelbrot.h"
#include "timing.h"
#include "utilities.h"

static const double CELL_SIDE_LENGTH = 0.001;

int main(int argc, char *argv[])
{
    int num_cores = get_environment_value("SLURM_CPUS_PER_TASK");
    if (num_cores == -1)
    {
        fprintf(stderr, "SLURM_CPUS_PER_TASK not set");
        return 0;
    }

    // set seed
    unsigned seed = 144545;
    dsrand(seed);

    // initialize timing measures
    double start_wc_time = 0.0, end_wc_time = 0.0;
    double start_cpu_time = 0.0, end_cpu_time = 0.0;

    int number_of_cells_inside_mandelbrot_set = 0, total_iterations = 0;

    // start the timing
    timing(&start_wc_time, &start_cpu_time);
    for (size_t n = 0; n < 2500; n++)
    {
        double current_bottom_left_x = -2.0 + CELL_SIDE_LENGTH * n;
        double max_x = current_bottom_left_x + CELL_SIDE_LENGTH;
        for (size_t m = 0; m < 1250; m++)
        {
            double current_bottom_left_y = 0.0 + CELL_SIDE_LENGTH * m;
            double max_y = current_bottom_left_y + CELL_SIDE_LENGTH;
            double random_x = get_random_double_in_bounds(current_bottom_left_x, max_x);
            double random_y = get_random_double_in_bounds(current_bottom_left_y, max_y);

            int increment = mandelbrot_iteration(random_x, random_y);

            if (increment == 1)
            {
                write_success_to_file("test.csv", "seq", random_x, random_y, n, m);
            }

            number_of_cells_inside_mandelbrot_set += increment;
            total_iterations += 1;
        }
    }
    // end the timing
    timing(&end_wc_time, &end_cpu_time);

    double elapsed_wc_time = end_wc_time - start_wc_time;
    double elapsed_cpu_time = end_cpu_time - start_cpu_time;

    int number_of_cells_outside_mandelbrot_set = total_iterations - number_of_cells_inside_mandelbrot_set;
    double area_of_grid = 1.25 * 2.5;
    double ratio = (double)number_of_cells_inside_mandelbrot_set / (number_of_cells_inside_mandelbrot_set + number_of_cells_outside_mandelbrot_set);
    double area = 2.0 * area_of_grid * ratio;

    write_data_to_file("out/serial.csv",
                       "serial",
                       num_cores,
                       1,
                       seed,
                       elapsed_wc_time,
                       area);

    printf("\narea estimate = %f\n", area);
    printf("Num inside: %d\n", number_of_cells_inside_mandelbrot_set);
    printf("Num outside: %d\n", number_of_cells_outside_mandelbrot_set);
    printf("elapsed wall clock time = %f\n", elapsed_wc_time);
    printf("elapsed cpu time = %f\n", elapsed_cpu_time);

    return 0;
}