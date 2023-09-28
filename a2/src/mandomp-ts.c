#include <omp.h>
#include <stdio.h>

#include "complex.h"
#include "drand.h"
#include "mandelbrot.h"
#include "timing.h"
#include "utilities.h"
#include "writer.h"

static const double CELL_SIDE_LENGTH = 0.001;

int main(int argc, char *argv[])
{
    int num_cores = get_environment_value("SLURM_CPUS_PER_TASK");
    if (num_cores == -1)
    {
        fprintf(stderr, "SLURM_CPUS_PER_TASK not set");
        return 0;
    }

    int num_threads = get_environment_value("OMP_NUM_THREADS");
    if (num_threads == -1)
    {
        fprintf(stderr, "OMP_NUM_THREADS not set");
        return 0;
    }

    // initialize timing measures
    double start_wc_time = 0.0, end_wc_time = 0.0;
    double start_cpu_time = 0.0, end_cpu_time = 0.0;

    int number_of_cells_inside_mandelbrot_set = 0, total_iterations = 0;
    int max_x = (int)(2.5 / CELL_SIDE_LENGTH);
    int max_y = (int)(1.25 / CELL_SIDE_LENGTH);

    
    unsigned seed = 144545;
    // start the timing
    timing(&start_wc_time, &start_cpu_time);
    #pragma omp parallel shared(number_of_cells_inside_mandelbrot_set, total_iterations)
    {
        dsrand_ts(seed);  // Ensure that this is set for each thread

        int number_of_cells_inside_mandelbrot_set_th = 0, total_iterations_th = 0;

        #pragma omp for
        for (size_t i = 0; i <= max_x; ++i)
        {
            double current_bottom_left_x = -2.0 + i * CELL_SIDE_LENGTH;
            double cell_max_x = current_bottom_left_x + CELL_SIDE_LENGTH;

            for (size_t j = 0; j <= max_y; ++j)
            {
                double current_bottom_left_y = 0.0 + j * CELL_SIDE_LENGTH;
                double cell_max_y = current_bottom_left_y + CELL_SIDE_LENGTH;
                double random_x = get_random_double_in_bounds_ts(current_bottom_left_x, cell_max_x);
                double random_y = get_random_double_in_bounds_ts(current_bottom_left_y, cell_max_y);

                int increment = mandelbrot_iteration(random_x, random_y);

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

    int number_of_cells_outside_mandelbrot_set = total_iterations - number_of_cells_inside_mandelbrot_set;
    double area_of_grid = 1.25 * 2.5;
    double ratio = (double)number_of_cells_inside_mandelbrot_set / (number_of_cells_inside_mandelbrot_set + number_of_cells_outside_mandelbrot_set);
    double area = 2.0 * area_of_grid * ratio;

    write_data_to_file("out/omp-2.csv",
                       "omp-ts",
                       num_cores,
                       num_threads,
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