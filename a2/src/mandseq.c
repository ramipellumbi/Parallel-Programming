#include <stdio.h>

#include "complex.h"
#include "drand.h"
#include "timing.h"

static const double CELL_SIDE_LENGTH = 0.001;

/**
 * @param
 * @param
 * @returns 1 if c is in the mandelbrot set; else 0
 */
static inline int mandelbrot_iteration(double c_re, double c_im)
{
    // z := z_re + i * z_im is updated via the Mandelbrot update rule:
    // z <- z^2 + c, where c := x + iy
    double z_re = 0.0, z_im = 0.0;

    // used to permit the invariant that z is the result of a Mandelbrot update
    double temp_re, temp_im;

    size_t MAX_ITERATIONS = 25000;
    size_t UNROLL_COUNT = 10;

    for (size_t i = 0; i < MAX_ITERATIONS / UNROLL_COUNT; i += 5)
    {
        for (size_t i = 0; i < UNROLL_COUNT; ++i)
        {
            // compute z^2 and store the result of the complex multiply in temp_re, temp_im
            complex_multiply(z_re, z_im, z_re, z_im, &temp_re, &temp_im);
            // compute z^2 + c and store the result back in z_re,z_im
            complex_add(temp_re, temp_im, c_re, c_im, &z_re, &z_im);
        }

        if (complex_magnitude_squared(z_re, z_im) > 4.0)
        {
            return 0;
        }
    }
    return 1;
}

static inline double get_random_double_in_bounds(double min, double max)
{
    return min + drand() * (max - min);
}

int main(int argc, char *argv[])
{
    // initialize timing measures
    double start_wc_time = 0.0, end_wc_time = 0.0;
    double start_cpu_time = 0.0, end_cpu_time = 0.0;

    dsrand(10);

    int number_of_cells_inside_mandelbrot_set = 0, total_iterations = 0;

    // start the timing
    timing(&start_wc_time, &start_cpu_time);

    for (double current_bottom_left_x = -2.0; current_bottom_left_x <= 0.5; current_bottom_left_x += CELL_SIDE_LENGTH)
    {
        double max_x = current_bottom_left_x + CELL_SIDE_LENGTH;
        for (double current_bottom_left_y = 0.0; current_bottom_left_y <= 1.25; current_bottom_left_y += CELL_SIDE_LENGTH)
        {
            double max_y = current_bottom_left_y + CELL_SIDE_LENGTH;
            double random_x = get_random_double_in_bounds(current_bottom_left_x, max_x);
            double random_y = get_random_double_in_bounds(current_bottom_left_y, max_y);

            int increment = mandelbrot_iteration(random_x, random_y);

            number_of_cells_inside_mandelbrot_set += increment;
            total_iterations += 1;
        }
    }
    // end the timing
    timing(&end_wc_time, &end_cpu_time);

    double elapsed_wc_time = end_wc_time - start_wc_time;
    double elapsed_cpu_time = end_cpu_time - start_cpu_time;

    double number_of_cells_outside_mandelbrot_set = total_iterations - number_of_cells_inside_mandelbrot_set;
    double area_of_grid = 1.25 * 2.5;
    double ratio = (double)number_of_cells_inside_mandelbrot_set / (number_of_cells_inside_mandelbrot_set + number_of_cells_outside_mandelbrot_set);
    double area = 2.0 * area_of_grid * ratio;

    printf("\narea estimate = %f\n", area);
    printf("Num inside: %d\n", number_of_cells_inside_mandelbrot_set);
    printf("Num outside: %d\n", number_of_cells_outside_mandelbrot_set);
    printf("elapsed wall clock time = %f\n", elapsed_wc_time);
    printf("elapsed cpu time = %f\n", elapsed_cpu_time);

    return 0;
}