#include <complex.h>    
#include <stdio.h>

#include "drand.h"
#include "timing.h"

static const double CELL_SIDE_LENGTH = 0.001;

// | a + bi |^2 = a^2 + b^
static inline double complex_magnitude_squared(double a, double b) {
    return a * a + b * b;
}

// (a + bi) * (c + di) = (ac - bd) + (bc + ad)i
// assign (ac - bd) to re and (bc + ad) to im
static inline void complex_multiply(double a, double b, double c, double d, double *re, double *im) {
    *re = a * c - b * d;
    *im = b * c + a * d;
}

// (a + bi) + (c + di) = (a + c) + (b + d)i
// assign (a + c) to re and (b + d) to im
static inline void complex_add(double a, double b, double c, double d, double *re, double *im) {
    *re = a + c;
    *im = b + d;
}

static inline int mandelbrot_iteration(double x, double y)
{
    // z is updated via the Mandelbrot update rule z <- z^2 + c, where c = x + iy
    double z_x = 0.0, z_y = 0.0;
    // used to allow for the invariant that z := z_x + i z_y is always the result of a valid Mandelbrot update for c := x + iy
    double temp_x, temp_y;

    for (size_t i = 0; i < 25000; i += 5)
    {
        for (size_t i = 0; i < 5; ++i)
        {
            // compute z^2 and store the result of the complex multiply in temp_x, temp_y
            complex_multiply(z_x, z_y, z_x, z_y, &temp_x, &temp_y);
            // compute z^2 + c and store the result back in z_x, z_y
            complex_add(temp_x, temp_y, x, y, &z_x, &z_y);
        }

        if (complex_magnitude_squared(z_x, z_y) > 4.0)
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

    for(double current_bottom_left_x = -2.0; current_bottom_left_x <= 0.5; current_bottom_left_x += CELL_SIDE_LENGTH)
    {
        double max_x = current_bottom_left_x + CELL_SIDE_LENGTH;
        for(double current_bottom_left_y = 0.0; current_bottom_left_y <= 1.25; current_bottom_left_y += CELL_SIDE_LENGTH)
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