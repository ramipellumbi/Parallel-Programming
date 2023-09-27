#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#include "cell.h"
#include "drand.h"
#include "point.h"
#include "rand_utils.h"
#include "timing.h"

// the complex plane is isomorphic to the plane.
bool mandelbrot_iteration(Point c, size_t MAX_ITERATIONS)
{
    Point z = (Point){0.0, 0.0};
    size_t i;

    while (i < MAX_ITERATIONS)
    {
        double magnitude = magnitude_of_point_squared(z);
        if (magnitude > 4.0)
        {
            return false;
        }
        z = add_points(multiply_points(z, z), c);
        i += 1;
    }

    return true;
}

int main(int argc, char *argv[])
{
    // initialize timing measures
    double start_wc_time = 0.0, end_wc_time = 0.0;
    double start_cpu_time = 0.0, end_cpu_time = 0.0;

    dsrand(12345);
    Cell cell;
    double current_bottom_left_x = -2.0;
    double current_bottom_left_y = 0.0;

    int number_of_cells_inside_mandelbrot_set = 0;
    int number_of_cells_outside_mandelbrot_set = 0;

    // start the timing
    timing(&start_wc_time, &start_cpu_time);

    // at each bottom x coordinate
    while (current_bottom_left_x <= 0.5)
    {
        // go through all the y boxes at this x
        while (current_bottom_left_y <= 1.25)
        {
            create_cell(&cell, current_bottom_left_x, current_bottom_left_y);
            Point c = generate_random_point_in_cell(&cell);

            bool is_in_set = mandelbrot_iteration(c, 25000);
            if (is_in_set)
            {
                number_of_cells_inside_mandelbrot_set += 1;
            }
            else
            {
                number_of_cells_outside_mandelbrot_set += 1;
            }

            current_bottom_left_y += CELL_SIDE_LENGTH;
        }

        // go to next box shifted CELL_SIDE_LENGTH to the right
        current_bottom_left_x += CELL_SIDE_LENGTH;
        // reset y
        current_bottom_left_y = 0.0;
    }
    // end the timing
    timing(&end_wc_time, &end_cpu_time);

    double elapsed_wc_time = end_wc_time - start_wc_time;
    double elapsed_cpu_time = end_cpu_time - start_cpu_time;

    double area_of_grid = 1.25 * 2.5;
    double area = 2.0 * area_of_grid * (number_of_cells_inside_mandelbrot_set / (number_of_cells_inside_mandelbrot_set + number_of_cells_outside_mandelbrot_set));

    printf("\narea estimate = %f\n", area);
    printf("Random num %f\n", drand());
    printf("Num inside: %d\n", number_of_cells_inside_mandelbrot_set);
    printf("Num outside: %d\n", number_of_cells_outside_mandelbrot_set);
    printf("elapsed wall clock time = %f\n", elapsed_wc_time);
    printf("elapsed cpu time = %f\n", elapsed_cpu_time);

    return 0;
}