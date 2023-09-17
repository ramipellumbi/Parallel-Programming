#include "integrate.h"

double integrate_midpoint_rule(double start,
                               double end,
                               double num_steps,
                               double (*func)(double))
{
    // width of the interval at each midpoint
    double width_of_interval = 1.0 / num_steps;
    // initial midpoint is halfway between start and 1/N
    double current_midpoint = (start + (1.0 / num_steps)) / 2.0;
    // function value at the midpoint - initialized to the value at inital midpoint
    double current_function_value = func(current_midpoint);

    // initialize the sum value to the first rectangles area
    double sum = current_function_value * width_of_interval;

    // compute each midpoint x_1, ... , x_{N-1}, function value at that midpoint,
    // and add the rectangle area to the sum
    for (int i = 1; i < num_steps; i++)
    {
        current_midpoint += width_of_interval;
        current_function_value = func(current_midpoint);
        sum += current_function_value * width_of_interval;
    }

    return sum;
}