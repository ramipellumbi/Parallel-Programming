#include "integrate.h"

double integrate_trapezoidal_rule(double start,
                                  double end,
                                  double num_steps,
                                  double (*func)(double))
{
    double width = (end - start) / num_steps;
    double x;

    double inital_sum = 0.5 * (func(start) + func(end));
    for (int i = 1; i < num_steps; i++)
    {
        x = start + (double)i * width;
        inital_sum += func(x);
    }

    return inital_sum * width;
}

double midpoint_integration_rule(double start,
                                 double end,
                                 double num_steps,
                                 double (*func)(double))
{
    double width_of_interval = 1.0 / num_steps;
    double current_midpoint = start + (1.0 / (2.0 * num_steps));
    double current_function_value = func(current_midpoint);

    double sum = current_function_value * width_of_interval;

    for (int i = 1; i < num_steps; i++)
    {
        current_midpoint += width_of_interval;
        current_function_value = func(current_midpoint);
        sum += current_function_value * width_of_interval;
    }

    return sum;
}