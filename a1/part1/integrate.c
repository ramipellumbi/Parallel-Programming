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