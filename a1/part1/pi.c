/**
 * Benchmark a code that approximates \pi by numerically integrating the function
 * f(x) = 1.0 / (1.0 + x^2) from 0 to 1 and multiplying the result by 4
 *
 * All floating point operations are performed in double precision
 */

#include "integrate.h"
#include <stdio.h>
#include "timing.h"

static const double NUM_STEPS = 1000000;
static const double START_X = 0.0;
static const double END_X = 1.0;
static const double SCALING_FACTOR = 4.0;

double function_to_integrate(double x);

double function_to_integrate(double x)
{
    // 1 divide, 1 add, 1 multiply
    return 1.0 / (1.0 + x * x);
}

int main(int argc, char *argv[])
{
    // initialize timing measures
    double start_wc_time = 0.0, end_wc_time = 0.0;
    double start_cpu_time = 0.0, end_cpu_time = 0.0;

    // start the timing
    timing(&start_wc_time, &start_cpu_time);
    double result = integrate_midpoint_rule(START_X,
                                            END_X,
                                            NUM_STEPS,
                                            function_to_integrate);
    timing(&end_wc_time, &end_cpu_time);
    double pi = result * SCALING_FACTOR;
    // end the timing

    printf("\npi = %f\n", pi);
    printf("elapsed wall clock time = %f\n", end_wc_time - start_wc_time);
    printf("elapsed cpu time = %f\n", end_cpu_time - start_cpu_time);
    return 0;
}