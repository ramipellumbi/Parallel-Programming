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

// 3 flops
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
    // Performing 5 FLOPS per iteration, and NUM_STEPS iterations.
    double result = integrate_midpoint_rule(START_X,
                                            END_X,
                                            NUM_STEPS,
                                            function_to_integrate);
    // end the timing
    timing(&end_wc_time, &end_cpu_time);
    double pi = result * SCALING_FACTOR;

    double elapsed_wc_time = end_wc_time - start_wc_time;
    double elapsed_cpu_time = end_cpu_time - start_cpu_time;

    // Have CYCLE_TIME = 1 / PROCESSOR_FREQUENCY
    // Have NUMBER_OF_CYCLES = WALL_CLOCK_TIME / CYCLE_TIME = WALL_CLOCK_TIME * PROCESSOR_FREQUENCY
    // Assuming that division dominates the computation time of each iteration,
    // we can estimate the number of cycles spent on division by dividing the
    // total number of cycles by the total number of flops.

    double total_number_of_flops = 5 * NUM_STEPS;
    double cpu_frequency_hertz = 3.7e9;
    double cycle_time = 1.0 / cpu_frequency_hertz;
    double total_number_of_cycles = elapsed_wc_time * cpu_frequency_hertz;

    double mega_flops_per_second = total_number_of_flops / elapsed_wc_time / 1.0e6;
    double estimated_divide_latency = total_number_of_cycles / NUM_STEPS;

    printf("\npi estimate = %f\n", pi);
    printf("elapsed wall clock time = %f\n", elapsed_wc_time);
    printf("elapsed cpu time = %f\n", elapsed_cpu_time);
    printf("Cycle time = %.17g\n", cycle_time);
    printf("Total number of cycles = %f\n", total_number_of_cycles);
    printf("Estimated divide latency = %f\n", estimated_divide_latency);
    return 0;
}