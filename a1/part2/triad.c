#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <timing.h>
#include <time.h>
#include "utilities.h"

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s <k-value>\n", argv[0]);
        return 1;
    }

    int k = atoi(argv[1]);
    if (k < 0)
    {
        fprintf(stderr, "Using %d which is not a positive integer", k);
        return 1;
    }

    // seed random number generator
    srand(time(NULL));

    // initialize timing measures
    double start_wc_time = 0.0, end_wc_time = 0.0;
    double start_cpu_time = 0.0, end_cpu_time = 0.0;

    size_t num_elements = floor(pow(2.1, k));

    double *a = (double *)allocate_double_array(num_elements);
    double *b = (double *)allocate_double_array(num_elements);
    double *c = (double *)allocate_double_array(num_elements);
    double *d = (double *)allocate_double_array(num_elements);

    initialize_array_with_random_numbers(a, num_elements);
    initialize_array_with_random_numbers(b, num_elements);
    initialize_array_with_random_numbers(c, num_elements);
    initialize_array_with_random_numbers(d, num_elements);

    int number_of_repetitions = 1;
    double runtime = 0.0;
    while (runtime < 1.0)
    {
        // start the timing
        timing(&start_wc_time, &start_cpu_time);
        for (int r = 0; r < number_of_repetitions; ++r)
        {
            for (size_t i = 0; i < num_elements; i++)
            {
                a[i] = b[i] + c[i] * d[i];
            }

            if (a[num_elements >> 1] < 0)
            {
                dummy(a, b, c, d);
            }
        }
        // end the timing
        timing(&end_wc_time, &end_cpu_time);
        runtime = end_wc_time - start_wc_time;
        number_of_repetitions *= 2;
    }
    number_of_repetitions /= 2;

    double elapsed_wc_time = end_wc_time - start_wc_time;
    double elapsed_cpu_time = end_cpu_time - start_cpu_time;

    double total_mega_flops = (2.0 * (double)num_elements * number_of_repetitions) / elapsed_wc_time / 1.0e6;

    printf("\nnumber of array elements: %d", num_elements);
    printf("\n%d repetitions performed", number_of_repetitions);
    printf("\nelapsed wall clock time = %f\n", elapsed_wc_time);
    printf("elapsed cpu time = %f\n", elapsed_cpu_time);
    printf("estimated MFLOPS: %f", total_mega_flops);
    write_data_to_file("data/part_2_data.csv", total_mega_flops, num_elements);

    free(a);
    free(b);
    free(c);
    free(d);

    return 0;
}