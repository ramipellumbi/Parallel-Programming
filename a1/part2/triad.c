#include <math.h>
#include <stdlib.h>
#include <timing.h>
#include <time.h>

void *allocate_double_array(size_t num_elements)
{
    void *array = malloc(num_elements * sizeof(double));
    if (array == NULL)
    {
        fprintf(stderr, "Failed to allocate memory for double array.\n");
        exit(1);
    }

    return array;
}

void initialize_array_with_random_numbers(double *array, size_t num_elements)
{
    for (size_t i = 0; i < num_elements; ++i)
    {
        array[i] = ((double)rand() / (double)RAND_MAX) * 100.0;
    }
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s <k-value>\n", argv[0]);
        return 1;
    }

    int k = atoi(argv[1]);
    if (k < 0 || k > 25)
    {
        fprintf(stderr, "Using %d which is not in range 1...25", k);
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

    int repeat = 1;
    double runtime = 0.0;
    while (runtime < 1.0)
    {
        // start the timing
        timing(&start_wc_time, &start_cpu_time);
        for (int r = 0; r < repeat; ++r)
        {
            for (size_t i = 0; i < num_elements; i++)
            {
                a[i] = b[i] + c[i] * d[i];
            }
            if (a[num_elements >> 1] < 0)
            {
                printf("running dummy");
                dummy(a, b, c, d);
            }
        }
        // end the timing
        timing(&end_wc_time, &end_cpu_time);
        runtime = end_wc_time - start_wc_time;
        repeat *= 2;
    }
    repeat /= 2;

    printf("\nelapsed wall clock time = %f\n", end_wc_time - start_wc_time);
    printf("elapsed cpu time = %f\n", end_cpu_time - start_cpu_time);
    free(a);
    free(b);
    free(c);
    free(d);

    return 0;
}