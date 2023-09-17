#include <stdlib.h>
#include <stdio.h>

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

void write_data_to_file(const char *filename, double mega_flops, size_t num_elements)
{
    FILE *fp = fopen(filename, "a");
    if (fp == NULL)
    {
        fprintf(stderr, "Could not open or create file: %s\n", filename);
        return;
    }
    fprintf(fp, "MFLOPS: %f, Elements: %d\n", mega_flops, num_elements);
    fclose(fp);
}