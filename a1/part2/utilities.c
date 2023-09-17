#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

void *allocate_double_array(size_t num_elements)
{
    // allocate uninitialized memory
    void *array = malloc(num_elements * sizeof(double));

    // exit if allocation failed
    if (array == NULL)
    {
        fprintf(stderr, "Failed to allocate memory for double array.\n");
        exit(1);
    }

    return array;
}

void initialize_array_with_random_numbers(double *array, size_t num_elements)
{
    // for each element in the array, store a random double in [0,100]
    for (size_t i = 0; i < num_elements; ++i)
    {
        array[i] = ((double)rand() / (double)RAND_MAX) * 100.0;
    }
}

void write_data_to_file(const char *filename, double mega_flops, size_t num_elements)
{
    // check if the file exists
    FILE *check_file = fopen(filename, "r");
    bool does_file_exist = false;
    if (check_file != NULL)
    {
        fclose(check_file);
        does_file_exist = true;
    }

    // open the file with means to append
    FILE *fp = fopen(filename, "a");
    if (fp == NULL)
    {
        fprintf(stderr, "Could not open or create file: %s\n", filename);
        return;
    }

    // if the file does not exist, add the header row
    if (!does_file_exist)
    {
        fprintf(fp, "MFLOPS,N\n");
    }

    // add data to row
    fprintf(fp, "%f,%zu\n", mega_flops, num_elements);
    fclose(fp);
}