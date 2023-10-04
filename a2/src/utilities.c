#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

double compute_mandelbrot_area_estimate(int n_i, int total_it)
{
    int n_o = total_it - n_i;
    double area_of_grid = 1.25 * 2.5;
    double ratio = (double)n_i / (n_i + n_o);
    double area = 2.0 * area_of_grid * ratio;

    return area;
}

int get_environment_value(const char *env_name)
{
    char *value_str = getenv(env_name);
    if (value_str == NULL)
    {
        return -1;
    }
    int value = atoi(value_str);

    return value;
}

void write_data_to_file(const char *filename,
                        const char *program,
                        int num_cores,
                        int num_threads,
                        unsigned seed,
                        float wc_time,
                        float area)
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
        fprintf(fp, "num_cores,num_threads,seed,program,wc_time,area\n");
    }

    // add data to row
    fprintf(fp, "%d,%d,%u,%s,%f,%f\n", num_cores, num_threads, seed, program, wc_time, area);
    fclose(fp);
}