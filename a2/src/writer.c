#include <stdbool.h>
#include <stdio.h>
#include "writer.h"

void write_data_to_file(const char *filename,
                        const char *program,
                        int num_cores,
                        int num_threads,
                        unsigned seed,
                        double wc_time,
                        double area)
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