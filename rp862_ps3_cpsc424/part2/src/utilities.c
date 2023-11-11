#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

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
                        int N,
                        int np,
                        double exe_time,
                        double f_norm)
{
    int num_cores = get_environment_value("SLURM_CPUS_PER_TASK");
    if (num_cores == -1)
    {
        fprintf(stderr, "SLURM_CPUS_PER_TASK not set");
    }

    int num_tasks_per_node = get_environment_value("SLURM_NTASKS_PER_NODE");
    if (num_tasks_per_node == -1)
    {
        fprintf(stderr, "SLURM_NTASKS_PER_NODE not set");
    }

    int num_tasks_per_socket = get_environment_value("SLURM_NTASKS_PER_SOCKET");
    if (num_tasks_per_socket == -1)
    {
        fprintf(stderr, "SLURM_NTASKS_PER_SOCKET not set");
    }

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
        fprintf(fp, "program,num_cores,np,num_tasks_per_node,num_tasks_per_socket,N,exe_time,f_norm\n");
    }

    // add data to row
    fprintf(fp, "\"%s\",%d,%d,%d,,%d,%d,%f,%.12f\n", program, num_cores, np, num_tasks_per_node, num_tasks_per_socket, N, exe_time, f_norm);
    fclose(fp);
}