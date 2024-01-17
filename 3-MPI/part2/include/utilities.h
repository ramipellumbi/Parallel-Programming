#ifndef UTILITIES_H
#define UTILITIES_H

int get_environment_value(const char *env_name);

/**
 * write runs to csv for tabulation
 */
void write_data_to_file(const char *filename,
                        const char *program,
                        int N,
                        int np,
                        double exe_time,
                        double f_norm);

#endif