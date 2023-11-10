#ifndef UTILITIES_H
#define UTILITIES_H

int get_environment_value(const char *);

/**
 * write runs to csv for tabulation
 */
void write_data_to_file(const char *filename,
                        const char *program,
                        unsigned seed,
                        float wc_time,
                        float area);

#endif