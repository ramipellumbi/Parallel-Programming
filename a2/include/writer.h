#ifndef WRITER_H
#define WRITER_H

/**
 * write runs to csv for tabulation
 */
void write_data_to_file(const char *filename,
                        const char *program,
                        int num_cores,
                        int num_threads,
                        unsigned seed,
                        double wc_time,
                        double area);

#endif
