#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdlib.h>
#include <stdio.h>

/**
 * Dynamically allocates memory for an array of `double`'s of the 
 * specified size. It will exit the program if the memory allocation fails.
 *
 * @param num_elements The number of elements in the array.
 * @return A pointer to the allocated array.
 */
void *allocate_double_array(size_t num_elements);

/**
 * Fills the provided array with random numbers. Each number is generated 
 * to be a double-precision floating-point value between 0 and 100.
 *
 * @param array The array to fill with random numbers.
 * @param num_elements The number of elements in the array.
 */
void initialize_array_with_random_numbers(double *array, size_t num_elements);

/**
 * This function appends the MegaFLOPS and the number of elements to the specified 
 * text file. Each new record is written on a new line. If the file does not exist, 
 * it will be created. If the file cannot be opened or created, an error message 
 * will be printed to stderr.
 *
 * @param filename The name of the file to write to.
 * @param mega_flops The MFLOPS to be written.
 * @param num_elements The number of elements to be written.
 */
void write_data_to_file(const char *filename, double mega_flops, size_t num_elements);

#endif 
