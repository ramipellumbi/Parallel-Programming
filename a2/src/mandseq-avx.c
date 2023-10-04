#include <immintrin.h>
#include <nmmintrin.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include "complex.h"
#include "drand.h"
#include "mandelbrot.h"
#include "timing.h"
#include "utilities.h"

static const NUM_X_ITERATIONS = 2500;
static const NUM_Y_ITERATIONS = 1250;
static const double CELL_SIDE_LENGTH = 0.001;
static const int PACKING_SIZE = 8;

int sum_of_bits_in_mmask16(__mmask8 mask)
{
    int mask_as_int = (int)mask;
    int count = 0;
    while (mask_as_int)
    {
        count += mask_as_int & 1;
        mask_as_int >>= 1;
    }
    return count;
}

int main(int argc, char *argv[])
{
    int num_cores = get_environment_value("SLURM_CPUS_PER_TASK");
    if (num_cores == -1)
    {
        fprintf(stderr, "SLURM_CPUS_PER_TASK not set");
        return 0;
    }

    // set seed
    unsigned seed = 144545;
    dsrand(seed);

    // initialize a random number for each cell immediately
    // each cell needs 2 random numbers so we have 2 * num_cells in the packing loop
    // precomputed numbers
    int num_elements = 2500 * (1250 / PACKING_SIZE * PACKING_SIZE) * 2;
    double precomputed_randoms[num_elements];
    for (int i = 0; i < num_elements; ++i)
    {
        precomputed_randoms[i] = drand();
    }

    // initialize timing measures
    double start_wc_time = 0.0, end_wc_time = 0.0;
    double start_cpu_time = 0.0, end_cpu_time = 0.0;

    int number_of_cells_inside_mandelbrot_set = 0, total_iterations = 0;

    // start the timing
    timing(&start_wc_time, &start_cpu_time);

    // index to track where we are in the random index array
    int random_index = 0;

    // We have 8 packed doubles in a 512 bit register
    // This stores the corner offsets - 0, 0.001, 0.002, ..., 0.007
    // so that at every x value we can assess 8 y values at a time
    __m512d pxs_deltas512 = _mm512_mul_pd(
        _mm512_set_pd(0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0),
        _mm512_set1_pd(CELL_SIDE_LENGTH));

    // for each x value
    for (size_t n = 0; n < 2500; n++)
    {
        double current_bottom_left_x = -2.0 + CELL_SIDE_LENGTH * n;
        double max_x = current_bottom_left_x + CELL_SIDE_LENGTH;

        // let's compute 8 y-values simultaneously
        for (size_t m = 0; m < 1250 / PACKING_SIZE * PACKING_SIZE; m += PACKING_SIZE)
        {
            // grab 8 random numbers for the 8 needed random x coordinates
            __m512d random_numbers_x = _mm512_loadu_pd(&precomputed_randoms[random_index]);
            // grab 8 random numbers for the 8 needed random y coordinates
            __m512d random_numbers_y = _mm512_loadu_pd(&precomputed_randoms[random_index + PACKING_SIZE]);

            // get the 8 bottom left corners for this iteration of the loop
            __m512d bottom_left_y_values = _mm512_add_pd(
                _mm512_set1_pd(0.0),
                _mm512_add_pd(_mm512_set1_pd(m * CELL_SIDE_LENGTH), pxs_deltas512));
            // get the 8 top left corners for this iteration of the loop
            __m512d top_left_y_values = _mm512_add_pd(
                _mm512_set1_pd(CELL_SIDE_LENGTH),
                bottom_left_y_values);

            // compute the random x coordinates for the 8 cells
            __m512d x_values = _mm512_fmadd_pd(
                random_numbers_x,
                _mm512_set1_pd(max_x - current_bottom_left_x),
                _mm512_set1_pd(current_bottom_left_x));

            // compute the random y coordinates for the 8 cells
            __m512d y_values = _mm512_fmadd_pd(
                random_numbers_y,
                _mm512_sub_pd(top_left_y_values, bottom_left_y_values),
                bottom_left_y_values);

            // stores the number of iterations each c went through
            __m512i iters = _mm512_setzero_si512();

            // These are the 8 z values
            __m512d z_re = _mm512_set1_pd(0.0);
            __m512d z_im = _mm512_set1_pd(0.0);

            // 8 bit mask. 1 means that the c value at that index has diverged.
            __mmask8 diverged_indices = 0;

            // Assess the 8 c values concurrently
            for (size_t iteration = 0; iteration < 25000; iteration++)
            {
                // componentwise: z_re = z_re * z_re + z_im * z_im + c_re
                __m512d xsn = _mm512_add_pd(_mm512_sub_pd(
                                                _mm512_mul_pd(z_re, z_re),
                                                _mm512_mul_pd(z_im, z_im)),
                                            x_values);

                // componentwise: z_im = 2 * z_re * z_im + c_im
                __m512d ysn = _mm512_add_pd(_mm512_mul_pd(
                                                _mm512_mul_pd(
                                                    _mm512_set1_pd(2.0f),
                                                    z_re),
                                                z_im),
                                            y_values);

                // we update only for the indices that have not diverged
                __m512d z_re = _mm512_mask_mov_pd(z_re, ~diverged_indices, xsn);
                __m512d z_im = _mm512_mask_mov_pd(z_im, ~diverged_indices, ysn);

                // compute the magnitude squared componentwise (could optimize
                // this to only do so for non diverged indices)
                __m512d magnitude_squared = _mm512_add_pd(
                    _mm512_mul_pd(z_re, z_re),
                    _mm512_mul_pd(z_im, z_im));

                // Generate a mask for numbers that have diverged (magnitude squared > 4)
                __mmask16 maskDiverged = _mm512_cmp_pd_mask(magnitude_squared,
                                                            _mm512_set1_pd(4.0),
                                                            _CMP_GT_OS);

                // Update diverged_indices using bitwise OR operation
                diverged_indices |= maskDiverged;

                // Update the iteration counter, incrementing only where diverged_indices is zero
                iters = _mm512_mask_add_epi32(iters,
                                              ~diverged_indices,
                                              iters,
                                              _mm512_set1_epi32(1));

                // Break if all indices have diverged
                if (diverged_indices == 0xFF)
                {
                    break;
                }
            }

            // the 1's in this mask are the iterations that did NOT diverge
            __mmask16 maskNeverDiverged = ~diverged_indices;
            int count = sum_of_bits_in_mmask16(maskNeverDiverged);
            number_of_cells_inside_mandelbrot_set += count;

            random_index += PACKING_SIZE * 2;
            total_iterations += PACKING_SIZE;
        }

        // cleanup by performing the naive serial implementation
        for (size_t m = 1250 / PACKING_SIZE * PACKING_SIZE; m < 1250; m += 1)
        {
            double current_bottom_left_y = 0.0 + CELL_SIDE_LENGTH * m;
            double max_y = current_bottom_left_y + CELL_SIDE_LENGTH;

            double c_re = current_bottom_left_x + drand() * (max_x - current_bottom_left_x);
            double c_im = current_bottom_left_y + drand() * (max_y - current_bottom_left_y);

            double z_re = 0.0;
            double z_im = 0.0;

            int counter = 1;
            for (size_t iteration = 0; iteration < 25000; iteration++)
            {
                double xn = z_re * z_re - z_im * z_im + c_re;
                z_im = (z_re + z_re) * z_im + c_im;
                z_re = xn;
                if (z_re * z_re + z_im * z_im > 4.0)
                {
                    counter = 0;
                    break;
                }
            }

            number_of_cells_inside_mandelbrot_set += counter;
            total_iterations++;
        }
    }
    // end the timing
    timing(&end_wc_time, &end_cpu_time);

    double elapsed_wc_time = end_wc_time - start_wc_time;
    double elapsed_cpu_time = end_cpu_time - start_cpu_time;
    printf("\nTotal iterations: %d", total_iterations);

    int number_of_cells_outside_mandelbrot_set = total_iterations - number_of_cells_inside_mandelbrot_set;
    double area_of_grid = 1.25 * 2.5;
    double ratio = (double)number_of_cells_inside_mandelbrot_set / (number_of_cells_inside_mandelbrot_set + number_of_cells_outside_mandelbrot_set);
    double area = 2.0 * area_of_grid * ratio;

    write_data_to_file("out/serial.csv",
                       "serial-avx",
                       num_cores,
                       1,
                       seed,
                       elapsed_wc_time,
                       area);

    printf("\narea estimate = %f\n", area);
    printf("Num inside: %d\n", number_of_cells_inside_mandelbrot_set);
    printf("Num outside: %d\n", number_of_cells_outside_mandelbrot_set);
    printf("elapsed wall clock time = %f\n", elapsed_wc_time);
    printf("elapsed cpu time = %f\n", elapsed_cpu_time);

    return 0;
}