#include <immintrin.h>
#include <nmmintrin.h>

#ifndef MANDELBROT_H
#define MANDELBROT_H

/**
 * Determines if a complex number c = c_re + i*c_im belongs to the Mandelbrot set.
 *
 * @param c_re The real part of the complex number c.
 * @param c_im The imaginary part of the complex number c.
 * @param max_iterations Maximum number of iterations to decide whether the number is in the Mandelbrot set or not.
 * @returns Returns 1 if c is in the Mandelbrot set, otherwise 0.
 */
static inline int mandelbrot_iteration(double c_re, double c_im, size_t max_iterations)
{
    // Initialize z to 0. z is a complex number represented by z_re and z_im.
    // The Mandelbrot sequence is generated based on the iterative formula: z = z^2 + c
    double z_re = 0.0, z_im = 0.0;

    // stores the current magnitude of z through the iterations; starts at 0
    double magnitude_squared = 0.0;

    for (size_t i = 0; i < max_iterations; i++)
    {
        double temp_re = z_re * z_re - z_im * z_im + c_re;
        z_im = (z_re + z_re) * z_im + c_im;
        z_re = temp_re;
        magnitude_squared = z_re * z_re + z_im * z_im;

        if (magnitude_squared > 4.0)
        {
            return 0;
        }
    }

    return 1;
}

/**
 * Determines if a set of complex number belongs to the Mandelbrot set.
 *
 * Each __m512d contains 8 double precision floating point numbers.
 *
 * @param x_values The real parts of the complex numbers.
 * @param c_im The imaginary parts of the complex numbers.
 * @param max_iterations Maximum number of iterations to decide whether the number is in the Mandelbrot set or not.
 * @returns Returns an 8 bit mask which has a 1 on the indices that diverged
 */
static inline __mmask8 mandelbrot_iteration_avx(__m512d x_values,
                                                __m512d y_values,
                                                size_t max_iterations)
{
    // These are the 8 z values
    __m512d z_re = _mm512_set1_pd(0.0);
    __m512d z_im = _mm512_set1_pd(0.0);

    // 8 bit mask. 1 means that the c value at that index has diverged.
    __mmask8 diverged_indices = 0;

    // Assess the 8 c values concurrently
    for (size_t iteration = 0; iteration < max_iterations; iteration++)
    {
        // componentwise: z_re = z_re * z_re + z_im * z_im + c_re
        __m512d xsn = _mm512_add_pd(_mm512_sub_pd(
                                        _mm512_mul_pd(z_re, z_re),
                                        _mm512_mul_pd(z_im, z_im)),
                                    x_values);

        // componentwise: z_im = 2 * z_re * z_im + c_im
        __m512d ysn = _mm512_add_pd(_mm512_mul_pd(
                                        _mm512_add_pd(
                                            z_re,
                                            z_re),
                                        z_im),
                                    y_values);

        // Update only those positions where diverged_indices is zero
        __m512d new_z_re = _mm512_mask_mov_pd(z_re, ~diverged_indices, xsn);
        __m512d new_z_im = _mm512_mask_mov_pd(z_im, ~diverged_indices, ysn);

        // Use the new values
        // I SPENT HOURS DEBUGGING THIS YOU CAN NOT
        // JUST SET z_re, z_im WHERE YOU SET new_z_re, new_z_im
        z_re = new_z_re;
        z_im = new_z_im;

        // compute the magnitude squared componentwise (could optimize
        // this to only do so for non diverged indices)
        __m512d magnitude_squared = _mm512_add_pd(
            _mm512_mul_pd(z_re, z_re),
            _mm512_mul_pd(z_im, z_im));

        // Generate a mask for numbers that have diverged (magnitude squared > 4)
        __mmask8 maskDiverged = _mm512_cmp_pd_mask(magnitude_squared,
                                                   _mm512_set1_pd(4.0),
                                                   _CMP_GT_OS);

        // Update diverged_indices using bitwise OR operation
        diverged_indices |= maskDiverged;

        // Break if all indices have diverged
        if (diverged_indices == 0xFF)
        {
            return diverged_indices;
        }
    }

    return diverged_indices;
}

#endif