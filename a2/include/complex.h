#ifndef COMPLEX_H
#define COMPLEX_H

/**
 * Adds two complex numbers and stores the result in specified pointers.
 *
 * The function implements the formula:
 * (z1_re + z1_im * i) + (z2_re + z2_im * i) =
 * (z1_re + z2_re) + (z1_im + z2_im) * i
 *
 * @param z1_re Real part of the first complex number.
 * @param z1_im Imaginary part of the first complex number.
 * @param z2_re Real part of the second complex number.
 * @param z2_im Imaginary part of the second complex number.
 * @param result_re Pointer to store the real part of the result.
 * @param result_im Pointer to store the imaginary part of the result.
 */
static inline void complex_add(double z1_re, double z1_im,
                               double z2_re, double z2_im,
                               double *result_re, double *result_im)
{
    // Calculate the real part of the result: z1_re + z2_re
    *result_re = z1_re + z2_re;

    // Calculate the imaginary part of the result: z1_im + z2_im
    *result_im = z1_im + z2_im;
}

/**
 * Calculates the squared magnitude of a complex number.
 *
 * The squared magnitude is computed as `(z_re * z_re) + (z_im * z_im)`,
 *
 * @param z_re Real part of the complex number.
 * @param z_im Imaginary part of the complex number.
 * @return The squared magnitude of the complex number.
 */
static inline double complex_magnitude_squared(double z_re, double z_im)
{
    // Multiply the real part by itself: z_re * z_re
    // Multiply the imaginary part by itself: z_im * z_im
    // Sum the two products to get the squared magnitude
    return z_re * z_re + z_im * z_im;
}

/**
 * Multiplies two complex numbers and stores the result in specified pointers.
 *
 * The function implements the formula:
 * (z1_re + z1_im * i) * (z2_re + z2_im * i) =
 * (z1_re * z2_re - z1_im * z2_im) + (z1_im * z2_re + z1_re * z2_im) * i
 *
 * @param z1_re Real part of the first complex number.
 * @param z1_im Imaginary part of the first complex number.
 * @param z2_re Real part of the second complex number.
 * @param z2_im Imaginary part of the second complex number.
 * @param result_re Pointer to store the real part of the result.
 * @param result_im Pointer to store the imaginary part of the result.
 */
static inline void complex_multiply(double z1_re, double z1_im,
                                    double z2_re, double z2_im,
                                    double *result_re, double *result_im)
{
    // Calculate the real part of the result: z1_re * z2_re - z1_im * z2_im
    *result_re = z1_re * z2_re - z1_im * z2_im;

    // Calculate the imaginary part of the result: z1_im * z2_re + z1_re * z2_im
    *result_im = z1_im * z2_re + z1_re * z2_im;
}

#endif