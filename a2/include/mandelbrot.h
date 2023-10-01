#ifndef MANDELBROT_H
#define MANDELBROT_H

/**
 * Determines if a complex number c = c_re + i*c_im belongs to the Mandelbrot set.
 * The function uses loop unrolling for performance optimization.
 *
 * @param c_re The real part of the complex number c.
 * @param c_im The imaginary part of the complex number c.
 * @returns Returns 1 if c is in the Mandelbrot set, otherwise 0.
 */
static inline int mandelbrot_iteration(double c_re, double c_im)
{
    // Initialize z to 0. z is a complex number represented by z_re and z_im.
    // The Mandelbrot sequence is generated based on the iterative formula: z = z^2 + c
    double z_re = 0.0, z_im = 0.0;

    // Temporary variables to hold intermediate results of complex operations.
    double temp_re, temp_im;

    // Maximum number of iterations to decide whether the number is in the Mandelbrot set or not.
    // If the magnitude of z goes above 2 and stays there, it's not in the set.
    size_t MAX_ITERATIONS = 25000;

    // The count of iterations to unroll in the inner loop for performance optimization.
    size_t UNROLL_COUNT = 10;

    // Outer loop: Process batches of UNROLL_COUNT iterations at a time
    for (size_t i = 0; i < MAX_ITERATIONS / UNROLL_COUNT; i += 1)
    {
        // This is a loop unrolling optimization to minimize the loop overhead.
        for (size_t i = 0; i < UNROLL_COUNT; ++i)
        {
            // Calculate z^2. Store result in temporary variables temp_re and temp_im
            complex_multiply(z_re, z_im, z_re, z_im, &temp_re, &temp_im);

            // Update z by adding c to z^2. Store result in z_re and z_im
            complex_add(temp_re, temp_im, c_re, c_im, &z_re, &z_im);
        }

        // Check if the magnitude squared of the current z value squared exceeds 4.
        if (complex_magnitude_squared(z_re, z_im) > 4.0)
        {
            // If it does, the number does not belong to the Mandelbrot set and we return 0.
            return 0;
        }
    }

    // If we've reached this point, the number is assumed to be in the Mandelbrot set.
    return 1;
}

#endif