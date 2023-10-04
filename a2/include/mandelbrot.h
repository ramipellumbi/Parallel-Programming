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

    // stores the current magnitude of z through the iterations; starts at 0
    double magnitude_squared = 0.0;

    // Maximum number of iterations to decide whether the number is in the Mandelbrot set or not.
    // If the magnitude of z goes above 2 and stays there, it's not in the set.
    size_t MAX_ITERATIONS = 25000;
    size_t current_iteration = 0;

    // Perform MAX_ITERATIONS - 1 + 1 iterations
    while (current_iteration < MAX_ITERATIONS)
    {
        double temp_re = z_re * z_re - z_im * z_im + c_re;
        z_im = (z_re + z_re) * z_im + c_im;
        z_re = temp_re;
        magnitude_squared = z_re * z_re + z_im * z_im;

        if (magnitude_squared > 4.0)
        {
            return 0;
        }

        current_iteration += 1;
    }

    return 1;
}

#endif