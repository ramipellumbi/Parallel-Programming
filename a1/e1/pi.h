#ifndef PI_H
#define PI_H

/**
 * Function to be integrated for the approximation.
*/
double function_to_integrate(double x);

/**
 * Compute integral with trapezoidal rule: https://en.wikipedia.org/wiki/Trapezoidal_rule 
 * 
 * Note: the domain is discretized into N equally spaced panels and thus we have a uniform grid.
*/
double integrate_trapezoidal_rule(double start, double end, double num_steps, double (*func)(double));

#endif 
