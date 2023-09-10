#ifndef INTEGRATE_H
#define INTEGRATE_H

/**
 * Numerically integrate a function using the trapezoidal rule: https://en.wikipedia.org/wiki/Trapezoidal_rule.
 *
 * @param start The start of the interval to integrate over
 * @param end The end of the interval to integrate over
 * @param num_steps The number of steps to use in the approximation
 * @param func The function to integrate
*/
double integrate_trapezoidal_rule(double start, double end, double num_steps, double (*func)(double));

#endif