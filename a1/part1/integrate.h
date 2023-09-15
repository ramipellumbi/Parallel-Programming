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

/**
 * Have the interval:
 *      |--------X--------|--------X--------|  ...   |--------X--------|
 *      0      1/(2N)    1/N      3/(2N)   2/N   (N-1)/N   (2N-1)/2N  N/N
 *
 * Where the:
 *   - | are the endpoints that are distance 1/N apart 
 *   - X are the midpoints that are distance 1/N apart.
 * 
 * 0 can represent any arbitrary start and 1 can represent any arbitrary end.
 * 
 * @param start The start of the interval to integrate over
 * @param end The end of the interval to integrate over
 * @param num_steps The number of steps to use in the approximation
 * @param func The function to integrate
 */
double integrate_midpoint_rule(double start, double end, double num_steps, double (*func)(double));

#endif