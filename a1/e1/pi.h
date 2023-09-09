#ifndef PI_H
#define PI_H

// Function to be integrated for the approximation
double function_to_integrate(double x);

// Implementing the Trapezoidal Rule for numerical integration
double integrate_trapezoidal_rule(double start, double end, double num_steps, double (*func)(double));

#endif 
