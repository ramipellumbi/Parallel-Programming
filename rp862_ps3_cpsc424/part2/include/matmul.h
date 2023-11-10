#ifndef MATMUL_H
#define MATMUL_H

/**
 * @param N number of columns of A / number of rows of B
 * @param A An M x N matrix (stored by rows)
 * @param B An N x M matrix (stored by columns)
 * @param C Result of the multiple (stored by rows)
 *
 * Does the matrix multiply of A and B and stores it in C
 */
double matmul(int N, double *A, double *B, double *C);

#endif