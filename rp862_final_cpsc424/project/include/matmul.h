#ifndef MATMUL_H
#define MATMUL_H

/**
 * Returns the wall clock time elapsed in the naive matrix multiplication between A and B
 *
 * @param A double[M*N] (row wise storage)
 * @param B double[N*K] (column wise storage)
 * @param C double[L] (row wise storage) where L >= M*K.
 * @param M Number of rows of A
 * @param N Number of columns of A (rows of B)
 * @param K Number or columns of B
 * @param offset Offset to store the result in C. If L = M*K, offset should be 0.
 */
double matrix_multiply_naive(double *A, double *B, double *C, int M, int N, int K, int offset);

/**
 * Returns the wall clock time elapsed in the matrix multiplication between A and B
 * where AVX instructions are used on the dot product
 *
 * @param A double[M*N] (row wise storage)
 * @param B double[N*K] (column wise storage)
 * @param C double[L] (row wise storage) where L >= M*K.
 * @param M Number of rows of A
 * @param N Number of columns of A (rows of B)
 * @param K Number or columns of B
 * @param offset Offset to store the result in C. If L = M*K, offset should be 0.
 */
double matrix_multiply_avx(double *A, double *B, double *C, int M, int N, int K, int offset);

#endif