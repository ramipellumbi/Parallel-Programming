#include "timing.h"
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

double matmul(int N, double *A, double *B, double *C)
{

    /*
      This is part of the serial version of matrix multiplication for CPSC424/524 Assignment #3.

      This matmul() function multiplies a block row of A times a block column of B
      to compute a rectangular submatrix of the resulting C matrix (where C = A times B).

      If the block row of A is MxN, and the block row of B is NxM, then the submatrix of C is MxM.
      A typical parallel code would use M = N/p (for p processes), assuming N is a multiple of p.

      Author: Andrew Sherman, Yale University

      Date: 10/20/2023

    */

    int i, j, k;
    int iA; // iA is a pointer into the A matrix (stored by rows)
    int jB; // jB is a pointer into the B matrix (stored by columns)
    int iC; // iC is a pointer into the C matrix (stored by rows)

    double wctime0, wctime1, cputime; // only used for timing

    timing(&wctime0, &cputime);

    // This loop computes the matrix-matrix product assuming that the matrices
    // A and C are stored row-by-row, while the matrix B is stored column-by-column
    iC = 0;
    for (i = 0; i < N; i++)
    {
        iA = i * N; // Initializes row pointer for row i of A to skip over all previous rows of A.
                    // Note: if we're working on row i, then the total number of double entries for
                    // rows 0 through i-1 is i*N. Thus, iA points to the first entry of row i.

        for (j = 0; j < N; j++, iC++)
        {
            // We're going to compute the dot product of row i of A and column j of B.
            // iC points to the entries of row i of C. (We compute C[i,j] as the dot
            // product of row i of A with column j of B.)

            jB = j * N; // Initializes column pointer for col j of B to skip over all previous cols of B.
                        // Note: if we're working on col j, then the total number of double entries for
                        // cols 0 through j-1 is j*N. Thus, jB points to the first entry of col j.

            C[iC] = 0.;
            for (k = 0; k < N; k++)
                C[iC] += A[iA + k] * B[jB + k]; // Each iteration on k adds in the product
                                                // of A[i,k] times B[k,j].
        }
    }

    timing(&wctime1, &cputime);
    return (wctime1 - wctime0);
}
