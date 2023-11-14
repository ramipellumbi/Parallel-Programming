#ifndef MATMUL_H
#define MATMUL_H

double matmul(int N, double *A, double *B, double *C);

void gemm(int M, int N, double *A, double *B, double *C, int offset);

void gemm_k(int M, int N, int K, double *A, double *B, double *C);

#endif