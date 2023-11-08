#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
double matmul(int, double*, double*, double*);

int main(int argc, char **argv) {

  /*
    This is the serial main program for CPSC424/524 Assignment #3.

    Author: Andrew Sherman, Yale University

    Date: 9/30/2023

  */

  int N, i, run;
  double *A, *B, *C, *Ctrue;
  long sizeNxN;

  // This array contains the sizes of the test cases
  int sizes[4]={1000,2000,4000,8000};

  // This array contains the file names for the true answers

  char files [4][75]={"/gpfs/gibbs/project/cpsc424/shared/assignments/assignment3/data/C-1000.dat",\
                      "/gpfs/gibbs/project/cpsc424/shared/assignments/assignment3/data/C-2000.dat",\
                      "/gpfs/gibbs/project/cpsc424/shared/assignments/assignment3/data/C-4000.dat",\
                      "/gpfs/gibbs/project/cpsc424/shared/assignments/assignment3/data/C-8000.dat"};
  
  double wctime, Fnorm;

  FILE *fptr;

  // Print a table heading
  printf("Matrix multiplication times:\n   N      TIME (secs)    F-norm of Error\n -----   -------------  -----------------\n");

  // Now run the four test cases
  for (run=0; run<4; run++) {
    N = sizes[run];

    sizeNxN = N*N; 

    A = (double *) calloc(sizeNxN, sizeof(double));
    B = (double *) calloc(sizeNxN, sizeof(double));
    C = (double *) calloc(sizeNxN, sizeof(double));

    srand(12345); // Use a standard seed value for reproducibility

    // This code loads A row by row, and it loads B column by column.
    // If you don't load the matrix entries in the right order, then your results 
    // won't match the correct answers.
    for (i=0; i<sizeNxN; i++) A[i] = ((double) rand()/(double)RAND_MAX);
    for (i=0; i<sizeNxN; i++) B[i] = ((double) rand()/(double)RAND_MAX);

    // Time the serial matrix multiplication computation
    wctime = matmul(N, A, B, C);

    free(A);
    free(B);

    

    // Remainder of the code checks the result against the correct answer (read into Ctrue)
       Ctrue = (double *) calloc(sizeNxN, sizeof(double));

    fptr = fopen(files[run],"r");
    fread(Ctrue, sizeof(double), sizeNxN, fptr);
    fclose(fptr);    
/*
    fptr = fopen("/home/cpsc424_ahs3/cpsc424/private/2023/a3-2023/trueC.dat","w");
    fwrite(C, sizeof(double), sizeNxN, fptr);
    fclose(fptr);    
*/
    //Compute the Frobenius norm of Ctrue-C
    Fnorm = 0.;
    for (i=0; i<sizeNxN; i++) Fnorm += (Ctrue[i]-C[i])*(Ctrue[i]-C[i]);
    Fnorm = sqrt(Fnorm);

    // Print a table row
    printf ("  %5d    %9.4f  %17.12f\n", N, wctime, Fnorm);

    free(Ctrue);  

    free(C);
  }

}
