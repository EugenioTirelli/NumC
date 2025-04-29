#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "basics.h"
#include "cbasics.h"
#include "operators.h"
#include "matrix_generator.h"
#include "signal_processing.h"

// gcc main.c basics.c operators.c -o main -lm

int main(int argc, char *argv[])
{
  Matrix *t = linspacef(0.0, 10.0, 10, 1, 0);
  print_matrix(t);
  Matrix* yt = sin_op(t);

  CMatrix* yf = dft1d(yt);
  print_cmatrix(yf);
 
  free_matrix(t);
  free_cmatrix(yf);


  printf("End of program.\n");
  return 0;
}
