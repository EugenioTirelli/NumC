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
  Matrix* test_mat = rndmat(2, (int[]){10, 4}, -10.0, 10.0);

  print_matrix(test_mat);
  printf("\n");

  swap2d_rows(test_mat, 0, 9);
  print_matrix(test_mat);
  printf("\n");

  swap2d_cols(test_mat, 2, 3);
  print_matrix(test_mat);
  printf("\n");

  printf("End of program.\n");
  return 0;
}
