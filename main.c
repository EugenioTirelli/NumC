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
  Matrix* a = create_matrix(2, (int[]){3,3});
  Matrix* test_matrix = ones(2, (int[]){10000, 10000});

  Matrix* res = sum2d(test_matrix, 1);

  double a_values[] = {1,2,3,4,5,6,7,8,9};
 
  memcpy(a->data, a_values, sizeof(a_values));

  Matrix* s2d_axis0 = sum2d(a, 0);
  Matrix* s2d_axis1 = sum2d(a, 1);

  print_matrix(s2d_axis0);
  print_matrix(s2d_axis1);

  free_matrix(test_matrix);
  free_matrix(s2d_axis0);
  free_matrix(s2d_axis1);


  printf("End of program.\n");
  return 0;
}
