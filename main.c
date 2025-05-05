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
  Matrix* a = create_matrix(2, (int[]){1, 6});
  Matrix* b = create_matrix(2, (int[]){1, 3});

  double a_values[] = {4, 1, 0, 7, 6, 5};
  double b_values[] = {1, 2, 3};
 
  memcpy(a->data, a_values, sizeof(a_values));
  memcpy(b->data, b_values, sizeof(b_values));


  Matrix* full = convolve1d(a, b, 1);
  Matrix* valid = convolve1d(a, b, 2);
  Matrix* same = convolve1d(a, b, 3);

  print_matrix(a);
  print_matrix(b);
  printf("\n");
  print_matrix(full);
  print_matrix(valid);
  print_matrix(same);

  free_matrix(a);
  free_matrix(b);
  free_matrix(full);
  free_matrix(valid);
  free_matrix(same);


  printf("End of program.\n");
  return 0;
}
