#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include "basics.h"
#include "cbasics.h"
#include "operators.h"
#include "matrix_generator.h"
#include "signal_processing.h"
#include "linalg.h"

// gcc main.c basics.c operators.c -o main -lm

int main(int argc, char *argv[])
{

  double Adata[] = {
    0,2,3,4,
    6,9,3,7,
    0, 2, 5, 0,
    3, 1, 2, 3};

  double data2[] = {
    6, -2, 2,
    -2, 3, -1,
    2, -1, 3};

  double bdata[] = {
    6, 5, 3
  };

  Matrix *A = create_matrix(2, (int[]){3,3});
  Matrix *b = create_matrix(2, (int[]){3, 1});
  memcpy(A->data, data2, sizeof(data2));
  memcpy(b->data, bdata, sizeof(bdata));
//  Matrix* test_mat = rndmat(2, (int[]){3000,3000}, -10.0, 10.0);

  double determinantA = det(A);
  printf("%lf\n", determinantA);

  printf("End of program.\n");
  return 0;
}
