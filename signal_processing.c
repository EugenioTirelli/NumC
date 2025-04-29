#include "signal_processing.h"
#include "basics.h"
#include "cbasics.h"
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

CMatrix *dft1d(const Matrix *m) {
  assert(m != NULL);

  CMatrix *df_matrix = create_cmatrix(m->dims, m->shape);

  int size = 1;
  for (int i = 0; i < m->dims; i++) {
    size *= m->shape[i];
  }

  for (int j = 0; j < size; j++) {
    double complex sum = 0.0 + 0.0 * I;
    for (int k = 0; k < size; k++) {
      double angle = 2 * M_PI * j * k / size;
      sum += m->data[k] * cexp(-I * angle);
    }
    df_matrix->cdata[j] = sum;
  }

  return df_matrix;
}
