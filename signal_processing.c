#include "signal_processing.h"
#include "basics.h"
#include "cbasics.h"
#include "operators.h"
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


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


// ---------- 1D Convolution Function ----------
//
// Mode Definitions:
//   mode 1: FULL convolution, length = a_size + b_size - 1
//   mode 2: VALID convolution, length = a_size - b_size + 1  (assuming a_size >= b_size)
//   mode 3: SAME  convolution, length = max(a_size, b_size)
//   mode 4: SUBSET extraction: returns only a subset (here, first 10 elements)
//
// Note: A precise implementation of VALID and SAME would require adjusting the indices
// used in the summation. Here we compute the full convolution first and then extract
// the relevant region for mode 2, 3 or 4.
Matrix* convolve1d(Matrix* a, Matrix* b, int mode) {
  assert(a != NULL && b != NULL);
  assert(a->shape[0] != 2 || a->shape[1] != 2);
  assert(b->shape[0] != 2 || b->shape[1] != 2);

  int a_size = get_matrix_size(a);
  int b_size = get_matrix_size(b);

  int full_length = a_size + b_size - 1;
  int desired_length = full_length;
  int offset = 0;

  switch (mode) {
    case 1:
      desired_length = full_length;
      offset = 0;
      break;
    case 2:
      desired_length = (a_size - b_size) + 1;
      offset = b_size - 1;
      break;
    case 3:
      desired_length = max2(a_size, b_size);
      offset = (full_length - desired_length) / 2;
      break;
    default:
      printf(stderr, "Please provide a mode.\n");
      return NULL;
      break;
  }

  double *temp = (double *)calloc(full_length, sizeof(double));
  if (!temp) return NULL;
  for (int i = 0; i < a_size; i++) {
    for (int j = 0; j < b_size; j++) {
      temp[i + j] += a->data[i] * b->data[j];
    }
  }

  int new_shape[2] = {1, desired_length};
  Matrix *out = create_matrix(2, new_shape);
  if (!out) {
    free(temp);
    return NULL;
  }

  for (int i = 0; i < desired_length; i++) {
    if ((i + offset) < full_length)
      out->data[i] = temp[i + offset];
    else
      out->data[i] = 0.0;
  }
  free(temp);
  return out;
}
