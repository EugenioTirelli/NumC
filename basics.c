#include <assert.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "basics.h"



Matrix *create_matrix(int dims, int* shape) {
  assert(dims > 1);

  Matrix* mat = (Matrix*) malloc(sizeof(Matrix));

  if (mat == NULL) {
    perror("Memory allocation failed.");
    exit(EXIT_FAILURE);
  }

  mat->dims = dims;
  mat->shape = (int*) malloc(sizeof(int) * dims);
  if (mat->shape == NULL) {
    perror("Memory allocation failed.");
    free(mat->shape);
    free(mat);
    exit(EXIT_FAILURE);
  }


  int size = 1;
  for (int i = 0; i < dims; i++) {
    mat->shape[i] = shape[i];
    size *= shape[i];
  }

  mat->data = (double*) calloc(size, sizeof(double));
  if (mat->data == NULL) {
    perror("Memory allocation failed.");
    free(mat->shape);
    free(mat);
    exit(EXIT_FAILURE);
  }

  return mat;
}


void free_matrix(Matrix* mat) {
  free(mat->data);
  free(mat->shape);
  free(mat);
}


void get_matrix_shape(Matrix *mat) {
  for (int i = 0; i < mat->dims; i++) {
    printf("%8d", mat->shape[i]);
  }
  printf("\n");
}

int get_matrix_size(Matrix *mat) {
  int size = 1;
  for (int i = 0; i < mat->dims; i++) {
    size *= mat->shape[i];
  }
  return size;
}

void print_matrix(const Matrix* mat) {
  if (mat->dims == 2) {
    for (int i = 0; i < mat->shape[0]; i++) {
      for (int j = 0; j < mat->shape[1]; j++) {
        printf("%lf ", mat->data[i * mat->shape[1] + j]);
      }
      printf("\n");
    }
  } else if (mat->dims == 3) {
    for (int i = 0; i < mat->shape[0]; i++) {
      for (int j = 0; j < mat->shape[1]; j++) {
        for (int k = 0; k < mat->shape[2]; k++) {
          printf("%lf ", mat->data[i * mat->shape[1] * mat->shape[2] + j * mat->shape[2] + k]);
        }
        printf("\n");
      }
      printf("\n");
    }
  } else {
    printf("Dimensione della matrice non supportata.\n");
  }
}


void fill_matrix(Matrix* m, const double n) {
  assert (m != NULL);
  int size = 1;
  for (int i = 0; i < m->dims; i++) {
    size *= m->shape[i];
  }
  for (int i = 0; i < size; i++) {
    m->data[i] = n;
  }
}


Matrix* copy_matrix(Matrix* m) {
  Matrix* cm = create_matrix(m->dims, m->shape);
  memcpy(cm, m, sizeof(Matrix));
  return cm;
}
