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
  free(mat->data);
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
  int size = get_matrix_size(m);
  for (int i = 0; i < size; i++) {
    m->data[i] = n;
  }
}


Matrix* copy_matrix(Matrix* m) {
  Matrix* cm = create_matrix(m->dims, m->shape);
  memcpy(cm, m, sizeof(Matrix));
  return cm;
}


SparseMatrix *create_sparse_matrix(int nrows, int ncols, int nzero_entries) {
  SparseMatrix* smat = (SparseMatrix*) malloc(sizeof(SparseMatrix));
  if (smat == NULL) {
    perror("Memory allocation failed.");
    exit(EXIT_FAILURE);
  }
  smat->row_indexes = (int*)calloc(nzero_entries, sizeof(int));
  smat->col_indexes = (int*)calloc(nzero_entries, sizeof(int));
  if (smat->row_indexes == NULL || smat->col_indexes == NULL) {
    perror("Memory allocation failed.");
    free(smat->row_indexes);
    free(smat->col_indexes);
    free(smat);
    exit(EXIT_FAILURE);
  }
  smat->data = (double*) calloc(nzero_entries, sizeof(double));
  if (smat->data == NULL) {
    perror("Memory allocation failed.");
    free(smat->data);
    free(smat);
    exit(EXIT_FAILURE);
  }
  return smat;
}


SparseMatrix* copy_sparse_matrix(SparseMatrix* m) {
  SparseMatrix* sm = create_sparse_matrix(m->num_rows,
                                          m->num_cols,
                                          m->num_of_nnz_entries);
  memcpy(sm, m, sizeof(SparseMatrix));
  return sm;
}


void free_sparse_matrix(SparseMatrix* sm) {
  free(sm->row_indexes);
  free(sm->col_indexes);
  free(sm->data);
  free(sm);
}

void get_sparse_matrix_shape(SparseMatrix* sm) {
  printf("%8d %8d %8d\n", sm->num_rows, sm->num_cols, sm->num_of_nnz_entries);
}


void print_sparse_matrix(SparseMatrix* sm) {



}
