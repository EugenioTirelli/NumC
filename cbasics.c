#include <assert.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "basics.h"
#include "cbasics.h"
#include <complex.h>


CMatrix *create_cmatrix(int dims, int* shape) {
  assert(dims > 1);

  CMatrix* mat = (CMatrix*) malloc(sizeof(Matrix));

  if (mat == NULL) {
    perror("Memory allocation failed.");
    exit(EXIT_FAILURE);
  }

  mat->cdims = dims;
  mat->cshape = (int*) malloc(sizeof(int) * dims);
  if (mat->cshape == NULL) {
    perror("Memory allocation failed.");
    free(mat->cshape);
    free(mat);
    exit(EXIT_FAILURE);
  }

  int size = 1;

  for (int i = 0; i < dims; i++) {
    mat->cshape[i] = shape[i];
    size *= shape[i];
  }

  mat->cdata = (complex double*) calloc(size, sizeof(double complex));
  if (mat->cdata == NULL) {
    perror("Memory allocation failed.");
    free(mat->cshape);
    free(mat);
    exit(EXIT_FAILURE);
  }

  return mat;
}


void free_cmatrix(CMatrix* mat) {
  free(mat->cdata);
  free(mat->cshape);
  free(mat);
}


void get_cmatrix_shape(CMatrix *mat) {
  for (int i = 0; i < mat->cdims; i++) {
    printf("%8d", mat->cshape[i]);
  }
  printf("\n");
}


void print_cmatrix(const CMatrix* mat) {
  if (mat->cdims == 2) {
    for (int i = 0; i < mat->cshape[0]; i++) {
      for (int j = 0; j < mat->cshape[1]; j++) {
        printf("%8lf%8lfj ", creal(mat->cdata[i * mat->cshape[1] + j]),
                                    cimag(mat->cdata[i * mat->cshape[1] + j])) ;
      }
      printf("\n");
    }
  } else if (mat->cdims == 3) {
    for (int i = 0; i < mat->cshape[0]; i++) {
      for (int j = 0; j < mat->cshape[1]; j++) {
        for (int k = 0; k < mat->cshape[2]; k++) {
          printf("%8lf%8lfj ", creal(mat->cdata[i * mat->cshape[1] * mat->cshape[2] + j * mat->cshape[2] + k]), 
                                     cimag(mat->cdata[i * mat->cshape[1] * mat->cshape[2] + j * mat->cshape[2] + k]));
        }
        printf("\n");
      }
      printf("\n");
    }
  } else {
    printf("Dimensione della matrice non supportata.\n");
  }
}


CMatrix* copy_cmatrix(CMatrix* m) {
  CMatrix* ccm = create_cmatrix(m->cdims, m->cshape);
  memcpy(ccm, m, sizeof(CMatrix));
  return ccm;
}
