#ifndef CBASICS_H
#define CBASICS_H

#include "basics.h"
#include <assert.h>
#include <complex.h>
#include <memory.h>
#include <stdio.h>

typedef struct CMatrix {
  int cdims;
  int *cshape;
  _Complex double *cdata;
} CMatrix;

CMatrix *create_cmatrix(int dims, int *shape);
void free_cmatrix(CMatrix *mat);
void get_cmatrix_shape(CMatrix *mat);
void print_cmatrix(const CMatrix *mat);
CMatrix *copy_cmatrix(CMatrix *m);

#endif
