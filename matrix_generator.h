#ifndef MATRIX_GENERATOR_H
#define MATRIX_GENERATOR_H

#include "basics.h"

Matrix *zeros(int, int *);
Matrix *ones(int, int *);
Matrix *diag_k(int, double);
Matrix *rndmat(int, int *, const double, const double);
Matrix *linspacef(double, double, int, int, int);
#endif
