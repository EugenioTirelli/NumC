#ifndef OPERATORS_H
#define OPERATORS_H

#include "basics.h"
#include <math.h>
#include <stdlib.h>

Matrix *sin_op(const Matrix *);
Matrix *cos_op(const Matrix *);
Matrix *abs_op(const Matrix *);
Matrix *exp_op(const Matrix *);
Matrix *const_mul(const Matrix *, const double);
Matrix *const_add(const Matrix *, const double);
Matrix *const_sub(const Matrix *, const double);
Matrix *const_div(const Matrix *, const double);
Matrix *rad2deg(Matrix *);
Matrix *deg2rad(Matrix *);
Matrix *transpose2d(const Matrix *);
Matrix *dot_prod(Matrix *, Matrix *);
double maxel(Matrix *);
double minel(Matrix *);
Matrix *reverse1d(Matrix *);
#endif
