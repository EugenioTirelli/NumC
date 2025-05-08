#ifndef OPERATORS_H
#define OPERATORS_H

#include "basics.h"
#include <math.h>
#include <stdlib.h>

Matrix *sin_op(Matrix *);
Matrix *cos_op(Matrix *);
Matrix *abs_op(Matrix *);
Matrix *exp_op(Matrix *);
Matrix *sqrt_op(Matrix *);
Matrix *loge_op(Matrix *);
Matrix *log10_op(Matrix *);
Matrix *log2_op(Matrix *);
Matrix *power_op(Matrix *);
double sum1d(Matrix *);
double prod1d(Matrix *);
Matrix *sum2d(Matrix *, int);
Matrix *prod2d(Matrix *, int);
Matrix *const_mul(Matrix *, const double);
Matrix *const_add(Matrix *, const double);
Matrix *const_sub(Matrix *, const double);
Matrix *const_div(Matrix *, const double);
Matrix *rad2deg(Matrix *);
Matrix *deg2rad(Matrix *);
Matrix *transpose2d(Matrix *);
Matrix *dot_prod(Matrix *, Matrix *);
double maxel(Matrix *);
double minel(Matrix *);
Matrix *reverse1d(Matrix *);
#endif
