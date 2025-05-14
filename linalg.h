#ifndef LINALG_H
#define LINALG_H

#include "basics.h"
#include <stdbool.h>

bool is_squared_matrix(Matrix *);
void lu_decomposition(Matrix *, Matrix *, Matrix *);
Matrix *forward_LU_subs(Matrix *, Matrix *);
Matrix *backward_LU_subs(Matrix *, Matrix *);
Matrix *solve_lin_sys(Matrix *, Matrix *);
double det(Matrix *);
double norm(Matrix *);
Matrix *inv(Matrix *);
Matrix *add_mat(Matrix *, Matrix *);
Matrix *sub_mat(Matrix *, Matrix *);

#endif
