#ifndef BASICS_H
#define BASICS_H

/*
  Maybe create another instance of the Matrix Struct since the size
  is repeated too much in the code and loses in performance.
*/
typedef struct Matrix {
  /*
    A Matrix struct can only have 2 or 3 as dims value,
    since, for a 1D vector it has indeed rows or colums.
    For precision purposes it will be specified each time
    as a (1,N) or (N,1) vector.
  */
  int dims;
  int *shape;
  double *data;
} Matrix;

Matrix *create_matrix(int, int *);
void free_matrix(Matrix *);
void get_matrix_shape(Matrix *);
void print_matrix(const Matrix *);
void fill_matrix(Matrix *, const double);
Matrix *copy_matrix(Matrix *);

#endif
