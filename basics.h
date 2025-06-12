#ifndef BASICS_H
#define BASICS_H

#define max2(x, y) (((x) >= (y)) ? (x) : (y))
#define min2(x, y) (((x) <= (y)) ? (x) : (y))

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

typedef struct SparseMatrix {
  /*
     A SparseMatrix can be stored just by indexed and values
     since all its value except the ones in entry are zeros.
     This optimizes the memory allocation.
  */
  int num_rows;
  int num_cols;
  int num_of_nnz_entries;
  int *row_indexes;
  int *col_indexes;
  double *data;
} SparseMatrix;

Matrix *create_matrix(int, int *);
void free_matrix(Matrix *);

SparseMatrix *create_sparse_matrix(int, int, int);
SparseMatrix *copy_sparse_matrix(SparseMatrix *);
void print_sparse_matrix(SparseMatrix *);
void get_sparse_matrix_shape(SparseMatrix *);
void free_sparse_matrix(SparseMatrix *);

void get_matrix_shape(Matrix *);
int get_matrix_size(Matrix *);
void print_matrix(const Matrix *);
void fill_matrix(Matrix *, const double);
Matrix *copy_matrix(Matrix *);

int get_1darray_len(int *);

Matrix *SparseToMatrix(SparseMatrix *);
SparseMatrix *MatrixToSparse(Matrix *);

#endif
