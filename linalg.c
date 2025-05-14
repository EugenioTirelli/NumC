#include <omp.h>
#include "linalg.h"
#include "basics.h"
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix_generator.h"
#include "operators.h"


bool is_squared_matrix(Matrix* m) {
    assert(m->dims == 2);
    return (m->shape[0] == m->shape[1]);
}



void lu_decomposition(Matrix* A, Matrix* L, Matrix* U) {
  assert(A != NULL && L != NULL && U != NULL);
  assert(is_squared_matrix(A));
  assert(is_squared_matrix(L));
  assert(is_squared_matrix(U));

  int n = A->shape[0];
  double* Adata = A->data;
  double* Ldata = L->data;
  double* Udata = U->data;

  #pragma omp parallel for
  for (int j = 0; j < n; j++) {

    int pivot = j;
    double maxVal = fabs(Adata[j * n + j]);
    for (int i = j + 1; i < n; i++) {
      double currentVal = fabs(Adata[i * n + j]);
      if (currentVal > maxVal) {
        pivot = i;
        maxVal = currentVal;
      }
    }
    if (pivot != j) {
      swap2d_rows(A, j, pivot);
      for (int k = 0; k < j; k++) {
        double temp = Ldata[j * n + k];
        Ldata[j * n + k] = Ldata[pivot * n + k];
        Ldata[pivot * n + k] = temp;
      }
    }

    #pragma omp parallel for
    for (int i = j; i < n; i++) {
      double sum = 0.0;
      for (int k = 0; k < j; k++) {
        sum += Ldata[i * n + k] * Udata[k * n + j];
      }
      Ldata[i * n + j] = Adata[i * n + j] - sum;
    }

    if (fabs(Ldata[j * n + j]) < 1e-12) {
      fprintf(
        stderr,
        "Error in LU decomposotion. Matrix may have 0 determinant.\n");
      exit(EXIT_FAILURE);
    }

    #pragma omp parallel for
    for (int i = j + 1; i < n; i++) {
      double sum = 0.0;
      for (int k = 0; k < j; k++) {
        sum += Ldata[j * n + k] * Udata[k * n + i];
      }
      Udata[j * n + i] = (Adata[j * n + i] - sum) / Ldata[j * n + j];
    }
  }
}



Matrix* forward_LU_subs(Matrix *L, Matrix *b) {
  int n = L->shape[0];
  Matrix* y = create_matrix(b->dims, b->shape);
  #pragma omp parallel for
  for (int i = 0; i < n; i++) {
    double sum = 0.0;
    for (int k = 0; k < i; k++) {
      sum += L->data[i * n + k] * y->data[k];
    }
    y->data[i] = (b->data[i] - sum) / L->data[i * n + i];
  }
  return y;
}



Matrix* backward_LU_subs(Matrix *U, Matrix *y) {
  int n = U->shape[0];
  Matrix* x = create_matrix(y->dims, y->shape);

  #pragma omp parallel for
  for (int i = n - 1; i >= 0; i--) {
    double sum = 0.0;
    for (int k = i + 1; k < n; k++) {
      sum += U->data[i * n + k] * x->data[k];
    }
    x->data[i] = y->data[i] - sum;
  }
  return x;
}



Matrix* solve_lin_sys(Matrix *A, Matrix* b) {
  Matrix *L = create_matrix(A->dims, A->shape);
  Matrix *U = diag_k(A->shape[0], 1.0);
  Matrix *y = create_matrix(b->dims, b->shape);
  Matrix *x = create_matrix(b->dims, b->shape);

  lu_decomposition(A, L, U);
  y = forward_LU_subs(L, b);
  x = backward_LU_subs(U, y);

  free_matrix(L);
  free_matrix(U);
  free_matrix(y);

  return x;
}


double det(Matrix* A) {
  Matrix *L = create_matrix(A->dims, A->shape);
  Matrix *U = diag_k(A->shape[0], 1.0);

  int n = L->shape[0];
  lu_decomposition(A, L, U);
  double d = 1.0;

  for (int i = 0; i < n; i++) {
    d *= L->data[i * (n + 1)];
  }

  return d;
}



double norm(Matrix * m) {
  assert (m != NULL);
  assert (m->shape[0] == 1 || m->shape[1] == 1);
  double res = 0.0;
  int m_size = get_matrix_size(m);
  for (int i = 0; i < m_size; i++) {
    res += pow(m->data[i], 2);
  }
  return sqrt(res);
}



Matrix* inv(Matrix* A) {
  assert(A != NULL);
  assert(is_squared_matrix(A));
  int n = A->shape[0];
  Matrix* L = create_matrix(2, A->shape);
  Matrix* U = create_matrix(2, A->shape);
  Matrix* Inv = create_matrix(2, A->shape);
  lu_decomposition(A, L, U);

  #pragma omp parallel for
  for (int col = 0; col < n; col++) {
    Matrix* b = create_matrix(2, (int[2]){n, 1});
    for (int i = 0; i < n; i++) {
      b->data[i] = (i == col) ? 1.0 : 0.0;
    }
    Matrix* y = forward_LU_subs(L, b);
    Matrix* x = backward_LU_subs(U, y);

    for (int row = 0; row < n; row++) {
      Inv->data[row * n + col] = x->data[row];
    }
    free_matrix(b);
    free_matrix(y);
    free_matrix(x);
  }
  free_matrix(L);
  free_matrix(U);
  return Inv;
}



Matrix* add_mat(Matrix* a, Matrix* b) {
  assert(b != NULL);
  assert(a->dims == b->dims);
  assert(a->shape == b->shape);

  Matrix* c = create_matrix(a->dims, a->shape);
  int rows = a->shape[0];
  int cols = a->shape[1];

  #pragma omp parallel for
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      int index = i * cols + j;
      c->data[index] = a->data[index] + b->data[index];
    }
  }

  return c;
}


Matrix* sub_mat(Matrix* a, Matrix* b) {
  assert(b != NULL);
  assert(a->dims == b->dims);
  assert(a->shape == b->shape);

  Matrix* c = create_matrix(a->dims, a->shape);
  int rows = a->shape[0];
  int cols = a->shape[1];

  #pragma omp parallel for
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      int index = i * cols + j;
      c->data[index] = a->data[index] - b->data[index];
    }
  }

  return c;
}



