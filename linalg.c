#include "linalg.h"
#include "basics.h"
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <stdio.h>
#include "matrix_generator.h"
#include "operators.h"


bool is_squared_matrix(Matrix* m) {
    assert(m->dims == 2);
    return (m->shape[0] == m->shape[1]);
}

// void lu_decomposition(Matrix* A, Matrix* L, Matrix* U)
// {
//     assert(A != NULL && L != NULL && U != NULL); assert(is_squared_matrix(A));
//     assert(is_squared_matrix(L));
//     assert(is_squared_matrix(U));
// 
//     int n = A->shape[0];
//     // Access the raw pointers once to avoid repeated structure lookups
//     double* restrict Adata = A->data;
//     double* restrict Ldata = L->data;
//     double* restrict Udata = U->data;
// 
//     // -----------------------------------------------------------------------
//     // 1) Quick initialization for the first column of L and the first row of U
//     //    (Not strictly necessary in a standard Crout approach, but kept for
//     //    consistency with the original code's logic.)
//     // -----------------------------------------------------------------------
//     for (int i = 0; i < n; i++) {
//         // L(i, 0) = A(i, 0)
//         Ldata[i * n + 0] = Adata[i * n + 0];
// 
//         // For i > 0, U(0, i) = A(0, i) / L(0, 0), but original code uses A->data[i]
//         // so we replicate the logic: U(i, 0) if i >= 1.
//         // (This step typically belongs in the main loop, but is retained from the original.)
//         if (i >= 1) {
//             // U(0, i) in row-major is Udata[0*n + i], which is Udata[i].
//             // A(0, i) is Adata[0*n + i] => Adata[i].
//             Udata[i] = Adata[i] / Adata[0];
//         }
//     }
// 
//     // -----------------------------------------------------------------------
//     // 2) Main loop for Crout's LU Decomposition
//     //    For each column j from 0 to n-1:
//     //      a) L(i, j) = A(i, j) - SUM_{k=0..j-1} [L(i, k)*U(k, j)] for i>=j
//     //      b) U(j, i) = (A(j, i) - SUM_{k=0..j-1} [L(j, k)*U(k, i)]) / L(j, j) for i>j
//     // -----------------------------------------------------------------------
//     for (int j = 0; j < n; j++) {
// 
//         // 2a) Compute L(i, j) for i = j..n-1
//         for (int i = j; i < n; i++) {
//             double sum = 0.0;
//             // Pointers to row i in L and row i in A, to reduce repeated i*n calculations
//             double* Li = &Ldata[i * n];
//             double* Ai = &Adata[i * n];
// 
//             for (int k = 0; k < j; k++) {
//                 sum += Li[k] * Udata[k * n + j];  // L(i, k)*U(k, j)
//             }
//             // L(i, j) = A(i, j) - sum
//             Li[j] = Ai[j] - sum;
//         }
// 
//         // Check for zero pivot L(j, j)
//         if (Ldata[j * n + j] == 0.0) {
//             fprintf(stderr, "Zero pivot encountered at L[%d, %d]. LU decomposition fails.\n", j, j);
//             exit(EXIT_FAILURE);
//         }
// 
//         // 2b) Compute U(j, i) for i = j+1..n-1
//         {
//             // Pre-cache pointer to row j in L
//             double* Lj = &Ldata[j * n];
//             // Also pointer to row j in U, though we store each result in U(j, i).
//             for (int i = j + 1; i < n; i++) {
//                 double sum = 0.0;
//                 for (int k = 0; k < j; k++) {
//                     sum += Lj[k] * Udata[k * n + i]; // L(j, k)*U(k, i)
//                 }
//                 // U(j, i) = (A(j, i) - sum) / L(j, j)
//                 Udata[j * n + i] = (Adata[j * n + i] - sum) / Lj[j];
//             }
//         }
//     }
// }
void lu_decomposition(Matrix* A, Matrix* L, Matrix* U) {
  assert(A != NULL && L != NULL && U != NULL);
  assert(is_squared_matrix(A));
  assert(is_squared_matrix(L));
  assert(is_squared_matrix(U));

  int n = A->shape[0];
  double* Adata = A->data;
  double* Ldata = L->data;
  double* Udata = U->data;

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

