#include <complex.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "basics.h"
#include <math.h>


Matrix* sin_op(const Matrix* mat_) {
  if (!mat_) {
    fprintf(stderr, "Input matrix is NULL in sin_operator\n");
    return NULL;
  }

  Matrix* result = create_matrix(mat_->dims, mat_->shape);
  int size = 1;
  for (int i = 0; i < mat_->dims; i++) {
    size *= mat_->shape[i];
  }

  for (int i = 0; i < size; i++) {
    result->data[i] = sin(mat_->data[i]);
  }
  return result;
}

Matrix* cos_op(const Matrix* mat_) {
  if (!mat_) {
    fprintf(stderr, "Input matrix is NULL in cos_op\n");
    return NULL;
  }

  Matrix* result = create_matrix(mat_->dims, mat_->shape);
  int size = 1;
  for (int i = 0; i < mat_->dims; i++) {
    size *= mat_->shape[i];
  }

  for (int i = 0; i < size; i++) {
    result->data[i] = cos(mat_->data[i]);
  }
  return result;
}

Matrix* abs_op(const Matrix* mat_) {
  if (!mat_) {
    fprintf(stderr, "Input matrix is NULL in abs_op\n");
    return NULL;
  }

  Matrix* result = create_matrix(mat_->dims, mat_->shape);
  int size = 1;
  for (int i = 0; i < mat_->dims; i++) {
    size *= mat_->shape[i];
  }

  for (int i = 0; i < size; i++) {
    result->data[i] = fabs(mat_->data[i]); }
  return result;
}

Matrix* exp_op(const Matrix* mat_) {
  if (!mat_) {
    fprintf(stderr, "Input matrix is NULL in exp_op\n");
    return NULL;
  }

  Matrix* res_ = create_matrix(mat_->dims, mat_->shape);
  int size = 1;

  for (int i = 0; i < mat_->dims; i++) {
    size *= mat_->shape[i];
  }

  for (int i = 0; i < size; i++) {
    res_->data[i] = exp(mat_->data[i]);
  }

  return res_;
}

Matrix* const_mul(const Matrix* mat_, const double k) {
  if (!mat_) {
    fprintf(stderr, "Input matrix is NULL in const_mul\n");
    return NULL;
  }

  Matrix* ret_ = create_matrix(mat_->dims, mat_->shape);
  int size = 1;

  for (int i = 0; i < mat_->dims; i++) {
    size *= mat_->shape[i];
  }

  for (int i = 0; i < size; i++) {
    ret_->data[i] = mat_->data[i] * k;
  }

  return ret_;
}


Matrix* const_add(const Matrix* mat_, const double k) {
  if (!mat_) {
    fprintf(stderr, "Input matrix is NULL in const_add\n");
    return NULL;
  }

  Matrix* ret_ = create_matrix(mat_->dims, mat_->shape);
  int size = 1;

  for (int i = 0; i < mat_->dims; i++) {
    size *= mat_->shape[i];
  }

  for (int i = 0; i < size; i++) {
    ret_->data[i] = mat_->data[i] + k;
  }

  return ret_;
}

Matrix* const_sub(const Matrix* mat_, const double k) {
  if (!mat_) {
    fprintf(stderr, "Input matrix is NULL in const_sub\n");
    return NULL;
  }

  Matrix* ret_ = create_matrix(mat_->dims, mat_->shape);
  int size = 1;

  for (int i = 0; i < mat_->dims; i++) {
    size *= mat_->shape[i];
  }

  for (int i = 0; i < size; i++) {
    ret_->data[i] = mat_->data[i] - k;
  }

  return ret_;
}


Matrix* const_div(const Matrix* mat_, const double k) {
  assert (k > 0);

  if (!mat_) {
    fprintf(stderr, "Input matrix is NULL in const_add\n");
    return NULL;
  }

  Matrix* ret_ = create_matrix(mat_->dims, mat_->shape);
  int size = 1;

  for (int i = 0; i < mat_->dims; i++) {
    size *= mat_->shape[i];
  }

  for (int i = 0; i < size; i++) {
    ret_->data[i] = mat_->data[i] / k;
  }

  return ret_;
}

Matrix* rad2deg(Matrix* mat) {
  assert (mat != NULL);

  Matrix* deg_mat = create_matrix(mat->dims, mat->shape);

  int size = 1;
  for (int i = 0; i < mat->dims; i++) {
    size *= mat->shape[i];
  }

  for (int i = 0; i < size; i++) {
    deg_mat->data[i] = (180. / M_PI) * mat->data[i];
  }

  return deg_mat;
}


Matrix* deg2rad(Matrix* mat) {
  assert (mat != NULL);

  Matrix* rad_mat = create_matrix(mat->dims, mat->shape);

  int size = 1;
  for (int i = 0; i < mat->dims; i++) {
    size *= mat->shape[i];
  }

  for (int i = 0; i < size; i++) {
    rad_mat->data[i] = (M_PI / 180.) * mat->data[i];
  }

  return rad_mat;
}



Matrix* transpose2d(Matrix* mat_) {
  int tshape[] = {mat_->shape[1], mat_->shape[0]};

  Matrix* ret_ = create_matrix(mat_->dims, tshape);

  int size = 1;

  for (int i = 0; i < mat_->dims; i++) {
    size *= mat_->shape[i];
  }

  int rows = mat_->shape[1];
  int cols = mat_->shape[0];

  for (int i = 0; i <= rows; i++) {
    for (int j = 0; j <= cols; j++) {
      ret_->data[j * cols + i] = mat_->data[i * rows + j];
    }
  }

  return ret_;
}

Matrix* dot_prod(Matrix* a, Matrix* b) {
  assert(a != NULL);
  assert(b != NULL);

  assert(a->dims == 2);
  assert(b->dims == 2);

  if (a->shape[1] != b->shape[0]) {
    fprintf(stderr, "Cannot perform dot_prod; shapes mismatch.\n");
    return NULL;
  }

  int pshape[2] = {a->shape[0], b->shape[1]};
  Matrix* dprod = create_matrix(2, pshape);

  if (!dprod || !dprod->data) {
    fprintf(stderr, "Failed to allocate memory for dot product result.\n");
    return NULL;
  }

  int rowsA = a->shape[0];
  int colsA = a->shape[1];
  int colsB = b->shape[1];

  for (int i = 0; i < rowsA; i++) {
    for (int j = 0; j < colsB; j++) {
      double sum = 0.0;
      for (int k = 0; k < colsA; k++) {
        sum += a->data[i * colsA + k] * b->data[k * colsB + j];
      }
      dprod->data[i * colsB + j] = sum;
    }
  }

  return dprod;
}



double maxel(Matrix* m) {
    assert(m != NULL);
    int size = 1;
    for (int i = 0; i < m->dims; i++) {
        size *= m->shape[i];
    }

    double max_ = m->data[0];

    for (int i = 1; i < size; i++) {
        if (m->data[i] > max_) {
            max_ = m->data[i];
        }
    }

    return max_;
}


double minel(Matrix* m) {
    assert(m != NULL);
    int size = 1;
    for (int i = 0; i < m->dims; i++) {
        size *= m->shape[i];
    }

    double min_ = m->data[0];

    for (int i = 1; i < size; i++) {
        if (m->data[i] < min_) {
            min_ = m->data[i];
        }
    }

    return min_;
}

Matrix* reverse1d(Matrix* m) {
  assert (m != NULL);
  assert (m -> dims == 2);
  assert (m->shape[0] == 1 || m->shape[1] == 1);
  Matrix* rev_m = create_matrix(m->dims, m->shape);

  int size = 1;

  for (int i = 0; i < m->dims; i++) {
    size *= m->shape[i];
  }
  for (int i = 0; i < size; i++) {
    rev_m->data[i] = m->data[size - i - 1];
  }
  return rev_m;
}
