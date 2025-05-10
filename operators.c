#include <complex.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "basics.h"
#include <math.h>
#include "linalg.h"


Matrix* sin_op(Matrix* mat_) {
  if (!mat_) {
    fprintf(stderr, "Input matrix is NULL in sin_op\n");
    return NULL;
  }

  Matrix* result = create_matrix(mat_->dims, mat_->shape);
  int size = get_matrix_size(mat_);

  for (int i = 0; i < size; i++) {
    result->data[i] = sin(mat_->data[i]);
  }
  return result;
}

Matrix* cos_op(Matrix* mat_) {
  if (!mat_) {
    fprintf(stderr, "Input matrix is NULL in cos_op\n");
    return NULL;
  }

  Matrix* result = create_matrix(mat_->dims, mat_->shape);
  int size = get_matrix_size(mat_);

  for (int i = 0; i < size; i++) {
    result->data[i] = cos(mat_->data[i]);
  }
  return result;
}

Matrix* abs_op(Matrix* mat_) {
  if (!mat_) {
    fprintf(stderr, "Input matrix is NULL in abs_op\n");
    return NULL;
  }

  Matrix* result = create_matrix(mat_->dims, mat_->shape);
  int size = get_matrix_size(mat_);

  for (int i = 0; i < size; i++) {
    result->data[i] = fabs(mat_->data[i]); }
  return result;
}

Matrix* exp_op(Matrix* mat_) {
  if (!mat_) {
    fprintf(stderr, "Input matrix is NULL in exp_op\n");
    return NULL;
  }

  Matrix* res_ = create_matrix(mat_->dims, mat_->shape);
  int size = get_matrix_size(mat_);

  for (int i = 0; i < size; i++) {
    res_->data[i] = exp(mat_->data[i]);
  }

  return res_;
}

Matrix* sqrt_op(Matrix* mat_) {
  if (!mat_) {
    fprintf(stderr, "Input matrix is NULL in sqrt_op\n");
    return NULL;
  }

  Matrix* res_ = create_matrix(mat_->dims, mat_->shape);
  int size = get_matrix_size(mat_);

  for (int i = 0; i < size; i++) {
    res_->data[i] = sqrt(mat_->data[i]);
  }

  return res_;
}


Matrix* loge_op(Matrix* m) {
  if (!m) {
    fprintf(stderr, "Input matrix is NULL in loge_op\n");
    return NULL;
  }

  Matrix* res_ = create_matrix(m->dims, m->shape);
  int size = get_matrix_size(m);

  for (int i = 0; i < size; i++) {
    res_->data[i] = log(m->data[i]);
  }

  return res_;
}


Matrix* log10_op(Matrix* m) {
  if (!m) {
    fprintf(stderr, "Input matrix is NULL in log10_op\n");
    return NULL;
  }

  Matrix* res_ = create_matrix(m->dims, m->shape);
  int size = get_matrix_size(m);

  for (int i = 0; i < size; i++) {
    res_->data[i] = log10(m->data[i]);
  }

  return res_;
}

Matrix* log2_op(Matrix* m) {
  if (!m) {
    fprintf(stderr, "Input matrix is NULL in log2_op\n");
    return NULL;
  }

  Matrix* res_ = create_matrix(m->dims, m->shape);
  int size = get_matrix_size(m);

  for (int i = 0; i < size; i++) {
    res_->data[i] = log2(m->data[i]);
  }

  return res_;
}


Matrix* power_op(Matrix* m, double p) {
  if (!m) {
    fprintf(stderr, "Input matrix is NULL in power_op\n");
  }

  Matrix *res = create_matrix(m->dims, m->shape);
  int size = get_matrix_size(res);

  for (int i = 0; i < size; i++) {
    res->data[i] = pow(m->data[i], p);
  }
  return res;
}



double sum1d(Matrix * m) {
  assert (m != NULL);
  assert (m->shape[0] == 1 || m->shape[1] == 1);
  double sum = 0;
  int m_size = get_matrix_size(m);

  for (int i = 0; i < m_size; i++) {
    sum += m->data[i];
  }
  return sum;
}

double prod1d(Matrix * m) {
  assert (m != NULL);
  assert (m->shape[0] == 1 || m->shape[1] == 1);
  double prod = 1;
  int m_size = get_matrix_size(m);

  for (int i = 0; i < m_size; i++) {
    prod += m->data[i];
  }
  return prod;
}

Matrix* sum2d(Matrix* m, int axis) {
  assert(m->shape[0] > 1 && m->shape[1] > 1);

  int nrows = (axis == 0) ? m->shape[0] : 1;
  int ncols = (axis == 1) ? m->shape[1] : 1;
  int new_shape[2] = {nrows, ncols};

  Matrix* sum_2dmat = create_matrix(2, new_shape);
  if (axis == 0) {
    for (int i = 0; i < m->shape[0]; i++) {
      double row_sum = 0.0;
      for (int j = 0; j < m->shape[1]; j++) {
        row_sum += m->data[i * m->shape[1] + j];
      }
      sum_2dmat->data[i] = row_sum;
    }
  }
  else if (axis == 1) {
    for (int j = 0; j < m->shape[1]; j++) {
      double col_sum = 0.0;
      for (int i = 0; i < m->shape[0]; i++) {
        col_sum += m->data[i * m->shape[1] + j];
      }
      sum_2dmat->data[j] = col_sum;
    }
  }
  return sum_2dmat;
}


Matrix* prod2d(Matrix* m, int axis) {
  assert(m->shape[0] > 1 && m->shape[1] > 1);

  int nrows = (axis == 0) ? m->shape[0] : 1;
  int ncols = (axis == 1) ? m->shape[1] : 1;
  int new_shape[2] = {nrows, ncols};

  Matrix* prod_2dmat = create_matrix(2, new_shape);
  if (axis == 0) {
    for (int i = 0; i < m->shape[0]; i++) {
      double row_prod = 1.0;
      for (int j = 0; j < m->shape[1]; j++) {
        row_prod *= m->data[i * m->shape[1] + j];
      }
      prod_2dmat->data[i] = row_prod;
    }
  }
  else if (axis == 1) {
    for (int j = 0; j < m->shape[1]; j++) {
      double col_prod = 1.0;
      for (int i = 0; i < m->shape[0]; i++) {
        col_prod *= m->data[i * m->shape[1] + j];
      }
      prod_2dmat->data[j] = col_prod;
    }
  }
  return prod_2dmat;
}


Matrix* const_mul(Matrix* mat_, const double k) {
  if (!mat_) {
    fprintf(stderr, "Input matrix is NULL in const_mul\n");
    return NULL;
  }

  Matrix* ret_ = create_matrix(mat_->dims, mat_->shape);
  int size = get_matrix_size(mat_);


  for (int i = 0; i < size; i++) {
    ret_->data[i] = mat_->data[i] * k;
  }

  return ret_;
}


Matrix* const_add(Matrix* mat_, const double k) {
  if (!mat_) {
    fprintf(stderr, "Input matrix is NULL in const_add\n");
    return NULL;
  }

  Matrix* ret_ = create_matrix(mat_->dims, mat_->shape);
  int size = get_matrix_size(mat_);

  for (int i = 0; i < size; i++) {
    ret_->data[i] = mat_->data[i] + k;
  }

  return ret_;
}

Matrix* const_sub(Matrix* mat_, const double k) {
  if (!mat_) {
    fprintf(stderr, "Input matrix is NULL in const_sub\n");
    return NULL;
  }

  Matrix* ret_ = create_matrix(mat_->dims, mat_->shape);
  int size = get_matrix_size(mat_);

  for (int i = 0; i < size; i++) {
    ret_->data[i] = mat_->data[i] - k;
  }

  return ret_;
}


Matrix* const_div(Matrix* mat_, const double k) {
  assert (k > 0);

  if (!mat_) {
    fprintf(stderr, "Input matrix is NULL in const_add\n");
    return NULL;
  }

  Matrix* ret_ = create_matrix(mat_->dims, mat_->shape);
  int size = get_matrix_size(mat_);

  for (int i = 0; i < size; i++) {
    ret_->data[i] = mat_->data[i] / k;
  }

  return ret_;
}

Matrix* rad2deg(Matrix* mat) {
  assert (mat != NULL);

  Matrix* deg_mat = create_matrix(mat->dims, mat->shape);

  int size = get_matrix_size(mat);

  for (int i = 0; i < size; i++) {
    deg_mat->data[i] = (180. / M_PI) * mat->data[i];
  }

  return deg_mat;
}


Matrix* deg2rad(Matrix* mat) {
  assert (mat != NULL);

  Matrix* rad_mat = create_matrix(mat->dims, mat->shape);

  int size = get_matrix_size(mat);

  for (int i = 0; i < size; i++) {
    rad_mat->data[i] = (M_PI / 180.) * mat->data[i];
  }

  return rad_mat;
}



Matrix* transpose2d(Matrix* mat_) {
  int tshape[] = {mat_->shape[1], mat_->shape[0]};

  Matrix* ret_ = create_matrix(mat_->dims, tshape);

  int size = get_matrix_size(mat_);

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
  int size = get_matrix_size(m);

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
  int size = get_matrix_size(m);

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

  int size = get_matrix_size(m);

  for (int i = 0; i < size; i++) {
    rev_m->data[i] = m->data[size - i - 1];
  }
  return rev_m;
}

void swap2d_rows(Matrix* m, int id_rows1, int id_rows2) {
  assert(m != NULL);
  assert(m->shape[0] > 1 && m->shape[1] > 1);

  assert(id_rows1 >= 0 && id_rows1 < m->shape[0]);
  assert(id_rows2 >= 0 && id_rows2 < m->shape[0]);

  int ncols = m->shape[1];

  Matrix* temp = create_matrix(2, (int[]){1,ncols});
  if (temp == NULL) {
    return;
  }

  for (int j = 0; j < ncols; j++) {
    temp->data[j] = m->data[id_rows1 * ncols + j];
  }

  for (int j = 0; j < ncols; j++) {
    m->data[id_rows1 * ncols + j] = m->data[id_rows2 * ncols + j];
  }

  for (int j = 0; j < ncols; j++) {
    m->data[id_rows2 * ncols + j] = temp->data[j];
  }

  free(temp);
}

void swap2d_cols(Matrix* m, int id_cols1, int id_cols2) {
  assert(m != NULL);
  assert(m->shape[0] > 1 && m->shape[1] > 1);

  assert(id_cols1 >= 0 && id_cols1 < m->shape[1]);
  assert(id_cols2 >= 0 && id_cols2 < m->shape[1]);

  int nrows = m->shape[0];
  int ncols = m->shape[1];

  Matrix* temp = create_matrix(2, (int[]){nrows, 1});
  if (temp == NULL) {
    return;
  }

  for (int j = 0; j < nrows; j++) {
    temp->data[j] = m->data[id_cols1 + ncols * j];
  }

  for (int j = 0; j < nrows; j++) {
    m->data[id_cols1 + ncols * j] = m->data[id_cols2 + ncols * j];
  }

  for (int j = 0; j < nrows; j++) {
    m->data[id_cols2 + ncols * j] = temp->data[j];
  }

  free(temp);
}


Matrix* get_matrix_diag(Matrix* m) {
  assert (is_squared_matrix(m));
  int size = m->shape[0];
  Matrix* diag_element_array = create_matrix(2, (int[]){1, size});

  for (int i = 0; i < size; i++) {
    diag_element_array->data[i] = m->data[(size + 1) * i];
  }
  return diag_element_array;
}








