#include "basics.h"
#include <time.h>
#include <assert.h>
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>


Matrix* zeros(int dims, int* shape) {
  assert (dims > 0);
  assert (shape != NULL);

  Matrix* ret_ = create_matrix(dims, shape);

  return ret_;
}

Matrix* ones(int dims, int* shape) {
  assert (dims > 0);
  assert (shape != NULL);

  int size = 1;
  for (int i = 0; i<dims; i++) {
    size *= shape[i];
  }

  Matrix* ret_ = create_matrix(dims, shape);

  for (size_t i = 0; i < size; i++)
  {
    ret_->data[i] = 1.0f;
  }

  return ret_;
}


Matrix* diag_k(int size, const double k) {
  assert(k>0);
  assert(size>0);
  int shape[] = {size, size};
  Matrix* ret_diag = create_matrix(2, shape);

  for (int i = 0; i<size; i++) {
    for (int j = 0; j<size; j++) {
      if (i == j) {
        ret_diag->data[i * size + i] = k;
      }
    }
  }
  return ret_diag;
}

Matrix* rndmat(int dims, int* shape, const double min_, const double max_) {
  assert(max_ > min_);  // Max should be greater than min
  assert(dims > 0);     // Dimensions should be positive
  assert(shape != NULL); // Shape array should not be NULL

  Matrix* rnd_matrix = create_matrix(dims, shape);
  if (!rnd_matrix || !rnd_matrix->data) {
    fprintf(stderr, "Failed to allocate memory.\n");
    return NULL;
  }

  int size = 1;
  for (int i = 0; i < dims; i++) {
    size *= shape[i];
  }

  for (int i = 0; i < size; i++) {
    double random = ((double)rand() / RAND_MAX);
    rnd_matrix->data[i] = min_ + random * (max_ - min_);
  }
  return rnd_matrix;
}

Matrix* linspacef(double start, double stop, int nums, int klast, int axis) {
  assert(stop > start);                  // Ensure valid range
  assert(nums > 1);                      // Ensure nums is at least 2
  assert(klast == 0 || klast == 1);      // klast must be 0 or 1
  assert(axis == 0 || axis == 1);        // axis must be 0 or 1

  int shape_[2];
  shape_[0] = (axis == 0) ? 1 : nums;
  shape_[1] = (axis == 0) ? nums : 1;

  Matrix *ret_ = create_matrix(2, shape_);

  double step_ = (stop - start) / (klast == 1 ? (nums - 1) : nums);

  for (int i = 0; i < nums; i++) {
    ret_->data[i] = start + i * step_;
  }

  return ret_;
}
