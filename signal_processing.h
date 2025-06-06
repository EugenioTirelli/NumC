#ifndef SIGNAL_PROCESSING_H
#define SIGNAL_PROCESSING_H

#include "basics.h"
#include "cbasics.h"
#include <complex.h>

CMatrix *dft1d(Matrix *);
Matrix *convolve1d(Matrix *, Matrix *, int);
#endif
