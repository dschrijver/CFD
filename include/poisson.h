#ifndef POISSON_H
#define POISSON_H

#include "datatypes.h"

void poisson_solve(double *scalar, double *rhs, double *r, double *p, double *q, double *residual, int *iters, BufferBag *buffers, ParamBag *params);

#endif