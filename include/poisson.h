#ifndef POISSON_H
#define POISSON_H

#include "datatypes.h"

void poisson_solve(double *scalar, double *rhs, double *residual, int *iters, PoissonBag *poisson, BufferBag *buffers, ParamBag *params);

#endif