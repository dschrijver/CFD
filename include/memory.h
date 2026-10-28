#ifndef MEMORY_H
#define MEMORY_H

#include "datatypes.h"

void allocate(double **field, ParamBag *params);
void allocate_poisson(PoissonBag *poisson, ParamBag *params);
void allocate_buffers(BufferBag *buffers, ParamBag *params);
void free_buffers(BufferBag *buffers);
void free_poisson(PoissonBag *poisson);

#endif