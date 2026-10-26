#ifndef INITIALIZE_H
#define INITIALIZE_H

#include "datatypes.h"

void initialize_MPI(ParamBag *params);
void initialize_HDF5(ParamBag *params);
void initialize_fields(double *u, double *v, double *w, double *pressure, ParamBag *params);

#endif