#ifndef BC_H
#define BC_H

#include "datatypes.h"

void velocity_bc(double *u, double *v, double *w, ParamBag *params);
void pressure_bc(double *pressure, ParamBag *params);
void pressure_poisson_bc(double *p, ParamBag *params);

#endif