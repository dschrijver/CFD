#ifndef PARAMS_H
#define PARAMS_H

#include "../include/datatypes.h"

static inline void set_params(ParamBag *params)
{
    params->NX = 100;
    params->NY = 100;
    params->NZ = 100;

    double Lx = 1.0;
    params->dx = Lx / (double)params->NX;

    params->dt = 1.0;

    params->poisson_epsilon = 1e-6;
}

#endif