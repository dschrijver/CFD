#ifndef PARAMS_H
#define PARAMS_H

#include "include/datatypes.h"

static inline void set_params(ParamBag *params)
{
    params->NX = 100;
    params->NY = 100;
    params->NZ = 100;

    params->h = 0.01;

    params->poisson_epsilon = 1e-8;
}

#endif