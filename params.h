#ifndef PARAMS_H
#define PARAMS_H

#include "include/datatypes.h"

static inline void set_params(ParamBag *params)
{
    params->NX = 1;
    params->NY = 100;
    params->NZ = 1;
    params->dx = 0.01;

    params->t_final = 1.0;
    params->dt = 1e-4;

    params->dt_output = 1e-2;
    params->dt_log = 1e-3;

    params->Re = 100.0;
    params->Fx = 1.0;
    params->Fy = 0.0;
    params->Fz = 0.0;

    params->poisson_epsilon = 1e-8;
}

#endif