#ifndef PARAMS_H
#define PARAMS_H

#include "include/datatypes.h"

static inline void set_params(ParamBag *params)
{
    params->NX = 64;
    params->NY = 64;
    params->NZ = 1;

    double Lx = 1.0;
    params->dx = Lx / (double)params->NX;

    params->t_final = 50.0;
    params->dt = 0.005;

    params->dt_output = 0.5;
    params->dt_log = 0.05;

    params->Re = 100.0;
    params->Fx = 0.0;
    params->Fy = 0.0;
    params->Fz = 0.0;

    params->poisson_epsilon = 1e-6;
}

#endif