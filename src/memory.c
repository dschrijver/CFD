#include <stdlib.h>

#include "../include/datatypes.h"

void allocate(double **scalar, ParamBag *params)
{
    // Scalar:
    // N physical nodes: 0.5, 1.5, ... N-0.5
    // 2 ghost nodes:    -0.5, N+0.5
    
    // U, V, W:
    // N+1 physical nodes: 0, 1, ..., N
    // 1 ghost node: -1 
    // Due to periodicity / noslip / freeslip, value at 0 = value at N
    // So treat N like a ghost node carrying the value of 0. 

    // Means scalar and vector fields have same memory structure and communication. 

    const int N = (params->NX_proc + 2) * (params->NY_proc + 2) * (params->NZ_proc + 2);

    *scalar = (double *)malloc(N * sizeof(double));
}