#include <stdlib.h>

#include "../include/datatypes.h"

void allocate_scalar(double **scalar, ParamBag *params)
{
    const int N = (params->NX_proc + 2) * (params->NY_proc + 2) * (params->NZ_proc + 2);

    *scalar = (double *)malloc(N * sizeof(double));
}