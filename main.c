#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "params.h"
#include "include/initialize.h"
#include "include/memory.h"
#include "include/utils.h"
#include "include/output.h"
#include "include/communicate.h"
#include "include/poisson.h"
#include "definitions.h"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    ParamBag *params = (ParamBag *)malloc(sizeof(ParamBag));
    set_params(params);

    initialize_MPI(params);
    initialize_HDF5(params);

    srand(time(NULL) + params->process_rank);

    double *pressure;
    allocate_scalar(&pressure, params);

    double *rhs;
    allocate_scalar(&rhs, params);

    // Used for Poisson solver
    double *r, *p, *q;
    allocate_scalar(&r, params);
    allocate_scalar(&p, params);
    allocate_scalar(&q, params);

    params->time = 0.0;
    params->n_output = 0;

    UNPACK_GRID

#ifdef INI_RHS_NOISE
    double avg_local = 0.0;
    FOR_DOMAIN
    {
        rhs[IDX(i, j, k)] = (double)rand() / (((double)RAND_MAX) + 1.0);
        avg_local += rhs[IDX(i, j, k)];
    }
    double avg;
    MPI_Allreduce(&avg_local, &avg, 1, MPI_DOUBLE, MPI_SUM, params->comm_cart);
    avg /= NX*NY*NZ;
    FOR_DOMAIN
    {
        rhs[IDX(i, j, k)] -= avg;
    }
#else
    FOR_DOMAIN
    {
        const double x = PHYS(i);
        const double y = PHYS(j);
        const double z = PHYS(k);
        rhs[IDX(i, j, k)] = sin(2.0 * PI * x / LX) * sin(2.0 * PI * y / LY) * sin(2.0 * PI * z / LZ)
                    + 0.5 * sin(10.0 * 2.0 * PI * x / LX) * sin(10.0 * 2.0 * PI * y / LY) * sin(10.0 * 2.0 * PI * z / LZ)
                    + 0.1 * sin(20.0 * 2.0 * PI * x / LX) * sin(20.0 * 2.0 * PI * y / LY) * sin(20.0 * 2.0 * PI * z / LZ);
    }
#endif

    output_scalar(rhs, "rhs", params);

    // Initial guess
    FOR_DOMAIN
    {
        pressure[IDX(i, j, k)] = 0.0;
    }
    communicate_scalar(pressure, params);

    double residual;
    int iters;
    poisson_solve(pressure, rhs, r, p, q, &residual, &iters, params);

    if (params->process_rank == 0)
        printf("Iterations: %d, residual: %.3e\n", iters, residual);

    output_scalar(pressure, "pressure", params);

    free(pressure);
    free(params);
    MPI_Finalize();
    return 0;
}