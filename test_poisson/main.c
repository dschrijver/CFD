#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "params.h"
#include "../include/initialize.h"
#include "../include/memory.h"
#include "../include/utils.h"
#include "../include/output.h"
#include "../include/communicate.h"
#include "../include/poisson.h"
#include "definitions.h"
#include "../include/bc.h"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    ParamBag *params = (ParamBag *)malloc(sizeof(ParamBag));
    set_params(params);

    initialize_MPI(params);
    initialize_HDF5(params);

    double *pressure;
    allocate(&pressure, params);

    double *rhs;
    allocate(&rhs, params);

    // Communication buffers
    BufferBag *buffers = (BufferBag *)malloc(sizeof(BufferBag));
    allocate_buffers(buffers, params);

    // Poisson fields
    PoissonBag *poisson = (PoissonBag *)malloc(sizeof(PoissonBag));
    allocate_poisson(poisson, params);

    double residual_poisson = 0.0;
    int iters_poisson = 0;

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
#elif defined(INI_RHS_COSINE)
    FOR_DOMAIN
    {
        const double x = PHYS_CENTER(i);
        const double y = PHYS_CENTER(j);
        const double z = PHYS_CENTER(k);
        rhs[IDX(i, j, k)] = cos(PI * x / LX) * cos(PI * y / LY) * cos(PI * z / LZ)
                    + 0.5 * cos(19 * PI * x / LX) * cos(19 * PI * y / LY) * cos(19 * PI * z / LZ)
                    + 0.1 * cos(39.0 * PI * x / LX) * cos(39.0 * PI * y / LY) * cos(39.0 * PI * z / LZ);
    }
#else //     INI_RHS_SINE
    FOR_DOMAIN
    {
        const double x = PHYS_CENTER(i);
        const double y = PHYS_CENTER(j);
        const double z = PHYS_CENTER(k);
        rhs[IDX(i, j, k)] = sin(2.0 * PI * x / LX) * sin(2.0 * PI * y / LY) * sin(2.0 * PI * z / LZ)
                    + 0.5 * sin(10.0 * 2.0 * PI * x / LX) * sin(10.0 * 2.0 * PI * y / LY) * sin(10.0 * 2.0 * PI * z / LZ)
                    + 0.1 * sin(20.0 * 2.0 * PI * x / LX) * sin(20.0 * 2.0 * PI * y / LY) * sin(20.0 * 2.0 * PI * z / LZ);
    }
#endif

    output(rhs, "rhs", params);

    communicate(pressure, buffers, params);
    pressure_bc(pressure, params);
    poisson_solve(pressure, rhs, &residual_poisson, &iters_poisson, poisson, buffers, params);

    if (params->process_rank == 0)
        printf("Iters: %d, Residual: %.3e\n", iters_poisson, residual_poisson);

    output(pressure, "pressure", params);

    free_poisson(poisson);

    free_buffers(buffers);

    free(pressure);

    free(params);
    MPI_Finalize();
    return 0;
}