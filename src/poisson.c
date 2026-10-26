#include <math.h>

#include "../include/datatypes.h"
#include "../include/utils.h"
#include "../include/communicate.h"
#include "../include/bc.h"

// Solves nabla^2 scalar = rhs with conjugate gradient method. 
void poisson_solve(double *scalar, double *rhs, double *r, double *p, double *q, double *residual, int *iters, BufferBag *buffers, ParamBag *params)
{
    UNPACK_GRID

    const double epsilon = params->poisson_epsilon;
    
    double r2_initial_local = 0.0;
    FOR_DOMAIN
    {
        const int idx0 = IDX(i, j, k);

        // At no-slip boundary conditions, zero normal gradients. 
        // At dirichlet boundary, set to specific value.
        const double d2dx2 = (scalar[IDX(i-1,j,k)] - 2.0*scalar[idx0] + scalar[IDX(i+1,j,k)]) / (dx*dx);
        const double d2dy2 = (scalar[IDX(i,j-1,k)] - 2.0*scalar[idx0] + scalar[IDX(i,j+1,k)]) / (dx*dx);
        const double d2dz2 = (scalar[IDX(i,j,k-1)] - 2.0*scalar[idx0] + scalar[IDX(i,j,k+1)]) / (dx*dx);

        const double lhs = d2dx2 + d2dy2 + d2dz2;

        r[idx0] = rhs[idx0] - lhs;
        r2_initial_local += r[idx0]*r[idx0];
        p[idx0] = r[idx0];
    }
    double r2_initial;
    MPI_Allreduce(&r2_initial_local, &r2_initial, 1, MPI_DOUBLE, MPI_SUM, params->comm_cart);

    double residual_initial = params->dt * sqrt(1.0/((double)NX*NY*NZ) * r2_initial);

    if (residual_initial < epsilon)
    {
        *residual = residual_initial;
        *iters = 0;
        return;
    }

    *iters = 0;
    double r2_previous = r2_initial;
    while (*iters < 1000)
    {
        communicate(p, buffers, params);
        pressure_poisson_bc(p, params);

        double alpha_denom_local = 0.0;
        FOR_DOMAIN
        {
            const int idx0 = IDX(i, j, k);

            // At no-slip boundary conditions, zero normal gradients. 
            // At dirichlet boundary, set to zero.
            // These two rules are defined by scalar_new = scalar_old + alpha p. 
            const double d2dx2 = (p[IDX(i-1,j,k)] - 2.0*p[idx0] + p[IDX(i+1,j,k)]) / (dx*dx);
            const double d2dy2 = (p[IDX(i,j-1,k)] - 2.0*p[idx0] + p[IDX(i,j+1,k)]) / (dx*dx);
            const double d2dz2 = (p[IDX(i,j,k-1)] - 2.0*p[idx0] + p[IDX(i,j,k+1)]) / (dx*dx);

            q[idx0] = d2dx2 + d2dy2 + d2dz2;

            alpha_denom_local += p[idx0]*q[idx0];
        }
        double alpha_denom;
        MPI_Allreduce(&alpha_denom_local, &alpha_denom, 1, MPI_DOUBLE, MPI_SUM, params->comm_cart);
        const double alpha = r2_previous / alpha_denom;

        FOR_DOMAIN
        {
            const int idx0 = IDX(i, j, k);

            scalar[idx0] += alpha*p[idx0];
            r[idx0] -= alpha*q[idx0];
        }

        double r2_local = 0.0;
        FOR_DOMAIN
        {
            const int idx0 = IDX(i, j, k);

            r2_local += r[idx0]*r[idx0];
        }
        double r2;
        MPI_Allreduce(&r2_local, &r2, 1, MPI_DOUBLE, MPI_SUM, params->comm_cart);

        *iters += 1;

        *residual = params->dt * sqrt(1.0/((double)NX*NY*NZ) * r2);

        if (*residual < epsilon)
            break;

        const double beta = r2 / r2_previous;

        FOR_DOMAIN
        {
            const int idx0 = IDX(i, j, k);

            p[idx0] = r[idx0] + beta*p[idx0];
        }

        r2_previous = r2;
    }
}