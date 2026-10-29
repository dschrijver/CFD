#include <math.h>

#include "../include/datatypes.h"
#include "../include/utils.h"
#include "../include/communicate.h"
#include "../include/bc.h"

// Solves nabla^2 scalar = rhs / dt with conjugate gradient method. 
void poisson_solve(double *scalar, double *rhs, double *residual, int *iters, PoissonBag *poisson, BufferBag *buffers, ParamBag *params)
{
    UNPACK_GRID
    UNPACK_POISSON

    const double epsilon = params->poisson_epsilon;
    const double dt = params->dt;

    double rr_initial_local = 0.0;
    double residual_local = 0.0;
    FOR_DOMAIN
    {
        const int idx0 = IDX(i, j, k);

        // At no-slip boundary conditions, zero normal gradients. 
        // At dirichlet boundary, set to specific value.
        const double d2dx2 = (scalar[IDX(i-1,j,k)] - 2.0*scalar[idx0] + scalar[IDX(i+1,j,k)]) / (dx*dx);
        const double d2dy2 = (scalar[IDX(i,j-1,k)] - 2.0*scalar[idx0] + scalar[IDX(i,j+1,k)]) / (dx*dx);
        const double d2dz2 = (scalar[IDX(i,j,k-1)] - 2.0*scalar[idx0] + scalar[IDX(i,j,k+1)]) / (dx*dx);

        const double lhs = d2dx2 + d2dy2 + d2dz2;

        r[idx0] = rhs[idx0] / dt - lhs;

        rr_initial_local += r[idx0] * r[idx0];
        const double residual_i = fabs(r[idx0]);
        if (residual_i > residual_local) residual_local = residual_i;

        p[idx0] = r[idx0];
    }
    double rr_initial;
    MPI_Allreduce(&rr_initial_local, &rr_initial, 1, MPI_DOUBLE, MPI_SUM, params->comm_cart);
    
    MPI_Allreduce(&residual_local, residual, 1, MPI_DOUBLE, MPI_MAX, params->comm_cart);
    *residual *= dt;

    if (*residual < epsilon)
    {
        *iters = 0;
        return;
    }

    *iters = 0;
    double rr_previous = rr_initial;
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
        const double alpha = rr_previous / alpha_denom;

        FOR_DOMAIN
        {
            const int idx0 = IDX(i, j, k);

            scalar[idx0] += alpha*p[idx0];
            r[idx0] -= alpha*q[idx0];
        }

        double rr_local = 0.0;
        double residual_local = 0.0;
        FOR_DOMAIN
        {
            const int idx0 = IDX(i, j, k);

            rr_local += r[idx0] * r[idx0];

            const double residual_i = fabs(r[idx0]);
            if (residual_i > residual_local) residual_local = residual_i;
        }
        double rr;
        MPI_Allreduce(&rr_local, &rr, 1, MPI_DOUBLE, MPI_SUM, params->comm_cart);
        MPI_Allreduce(&residual_local, residual, 1, MPI_DOUBLE, MPI_MAX, params->comm_cart);
        *residual *= params->dt;

        *iters += 1;

        if (*residual < epsilon)
            break;
            
        if (rr / rr_initial < epsilon)
            break;

        const double beta = rr / rr_previous;

        FOR_DOMAIN
        {
            const int idx0 = IDX(i, j, k);

            p[idx0] = r[idx0] + beta*p[idx0];
        }

        rr_previous = rr;
    }
}