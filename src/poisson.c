#include "../include/datatypes.h"
#include "../include/utils.h"
#include "../include/communicate.h"

// Solves nabla^2 scalar = rhs with conjugate gradient method. 
void poisson_solve(double *scalar, double *rhs, double *r, double *p, double *q, double *residual, int *iters, ParamBag *params)
{
    UNPACK_GRID

    const double epsilon = params->poisson_epsilon;
    
    double residual_initial_local = 0.0;
    FOR_DOMAIN
    {
        const int idx0 = IDX(i, j, k);

        // Periodic boundary conditions assumed. 
        // At no-slip boundary conditions, zero normal gradients. 
        // At dirichlet boundary, set to specific value.
        const double d2dx2 = (scalar[IDX(i-1,j,k)] - 2.0*scalar[idx0] + scalar[IDX(i+1,j,k)]) / (h*h);
        const double d2dy2 = (scalar[IDX(i,j-1,k)] - 2.0*scalar[idx0] + scalar[IDX(i,j+1,k)]) / (h*h);
        const double d2dz2 = (scalar[IDX(i,j,k-1)] - 2.0*scalar[idx0] + scalar[IDX(i,j,k+1)]) / (h*h);

        const double lhs = d2dx2 + d2dy2 + d2dz2;

        r[idx0] = rhs[idx0] - lhs;
        residual_initial_local += r[idx0]*r[idx0];
        p[idx0] = r[idx0];
    }
    double residual_initial;
    MPI_Allreduce(&residual_initial_local, &residual_initial, 1, MPI_DOUBLE, MPI_SUM, params->comm_cart);

    *iters = 0;
    double residual_previous = residual_initial;
    while (*iters < 100)
    {
        communicate_scalar(p, params);

        double alpha_denom_local = 0.0;
        FOR_DOMAIN
        {
            const int idx0 = IDX(i, j, k);

            // Periodic boundary conditions assumed.
            // At no-slip boundary conditions, zero normal gradients. 
            // At dirichlet boundary, set to zero.
            // These two rules are defined by scalar_new = scalar_old + alpha p. 
            const double d2dx2 = (p[IDX(i-1,j,k)] - 2.0*p[idx0] + p[IDX(i+1,j,k)]) / (h*h);
            const double d2dy2 = (p[IDX(i,j-1,k)] - 2.0*p[idx0] + p[IDX(i,j+1,k)]) / (h*h);
            const double d2dz2 = (p[IDX(i,j,k-1)] - 2.0*p[idx0] + p[IDX(i,j,k+1)]) / (h*h);

            q[idx0] = d2dx2 + d2dy2 + d2dz2;

            alpha_denom_local += p[idx0]*q[idx0];
        }
        double alpha_denom;
        MPI_Allreduce(&alpha_denom_local, &alpha_denom, 1, MPI_DOUBLE, MPI_SUM, params->comm_cart);
        const double alpha = residual_previous / alpha_denom;

        FOR_DOMAIN
        {
            const int idx0 = IDX(i, j, k);

            scalar[idx0] += alpha*p[idx0];
            r[idx0] -= alpha*q[idx0];
        }

        double residual_local = 0.0;
        FOR_DOMAIN
        {
            const int idx0 = IDX(i, j, k);

            residual_local += r[idx0]*r[idx0];
        }
        MPI_Allreduce(&residual_local, residual, 1, MPI_DOUBLE, MPI_SUM, params->comm_cart);

        *iters += 1;

        if (*residual / residual_initial < epsilon)
        {
            break;
        }

        const double beta = *residual / residual_previous;

        FOR_DOMAIN
        {
            const int idx0 = IDX(i, j, k);

            p[idx0] = r[idx0] + beta*p[idx0];
        }
    }
}