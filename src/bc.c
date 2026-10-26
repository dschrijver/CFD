#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/utils.h"

void velocity_bc(double *u, double *v, double *w, ParamBag *params)
{
    UNPACK_GRID

    (void)u; (void)v; (void)w;

    // -----------------------------------------------------------------------------
    // NO-SLIP WALLS
    // -----------------------------------------------------------------------------
#ifdef LEFT_NOSLIP
    if (i_start == 0) for (int j = j_start; j < j_end; j++) for (int k = k_start; k < k_end; k++)
    {
        u[IDX(i_start, j, k)] = 0.0;
        v[IDX(i_start - 1, j, k)] = -v[IDX(i_start, j, k)];
        w[IDX(i_start - 1, j, k)] = -w[IDX(i_start, j, k)];
    }
#endif
#ifdef RIGHT_NOSLIP
    if (i_end == NX) for (int j = j_start; j < j_end; j++) for (int k = k_start; k < k_end; k++)
    {
        u[IDX(i_end, j, k)] = 0.0;
        v[IDX(i_end, j, k)] = -v[IDX(i_end - 1, j, k)];
        w[IDX(i_end, j, k)] = -w[IDX(i_end - 1, j, k)];
    }
#endif
#ifdef BOTTOM_NOSLIP
    if (j_start == 0) for (int i = i_start; i < i_end; i++) for (int k = k_start; k < k_end; k++)
    {
        u[IDX(i, j_start - 1, k)] = -u[IDX(i, j_start, k)];
        v[IDX(i, j_start, k)] = 0.0;
        w[IDX(i, j_start - 1, k)] = -w[IDX(i, j_start, k)];
    }
#endif
#ifdef TOP_NOSLIP
    if (j_end == NY) for (int i = i_start; i < i_end; i++) for (int k = k_start; k < k_end; k++)
    {
        u[IDX(i, j_end, k)] = -u[IDX(i, j_end - 1, k)];
        v[IDX(i, j_end, k)] = 0.0;
        w[IDX(i, j_end, k)] = -w[IDX(i, j_end - 1, k)];
    }
#endif
#ifdef BACK_NOSLIP
    if (k_start == 0) for (int i = i_start; i < i_end; i++) for (int j = j_start; j < j_end; j++)
    {
        u[IDX(i, j, k_start - 1)] = -u[IDX(i, j, k_start)];
        v[IDX(i, j, k_start - 1)] = -v[IDX(i, j, k_start)];
        w[IDX(i, j, k_start)] = 0.0;
    }
#endif
#ifdef FRONT_NOSLIP
    if (k_end == NZ) for (int i = i_start; i < i_end; i++) for (int j = j_start; j < j_end; j++)
    {
        u[IDX(i, j, k_end)] = -u[IDX(i, j, k_end - 1)];
        v[IDX(i, j, k_end)] = -v[IDX(i, j, k_end - 1)];
        w[IDX(i, j, k_end)] = 0.0;
    }
#endif
}

void velocity_bc_after(double *u, double *v, double *w, ParamBag *params)
{
    UNPACK_GRID

    (void)u; (void)v; (void)w;

    // -----------------------------------------------------------------------------
    // NO-SLIP WALLS
    // -----------------------------------------------------------------------------
#ifdef LEFT_NOSLIP
    if (i_start == 0) for (int j = j_start; j < j_end; j++) for (int k = k_start; k < k_end; k++)
    {
        u[IDX(i_start, j, k)] = 0.0;
    }
#endif
#ifdef BOTTOM_NOSLIP
    if (j_start == 0) for (int i = i_start; i < i_end; i++) for (int k = k_start; k < k_end; k++)
    {
        v[IDX(i, j_start, k)] = 0.0;
    }
#endif
#ifdef BACK_NOSLIP
    if (k_start == 0) for (int i = i_start; i < i_end; i++) for (int j = j_start; j < j_end; j++)
    {
        w[IDX(i, j, k_start)] = 0.0;
    }
#endif
}

void pressure_bc(double *pressure, ParamBag *params)
{
    UNPACK_GRID
    
    (void)pressure;

    // -----------------------------------------------------------------------------
    // NO-SLIP WALLS
    // Implemented by setting normal pressure gradients to zero. 
    // -----------------------------------------------------------------------------
#ifdef LEFT_NOSLIP
    if (i_start == 0) for (int j = j_start; j < j_end; j++) for (int k = k_start; k < k_end; k++)
    {
        pressure[IDX(i_start - 1, j, k)] = pressure[IDX(i_start, j, k)];
    }
#endif
#ifdef RIGHT_NOSLIP
    if (i_end == NX) for (int j = j_start; j < j_end; j++) for (int k = k_start; k < k_end; k++)
    {
        pressure[IDX(i_end, j, k)] = pressure[IDX(i_end - 1, j, k)];
    }
#endif
#ifdef BOTTOM_NOSLIP
    if (j_start == 0) for (int i = i_start; i < i_end; i++) for (int k = k_start; k < k_end; k++)
    {
        pressure[IDX(i, j_start - 1, k)] = pressure[IDX(i, j_start, k)];
    }
#endif
#ifdef TOP_NOSLIP
    if (j_end == NY) for (int i = i_start; i < i_end; i++) for (int k = k_start; k < k_end; k++)
    {
        pressure[IDX(i, j_end, k)] = pressure[IDX(i, j_end - 1, k)];
    }
#endif
#ifdef BACK_NOSLIP
    if (k_start == 0) for (int i = i_start; i < i_end; i++) for (int j = j_start; j < j_end; j++)
    {
        pressure[IDX(i, j, k_start - 1)] = pressure[IDX(i, j, k_start)];
    }
#endif
#ifdef FRONT_NOSLIP
    if (k_end == NZ) for (int i = i_start; i < i_end; i++) for (int j = j_start; j < j_end; j++)
    {
        pressure[IDX(i, j, k_end)] = pressure[IDX(i, j, k_end - 1)];
    }
#endif
}

void pressure_poisson_bc(double *p, ParamBag *params)
{
    UNPACK_GRID

    (void)p;
    
    // -----------------------------------------------------------------------------
    // NO-SLIP WALLS
    // Implemented by setting normal p gradients to zero. 
    // -----------------------------------------------------------------------------
#ifdef LEFT_NOSLIP
    if (i_start == 0) for (int j = j_start; j < j_end; j++) for (int k = k_start; k < k_end; k++)
    {
        p[IDX(i_start - 1, j, k)] = p[IDX(i_start, j, k)];
    }
#endif
#ifdef RIGHT_NOSLIP
    if (i_end == NX) for (int j = j_start; j < j_end; j++) for (int k = k_start; k < k_end; k++)
    {
        p[IDX(i_end, j, k)] = p[IDX(i_end - 1, j, k)];
    }
#endif
#ifdef BOTTOM_NOSLIP
    if (j_start == 0) for (int i = i_start; i < i_end; i++) for (int k = k_start; k < k_end; k++)
    {
        p[IDX(i, j_start - 1, k)] = p[IDX(i, j_start, k)];
    }
#endif
#ifdef TOP_NOSLIP
    if (j_end == NY) for (int i = i_start; i < i_end; i++) for (int k = k_start; k < k_end; k++)
    {
        p[IDX(i, j_end, k)] = p[IDX(i, j_end - 1, k)];
    }
#endif
#ifdef BACK_NOSLIP
    if (k_start == 0) for (int i = i_start; i < i_end; i++) for (int j = j_start; j < j_end; j++)
    {
        p[IDX(i, j, k_start - 1)] = p[IDX(i, j, k_start)];
    }
#endif
#ifdef FRONT_NOSLIP
    if (k_end == NZ) for (int i = i_start; i < i_end; i++) for (int j = j_start; j < j_end; j++)
    {
        p[IDX(i, j, k_end)] = p[IDX(i, j, k_end - 1)];
    }
#endif
}