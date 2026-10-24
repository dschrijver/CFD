#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/utils.h"

void pressure_bc(double *pressure, ParamBag *params)
{
    UNPACK_GRID
    
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

    // -----------------------------------------------------------------------------
    // DIRICHLET BOUNDARIES
    // Implemented by setting average of boundary and ghost node to boundary value.
    // -----------------------------------------------------------------------------
#ifdef LEFT_PRESSURE
    if (i_start == 0) for (int j = j_start; j < j_end; j++) for (int k = k_start; k < k_end; k++)
    {
        pressure[IDX(i_start - 1, j, k)] = 2.0 * LEFT_PRESSURE - pressure[IDX(i_start, j, k)];
    }
#endif
#ifdef RIGHT_PRESSURE
    if (i_end == NX) for (int j = j_start; j < j_end; j++) for (int k = k_start; k < k_end; k++)
    {
        pressure[IDX(i_end, j, k)] = 2.0 * RIGHT_PRESSURE - pressure[IDX(i_end - 1, j, k)];
    }
#endif
#ifdef BOTTOM_PRESSURE
    if (j_start == 0) for (int i = i_start; i < i_end; i++) for (int k = k_start; k < k_end; k++)
    {
        pressure[IDX(i, j_start - 1, k)] = 2.0 * BOTTOM_PRESSURE - pressure[IDX(i, j_start, k)];
    }
#endif
#ifdef TOP_PRESSURE
    if (j_end == NY) for (int i = i_start; i < i_end; i++) for (int k = k_start; k < k_end; k++)
    {
        pressure[IDX(i, j_end, k)] = 2.0 * TOP_PRESSURE - pressure[IDX(i, j_end - 1, k)];
    }
#endif
#ifdef BACK_PRESSURE
    if (k_start == 0) for (int i = i_start; i < i_end; i++) for (int j = j_start; j < j_end; j++)
    {
        pressure[IDX(i, j, k_start - 1)] = 2.0 * BACK_PRESSURE - pressure[IDX(i, j, k_start)];
    }
#endif
#ifdef FRONT_PRESSURE
    if (k_end == NZ) for (int i = i_start; i < i_end; i++) for (int j = j_start; j < j_end; j++)
    {
        pressure[IDX(i, j, k_end)] = 2.0 * FRONT_PRESSURE - pressure[IDX(i, j, k_end - 1)];
    }
#endif
}

void pressure_poisson_bc(double *p, ParamBag *params)
{
    UNPACK_GRID
    
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

    // -----------------------------------------------------------------------------
    // DIRICHLET BOUNDARIES
    // Implemented by setting average of boundary and ghost node to zero.
    // -----------------------------------------------------------------------------
#ifdef LEFT_PRESSURE
    if (i_start == 0) for (int j = j_start; j < j_end; j++) for (int k = k_start; k < k_end; k++)
    {
        p[IDX(i_start - 1, j, k)] = -p[IDX(i_start, j, k)];
    }
#endif
#ifdef RIGHT_PRESSURE
    if (i_end == NX) for (int j = j_start; j < j_end; j++) for (int k = k_start; k < k_end; k++)
    {
        p[IDX(i_end, j, k)] = -p[IDX(i_end - 1, j, k)];
    }
#endif
#ifdef BOTTOM_PRESSURE
    if (j_start == 0) for (int i = i_start; i < i_end; i++) for (int k = k_start; k < k_end; k++)
    {
        p[IDX(i, j_start - 1, k)] = -p[IDX(i, j_start, k)];
    }
#endif
#ifdef TOP_PRESSURE
    if (j_end == NY) for (int i = i_start; i < i_end; i++) for (int k = k_start; k < k_end; k++)
    {
        p[IDX(i, j_end, k)] = -p[IDX(i, j_end - 1, k)];
    }
#endif
#ifdef BACK_PRESSURE
    if (k_start == 0) for (int i = i_start; i < i_end; i++) for (int j = j_start; j < j_end; j++)
    {
        p[IDX(i, j, k_start - 1)] = -p[IDX(i, j, k_start)];
    }
#endif
#ifdef FRONT_PRESSURE
    if (k_end == NZ) for (int i = i_start; i < i_end; i++) for (int j = j_start; j < j_end; j++)
    {
        p[IDX(i, j, k_end)] = -p[IDX(i, j, k_end - 1)];
    }
#endif
}