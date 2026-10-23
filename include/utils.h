#ifndef UTILS_H
#define UTILS_H

#define PI 3.14159265358979323846

#define IDX(i, j, k) ((NY_proc + 2) * (NZ_proc + 2) * (i + 1 - i_start) + (NZ_proc + 2) * (j + 1 - j_start) + (k + 1 - k_start))

#define FOR_DOMAIN \
    for (int i = i_start; i < i_end; i++) \
        for (int j = j_start; j < j_end; j++) \
            for (int k = k_start; k < k_end; k++)

#define UNPACK_GRID \
    const int i_start = params->i_start; \
    const int i_end = params->i_end; \
    const int j_start = params->j_start; \
    const int j_end = params->j_end; \
    const int k_start = params->k_start; \
    const int k_end = params->k_end; \
    const int NX = params->NX; \
    const int NY = params->NY; \
    const int NZ = params->NZ; \
    const int NX_proc = params->NX_proc; \
    const int NY_proc = params->NY_proc; \
    const int NZ_proc = params->NZ_proc; \
    const double LX = params->LX; \
    const double LY = params->LY; \
    const double LZ = params->LZ; \
    const double h = params->h; \
    (void)i_start; (void)i_end; (void)j_start; (void)j_end; (void)k_start; (void)k_end; (void)NX; (void)NY; (void)NZ; (void)NX_proc; (void)NY_proc; (void)NZ_proc; (void)LX; (void)LY; (void)LZ; (void)h;

#define PHYS(i) (((double)i + 0.5)*h)

#endif