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
    const double dx = params->dx; \
    (void)i_start; (void)i_end; (void)j_start; (void)j_end; (void)k_start; (void)k_end; (void)NX; (void)NY; (void)NZ; (void)NX_proc; (void)NY_proc; (void)NZ_proc; (void)LX; (void)LY; (void)LZ; (void)dx;

#define UNPACK_POSITIONS \
    const double x_face = PHYS_FACE(i); \
    const double y_face = PHYS_FACE(j); \
    const double z_face = PHYS_FACE(k); \
    const double x_center = PHYS_CENTER(i); \
    const double y_center = PHYS_CENTER(j); \
    const double z_center = PHYS_CENTER(k); \
    (void)x_face; (void)y_face; (void)z_face; (void)x_center; (void)y_center; (void)z_center;

#define UNPACK_BUFFERS \
    const int buffer_number_x = buffers->buffer_number_x; \
    const int buffer_number_y = buffers->buffer_number_y; \
    const int buffer_number_z = buffers->buffer_number_z; \
    double *send_buffer_x = buffers->send_buffer_x; \
    double *send_buffer_y = buffers->send_buffer_y; \
    double *send_buffer_z = buffers->send_buffer_z; \
    double *recv_buffer_x = buffers->recv_buffer_x; \
    double *recv_buffer_y = buffers->recv_buffer_y; \
    double *recv_buffer_z = buffers->recv_buffer_z;

#define UNPACK_POISSON \
    double *r = poisson->r; \
    double *p = poisson->p; \
    double *q = poisson->q;

#define PHYS_CENTER(i) (((double)i + 0.5)*dx)
#define PHYS_FACE(i) ((double)i*dx)

#define TIME(name, functions)                                        \
    if ((params->process_rank == 0) && (params->time > params->time_log - 1.5*params->dt)) \
    {                                                                \
        printf("%-69s", name);                                       \
    }                                                                \
    MPI_Barrier(MPI_COMM_WORLD);                                     \
    start_substep = MPI_Wtime();                                     \
    functions                                                        \
        MPI_Barrier(MPI_COMM_WORLD);                                 \
    duration_substep = MPI_Wtime() - start_substep;                  \
    if ((params->process_rank == 0) && (params->time > params->time_log - 1.5*params->dt)) \
    {                                                                \
        printf("[%8.5fs]\n", duration_substep);                      \
    }

#define TIME_OUTPUT(functions)                                                                \
    if ((params->process_rank == 0) && (params->time > params->time_log - 1.5*params->dt))                          \
    {                                                                                         \
        sprintf(output_info, "> \033[0;32mOutput to data_%d.h5...\033[0m", params->n_output); \
        printf("%-80s", output_info);                                                         \
    }                                                                                         \
    MPI_Barrier(MPI_COMM_WORLD);                                                              \
    start_substep = MPI_Wtime();                                                              \
    functions                                                                                 \
        MPI_Barrier(MPI_COMM_WORLD);                                                          \
    duration_substep = MPI_Wtime() - start_substep;                                           \
    if ((params->process_rank == 0) && (params->time > params->time_log - 1.5*params->dt))                          \
    {                                                                                         \
        printf("[%8.5fs]\n", duration_substep);                                               \
    }

#endif