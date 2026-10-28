#include <stdlib.h>
#include <string.h>

#include "../include/datatypes.h"

void allocate(double **field, ParamBag *params)
{
    // Scalar:
    // N physical nodes: 0.5, 1.5, ... N-0.5
    // 2 ghost nodes:    -0.5, N+0.5
    
    // U, V, W:
    // N+1 physical nodes: 0, 1, ..., N
    // 1 ghost node: -1 
    // Due to periodicity / noslip / freeslip, value at 0 = value at N
    // So treat N like a ghost node carrying the value of 0. 

    // Means scalar and vector fields have same memory structure and communication. 

    const int N = (params->NX_proc + 2) * (params->NY_proc + 2) * (params->NZ_proc + 2);

    *field = (double *)malloc(N * sizeof(double));
}

void allocate_buffers(BufferBag *buffers, ParamBag *params)
{
    buffers->buffer_number_x = params->NY_proc * params->NZ_proc;
    buffers->buffer_number_y = params->NX_proc * params->NZ_proc;
    buffers->buffer_number_z = params->NX_proc * params->NY_proc;

    buffers->send_buffer_x = (double *)malloc(buffers->buffer_number_x * sizeof(double));
    buffers->send_buffer_y = (double *)malloc(buffers->buffer_number_y * sizeof(double));
    buffers->send_buffer_z = (double *)malloc(buffers->buffer_number_z * sizeof(double));

    buffers->recv_buffer_x = (double *)malloc(buffers->buffer_number_x * sizeof(double));
    buffers->recv_buffer_y = (double *)malloc(buffers->buffer_number_y * sizeof(double));
    buffers->recv_buffer_z = (double *)malloc(buffers->buffer_number_z * sizeof(double));
}

void allocate_poisson(PoissonBag *poisson, ParamBag *params)
{
    const int N = (params->NX_proc + 2) * (params->NY_proc + 2) * (params->NZ_proc + 2);

    poisson->r = (double *)malloc(N * sizeof(double));
    poisson->p = (double *)malloc(N * sizeof(double));
    poisson->q = (double *)malloc(N * sizeof(double));
}

void free_poisson(PoissonBag *poisson)
{
    free(poisson->r);
    free(poisson->p);
    free(poisson->q);
    free(poisson);
}

void free_buffers(BufferBag *buffers)
{
    free(buffers->send_buffer_x);
    free(buffers->send_buffer_y);
    free(buffers->send_buffer_z);
    free(buffers->recv_buffer_x);
    free(buffers->recv_buffer_y);
    free(buffers->recv_buffer_z);
    free(buffers);
}