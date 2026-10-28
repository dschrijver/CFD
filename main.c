#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "params.h"
#include "include/initialize.h"
#include "include/memory.h"
#include "include/utils.h"
#include "include/output.h"
#include "include/communicate.h"
#include "include/poisson.h"
#include "definitions.h"
#include "include/bc.h"
#include "include/velocity.h"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    ParamBag *params = (ParamBag *)malloc(sizeof(ParamBag));
    set_params(params);

    initialize_MPI(params);
    initialize_HDF5(params);

    double *u, *v, *w;
    allocate(&u, params);
    allocate(&v, params);
    allocate(&w, params);

    double *u_predict, *v_predict, *w_predict;
    allocate(&u_predict, params);
    allocate(&v_predict, params);
    allocate(&w_predict, params);

    double *pressure;
    allocate(&pressure, params);

    double *divergence;
    allocate(&divergence, params);

    // Used for Poisson solver
    double *r_poisson, *s_poisson, *p_poisson, *q_poisson;
    allocate(&r_poisson, params);
    allocate(&s_poisson, params);
    allocate(&p_poisson, params);
    allocate(&q_poisson, params);

    // Communication buffers
    BufferBag *buffers = (BufferBag *)malloc(sizeof(BufferBag));
    allocate_buffers(buffers, params);

    // Poisson fields
    PoissonBag *poisson = (PoissonBag *)malloc(sizeof(PoissonBag));
    allocate_poisson(poisson, params);

    initialize_fields(u, v, w, pressure, divergence, params);

    params->time = 0.0;
    params->n_output = 0;

    output(u, "u", params);
    output(v, "v", params);
    output(w, "w", params);
    output(divergence, "divergence", params);
    output(pressure, "pressure", params);

    params->time_log = params->dt_log;
    params->time_output = params->dt_output;
    params->n_output = 1;

    double residual_poisson = 0.0;
    int iters_poisson = 0;

    char output_info[128];

    double start_substep = 0.0, duration_substep = 0.0;

    while (params->time < params->t_final)
    {
        // LOGGING TIME
        if (params->time > params->time_log - 1.5*params->dt)
        {
            if (params->process_rank == 0)
            {
                printf("================================================================================\n");
                printf("Time = %.3e\n", params->time);
                printf("--------------------------------------------------------------------------------\n");
            }
        }

        TIME("> Predict velocities...", 
            communicate(u, buffers, params);
            communicate(v, buffers, params);
            communicate(w, buffers, params);
            velocity_bc(u, v, w, params);
            predict_velocity_u(u, v, w, u_predict, params);
            predict_velocity_v(u, v, w, v_predict, params);
            predict_velocity_w(u, v, w, w_predict, params);
        )

        TIME("> Predicted divergence...", 
            communicate(u_predict, buffers, params);
            communicate(v_predict, buffers, params);
            communicate(w_predict, buffers, params);
            velocity_bc(u_predict, v_predict, w_predict, params);
            compute_divergence(u_predict, v_predict, w_predict, divergence, params);
        )

        TIME("> Solve pressure...", 
            communicate(pressure, buffers, params);
            pressure_bc(pressure, params);
            poisson_solve(pressure, divergence, &residual_poisson, &iters_poisson, poisson, buffers, params);
        )

        TIME("> Project velocities...", 
            communicate(pressure, buffers, params);
            pressure_bc(pressure, params);
            velocity_project(u_predict, v_predict, w_predict, u, v, w, pressure, params);
        )

        TIME("> Final divergence...", 
            communicate(u, buffers, params);
            communicate(v, buffers, params);
            communicate(w, buffers, params);
            velocity_bc(u, v, w, params);
            compute_divergence(u, v, w, divergence, params);
        )

        params->time += params->dt;

        if (params->time > params->time_output - 0.5*params->dt)
        {
            TIME_OUTPUT
            (
                output(u, "u", params);
                output(v, "v", params);
                output(w, "w", params);
                output(divergence, "divergence", params);
                output(pressure, "pressure", params);
            )

            params->n_output += 1;

            params->time_output += params->dt_output;
        }

        if (params->time > params->time_log - 0.5*params->dt)
        {
            if (params->process_rank == 0)
            {
                printf("--------------------------------------------------------------------------------\n");
                printf("Step completed!\n");
                printf("    Maximum divergence: = %.3e, CGM iters = %d\n", residual_poisson, iters_poisson);
            }
            params->time_log += params->dt_log;
        }
    }

    output(u, "u", params);
    output(v, "v", params);
    output(w, "w", params);
    output(pressure, "pressure", params);

    free_poisson(poisson);

    free_buffers(buffers);

    free(divergence);

    free(pressure);

    free(u_predict);
    free(v_predict);
    free(w_predict);

    free(u);
    free(v);
    free(w);

    free(params);
    MPI_Finalize();
    return 0;
}