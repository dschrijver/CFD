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
    double *r_poisson, *p_poisson, *q_poisson;
    allocate(&r_poisson, params);
    allocate(&p_poisson, params);
    allocate(&q_poisson, params);

    initialize_fields(u, v, w, pressure, params);

    params->time = 0.0;
    params->time_output = 0.0;
    params->time_log = 0.0;
    params->n_output = 0;

    double residual_poisson = 0.0;
    int iters_poisson = 0;
    while (params->time < params->t_final)
    {
        if (params->time > params->time_output - 0.5*params->dt)
        {
            output(u, "u", params);
            output(v, "v", params);
            output(w, "w", params);
            output(pressure, "pressure", params);

            params->n_output += 1;

            params->time_output += params->dt_output;
        }

        communicate(u, params);
        communicate(v, params);
        communicate(w, params);
        velocity_bc(u, v, w, params);
        predict_velocities(u, v, w, u_predict, v_predict, w_predict, params);

        communicate(u_predict, params);
        communicate(v_predict, params);
        communicate(w_predict, params);
        velocity_bc(u_predict, v_predict, w_predict, params);
        compute_divergence(u_predict, v_predict, w_predict, divergence, params);

        communicate(pressure, params);
        pressure_bc(pressure, params);
        poisson_solve(pressure, divergence, r_poisson, p_poisson, q_poisson, &residual_poisson, &iters_poisson, params);

        communicate(pressure, params);
        pressure_bc(pressure, params);
        velocity_project(u_predict, v_predict, w_predict, u, v, w, pressure, params);

        if (params->time > params->time_log - 0.5*params->dt)
        {
            if (params->process_rank == 0)
            {
                printf("--------------------------------------------------------------------------------\n");
                printf("Step completed!\n");
                printf("    time = %.3e\n", params->time);
                printf("    Poisson: residual = %.3e, iters = %d\n", residual_poisson, iters_poisson);
            }
            params->time_log += params->dt_log;
        }

        params->time += params->dt;
    }

    output(u, "u", params);
    output(v, "v", params);
    output(w, "w", params);
    output(pressure, "pressure", params);

    free(r_poisson);
    free(p_poisson);
    free(q_poisson);

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