#include "../include/datatypes.h"
#include "../include/utils.h"

void predict_velocity_u(double *u, double *v, double *w, double *u_predict, ParamBag *params)
{
    UNPACK_GRID

    const double dt = params->dt;
    const double Re = params->Re;

    FOR_DOMAIN
    {
        const int idx0 = IDX(i, j, k);

        // Convective term
        const double u_W = 0.5 * (u[IDX(i - 1, j, k)] + u[idx0]);
        const double u_E = 0.5 * (u[idx0] + u[IDX(i + 1, j, k)]);
        const double u_S = 0.5 * (u[IDX(i, j - 1, k)] + u[idx0]); 
        const double u_N = 0.5 * (u[idx0] + u[IDX(i, j + 1, k)]); 
        const double u_B = 0.5 * (u[IDX(i, j, k - 1)] + u[idx0]); 
        const double u_F = 0.5 * (u[idx0] + u[IDX(i, j, k + 1)]); 

        const double v_S = 0.5 * (v[IDX(i - 1, j, k)] + v[idx0]); 
        const double v_N = 0.5 * (v[IDX(i - 1, j + 1, k)] + v[IDX(i, j + 1, k)]); 

        const double w_B = 0.5 * (w[IDX(i - 1, j, k)] + w[idx0]); 
        const double w_F = 0.5 * (w[IDX(i - 1, j, k + 1)] + w[IDX(i, j, k + 1)]); 

        const double dx_uu = (u_E*u_E - u_W*u_W) / dx;
        const double dy_uv = (u_N*v_N - u_S*v_S) / dx;
        const double dz_uw = (u_F*w_F - u_B*w_B) / dx;

        // Viscous term
        const double dxx_u = (u[IDX(i - 1, j, k)] - 2.0*u[idx0] + u[IDX(i + 1, j, k)]) / (dx*dx);
        const double dyy_u = (u[IDX(i, j - 1, k)] - 2.0*u[idx0] + u[IDX(i, j + 1, k)]) / (dx*dx);
        const double dzz_u = (u[IDX(i, j, k - 1)] - 2.0*u[idx0] + u[IDX(i, j, k + 1)]) / (dx*dx);

        u_predict[idx0] = u[idx0] + dt * (1.0/Re*(dxx_u + dyy_u + dzz_u) - (dx_uu + dy_uv + dz_uw) + params->Fx);
    }
}

void predict_velocity_v(double *u, double *v, double *w, double *v_predict, ParamBag *params)
{
    UNPACK_GRID

    const double dt = params->dt;
    const double Re = params->Re;

    FOR_DOMAIN
    {
        const int idx0 = IDX(i, j, k);

        // Convective term
        const double v_W = 0.5 * (v[IDX(i - 1, j, k)] + v[idx0]); 
        const double v_E = 0.5 * (v[idx0] + v[IDX(i + 1, j, k)]); 
        const double v_S = 0.5 * (v[IDX(i, j - 1, k)] + v[idx0]);
        const double v_N = 0.5 * (v[idx0] + v[IDX(i, j + 1, k)]);
        const double v_B = 0.5 * (v[IDX(i, j, k - 1)] + v[idx0]); 
        const double v_F = 0.5 * (v[idx0] + v[IDX(i, j, k + 1)]); 

        const double u_W = 0.5 * (u[IDX(i, j - 1, k)] + u[idx0]); 
        const double u_E = 0.5 * (u[IDX(i + 1, j - 1, k)] + u[IDX(i + 1, j, k)]); 

        const double w_B = 0.5 * (w[IDX(i, j - 1, k)] + w[idx0]); 
        const double w_F = 0.5 * (w[IDX(i, j - 1, k + 1)] + w[IDX(i, j, k + 1)]); 

        const double dx_vu = (v_E*u_E - v_W*u_W) / dx;
        const double dy_vv = (v_N*v_N - v_S*v_S) / dx;
        const double dz_vw = (v_F*w_F - v_B*w_B) / dx;

        // Viscous term
        const double dxx_v = (v[IDX(i - 1, j, k)] - 2.0*v[idx0] + v[IDX(i + 1, j, k)]) / (dx*dx);
        const double dyy_v = (v[IDX(i, j - 1, k)] - 2.0*v[idx0] + v[IDX(i, j + 1, k)]) / (dx*dx);
        const double dzz_v = (v[IDX(i, j, k - 1)] - 2.0*v[idx0] + v[IDX(i, j, k + 1)]) / (dx*dx);

        v_predict[idx0] = v[idx0] + dt * (1.0/Re*(dxx_v + dyy_v + dzz_v) - (dx_vu + dy_vv + dz_vw) + params->Fy);
    }
}

void predict_velocity_w(double *u, double *v, double *w, double *w_predict, ParamBag *params)
{
    UNPACK_GRID

    const double dt = params->dt;
    const double Re = params->Re;

    FOR_DOMAIN
    {
        const int idx0 = IDX(i, j, k);

        // Convective term
        const double w_W = 0.5 * (w[IDX(i - 1, j, k)] + w[idx0]); 
        const double w_E = 0.5 * (w[idx0] + w[IDX(i + 1, j, k)]); 
        const double w_S = 0.5 * (w[IDX(i, j - 1, k)] + w[idx0]); 
        const double w_N = 0.5 * (w[idx0] + w[IDX(i, j + 1, k)]);
        const double w_B = 0.5 * (w[IDX(i, j, k - 1)] + w[idx0]);
        const double w_F = 0.5 * (w[idx0] + w[IDX(i, j, k + 1)]);

        const double u_W = 0.5 * (u[IDX(i, j, k - 1)] + u[idx0]); 
        const double u_E = 0.5 * (u[IDX(i + 1, j, k - 1)] + u[IDX(i + 1, j, k)]); 

        const double v_S = 0.5 * (v[IDX(i, j, k - 1)] + v[idx0]); 
        const double v_N = 0.5 * (v[IDX(i, j + 1, k - 1)] + v[IDX(i, j + 1, k)]); 

        const double dx_wu = (w_E*u_E - w_W*u_W) / dx;
        const double dy_wv = (w_N*v_N - w_S*v_S) / dx;
        const double dz_ww = (w_F*w_F - w_B*w_B) / dx;

        // Viscous term
        const double dxx_w = (w[IDX(i - 1, j, k)] - 2.0*w[idx0] + w[IDX(i + 1, j, k)]) / (dx*dx);
        const double dyy_w = (w[IDX(i, j - 1, k)] - 2.0*w[idx0] + w[IDX(i, j + 1, k)]) / (dx*dx);
        const double dzz_w = (w[IDX(i, j, k - 1)] - 2.0*w[idx0] + w[IDX(i, j, k + 1)]) / (dx*dx);

        w_predict[idx0] = w[idx0] + dt * (1.0/Re*(dxx_w + dyy_w + dzz_w) - (dx_wu + dy_wv + dz_ww) + params->Fz);
    }
}

void compute_divergence(double *u, double *v, double *w, double *divergence, ParamBag *params)
{
    UNPACK_GRID

    const double dt = params->dt;

    FOR_DOMAIN
    {
        const int idx0 = IDX(i, j, k);

        const double dx_u = (u[IDX(i + 1, j, k)] - u[idx0]) / dx;
        const double dy_v = (v[IDX(i, j + 1, k)] - v[idx0]) / dx;
        const double dz_w = (w[IDX(i, j, k + 1)] - w[idx0]) / dx;

        divergence[idx0] = (dx_u + dy_v + dz_w) / dt;
    }
}

void velocity_project(double *u_predict, double *v_predict, double *w_predict, double *u, double *v, double *w, double *pressure, ParamBag *params)
{
    UNPACK_GRID

    const double dt = params->dt;

    FOR_DOMAIN
    {
        const int idx0 = IDX(i, j, k);

        const double dx_p = (pressure[idx0] - pressure[IDX(i - 1, j, k)]) / dx;
        const double dy_p = (pressure[idx0] - pressure[IDX(i, j - 1, k)]) / dx;
        const double dz_p = (pressure[idx0] - pressure[IDX(i, j, k - 1)]) / dx;

        u[idx0] = u_predict[idx0] - dt*dx_p;
        v[idx0] = v_predict[idx0] - dt*dy_p;
        w[idx0] = w_predict[idx0] - dt*dz_p;
    }
}