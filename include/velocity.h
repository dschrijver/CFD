#ifndef VELOCITY_H
#define VELOCITY_H

#include "datatypes.h"

void predict_velocity_u(double *u, double *v, double *w, double *u_predict, ParamBag *params);
void predict_velocity_v(double *u, double *v, double *w, double *v_predict, ParamBag *params);
void predict_velocity_w(double *u, double *v, double *w, double *w_predict, ParamBag *params);
void compute_divergence(double *u, double *v, double *w, double *divergence, ParamBag *params);
void velocity_project(double *u_predict, double *v_predict, double *w_predict, double *u, double *v, double *w, double *pressure, ParamBag *params);
double compute_CFL(double *u, double *v, double *w, ParamBag *params);

#endif