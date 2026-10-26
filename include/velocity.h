#ifndef VELOCITY_H
#define VELOCITY_H

#include "datatypes.h"

void predict_velocities(double *u, double *v, double *w, double *u_predict, double *v_predict, double *w_predict, ParamBag *params);
void compute_divergence(double *u, double *v, double *w, double *divergence, ParamBag *params);
void velocity_project(double *u_predict, double *v_predict, double *w_predict, double *u, double *v, double *w, double *pressure, ParamBag *params);

#endif