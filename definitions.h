#ifndef DEFINITIONS_H
#define DEFINITIONS_H

// ---------------------------
//     Boundary conditions        
// ---------------------------
// #define XPERIODIC
// #define YPERIODIC
// #define ZPERIODIC

// #define LEFT_NOSLIP
// #define RIGHT_NOSLIP
// #define BOTTOM_NOSLIP
// #define TOP_NOSLIP
// #define BACK_NOSLIP
// #define FRONT_NOSLIP

#define LEFT_PRESSURE 0.0
#define RIGHT_PRESSURE 0.0
#define BOTTOM_PRESSURE 0.0
#define TOP_PRESSURE 0.0
#define BACK_PRESSURE 0.0
#define FRONT_PRESSURE 0.0

// #define INI_RHS_SINE // Standard option, to test periodic boundaries
// #define INI_RHS_COSINE // To test no-slip boundaries
#define INI_RHS_ODDSINE // To test Dirichlet boundaries
// #define INI_RHS_NOISE // To test solver robustness

#endif