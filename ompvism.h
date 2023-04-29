#ifndef OMPVISM_H
#define OMPVISM_H

#include "cfangpu.h"

void ompCFA(SoluteStruct &sol, double*** GB, GridData &grid);

void ompCFAcomp(SoluteStruct &sol, double**** CFAcomp,  GridData &grid);

void ompLJ(SoluteStruct &sol, double*** LJ, GridData &grid);

void ompLJ1d(SoluteStruct &sol, double* LJ1d, GridData &grid);

void omp_initsurf(SurfaceStruct &surf, SoluteStruct &sol, GridData &grid, bool tightfit);

double omp_LJsolsol(SoluteStruct &Asol,SoluteStruct &Bsol);

double omp_ELECsolsol(SoluteStruct &Asol,SoluteStruct &Bsol,double permittivity);

double omp_outsidebox_cfa(SoluteStruct &sol, GridData &boxgrid, int ntheta = 50, int nphi = 50, int nr = 50);

#endif //header guard