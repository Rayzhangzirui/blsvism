#ifndef VISM_ADDON_H
#define VISM_ADDON_H

#include "cfa.h"
#include "surf.h"
#include <sstream>
#include <array>

#ifdef OPENCL
#include "vismgpu.h"
extern clStruct OCL;
#endif

// un-average with water, use for epsilon and sigma that is alread averaged
inline double eps_unave(double eps){
   return eps*eps/WATEREPS;
}

inline double sigma_unave(double sigma){
   return 2*sigma-WATERSIG;
}


//sequential initialize surface
void seq_getinitsurf(SurfaceStruct &surf, SoluteStruct &sol, GridData &grid, bool tightfit);

//sequential initialize LJ
void seqLJ(SoluteStruct &sol, double ***LJ, GridData &grid);


// compute total LJ solute water in box
double OutsideIntegral(SurfaceStruct &surf, double ***f, GridData &grid, double scale);


void getinitsurf(SurfaceStruct &surf, SoluteStruct &sol, GridData &grid, bool tightfit);

void cal_LJ(SoluteStruct &sol, double ***LJ, GridData &grid);

double LJoutsidebox3D(SoluteStruct& sol, const GridData &grid);

double GLjSolWatAtom(SoluteStruct &sol, int ind, SurfaceStruct &surf, double cutoff, GridData &grid);


#endif //header guard
