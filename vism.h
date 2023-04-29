#ifndef VISM_H
#define VISM_H
#define CL_SILENCE_DEPRECATION //opencl deprecated in macOS 10.14 Mojave

#include "util.h"
#include "heap.h"
#include "kernel.h"
#include "Solute.h"
#include "globals.h"
#include "surf.h"


struct FlipInfo{
	int numflip;
	double time;
};

FlipInfo binaryflowheaplistinterfaceonly(SurfaceStruct &surf, SoluteStruct &sol, double ***lj, double ***elec, KernelStruct &kernel, GridData &grid);

double energychangelist(int *index, SurfaceStruct &surf, SoluteStruct &sol, double ***lj, double ***elec, KernelStruct &kernel, GridData &grid);

double EnergyChangeField(int *index, SurfaceStruct &surf, double ***field, GridData &grid);

double EnergyChangeSurf(int *index, SurfaceStruct &surf, SoluteStruct &sol, KernelStruct &kernel, GridData &grid);

// compute total energy of the system, after caclulating LJ field and flow
Energy VismEnergyInit(SoluteStruct &sol, double*** lj, double ***cfa, SurfaceStruct& surf, KernelStruct& kernel, GridData& grid, bool tightfit);

// compute total energy of the system, after caclulating LJ field and flow
std::tuple<Energy, FlipInfo> VismEnergyInitFlow(SoluteStruct &sol, double*** lj, double ***cfa, SurfaceStruct& surf, KernelStruct& kernel, GridData& grid, bool tightfit = true);

// compute surface energy
double GSurf(SurfaceStruct &surf, KernelStruct &kernel, GridData &grid);

Energy VismEnergyOnly(SoluteStruct &sol, double*** lj, double ***cfa, SurfaceStruct& surf, KernelStruct& kernel, GridData& grid);

char nexttointerface(double ***phi, int *tindex, GridData &grid);
#endif