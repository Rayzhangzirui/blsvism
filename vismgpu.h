#ifndef VISMGPU_H
#define VISMGPU_H

#include "oclstart.h"
#include "surf.h"
#include "Solute.h"

struct ParallelStruct
{
   int dim;
   cl_uint device;
   cl_context context;
   cl_command_queue command_queue;
   size_t *local_work_size;
   size_t *global_work_size;
};

//set worker size
void setparallel(ParallelStruct &parallel, int *nx, int dim);

struct clStruct
{
	cl_platform_id     platform;
	cl_device_id       device;
	ParallelStruct     parallel;
	cl_program         program;
	cl_kernel          getLJcl;//not used
	cl_kernel          getLJm2cl;
	cl_kernel          getphicl;//not used
	cl_kernel          getphi2cl;
	cl_kernel          getLJsolcl;
	cl_kernel          getGBcomponent;//not used
	cl_kernel          getCFAcomp;
	cl_kernel 			getoutsideGBkernel;
	clStruct();
};

// parallel calcualte solute-solute LJ 
double ocl_LJsolsol(SoluteStruct &sol, SoluteStruct &psol);

// parallel initialize surface of single solute
void ocl_initsurf(SurfaceStruct &surf, SoluteStruct &sol, GridData &grid, bool tightfit);

// helper for parallel LJ, read data from buffer
void readbufferdouble(double ***the3darray, cl_mem &thebuffer, int xdim, int ydim, int zdim, cl_command_queue &command_queue);

// helper for parallel GB, read data from buffer
void readbufferdouble4d(double ****the4darray, cl_mem &thebuffer, int xdim, int ydim, int zdim, cl_command_queue &command_queue);

// calculate LJ field of sol in parallel

void parLJ(SoluteStruct &sol, double*** LJ, GridData &grid) __attribute__ ((deprecated));

// calculate LJ field of sol in parallel
void parLJm2(SoluteStruct &sol,double*** LJ, GridData &grid);

// parallel initialize CFA field component
void parCFA(SoluteStruct &sol, double*** GB, GridData &grid);

// parallel initialize CFA field component
void parCFAcomp(SoluteStruct &sol, double**** CFAcomp, GridData &grid)__attribute__ ((deprecated));

// parallel initialize CFA field second method
void parCFAcomp2(SoluteStruct &sol, double**** CFAcomp, GridData &grid);

// outside box integration test
double paroutsidebox(int ntheta, int nphi, int nr, GridData &grid);

// outside box
double ocl_outsidebox_cfa(SoluteStruct &sol, int ntheta, int nphi, int nr, GridData &boxgrid);
double ocl_outsidebox_cfa(SoluteStruct &sol, GridData &grid);

// 1D LJ
void parLJ1d(SoluteStruct &sol, double* LJ1d, GridData &grid);


#endif