// kernel using vectors of value and subscripts, for more flexible experiment of area convergence
// not used in VISM study
// only used in kerneltest.cpp
#ifndef KERNEL2_H
#define KERNEL2_H

#include "util.h"
#include "Vector.h"
#include <functional>

extern const double A0_SIN2;
extern const function<double(double)> KERNEL_SIN2;

extern const double A0_IND;
extern const function<double(double)> KERNEL_IND;

extern const double A0_COS;
extern const function<double(double)> KERNEL_COS;



struct Kernel {
	GridData grid;
	double radius;
	int rad; //int half rad, kernel size 2rad + 1 
	int N;
	double ***K; //3d array of kernel
	std::vector<Sub> index; //indices in local coordinate
	std::vector<double> value; // value at indices

	Kernel(GridData &grid, int multiple_dx, int code = 1); // constructor
	Kernel(GridData &grid, double rad, int code = 1); // constructor
	Kernel(){};
	~Kernel();           //  and destructor.
	void set_zero();
};


void InitKernelDoubleRadius(Kernel &kernel, GridData &grid, double rad, std::function<double(double)> kernel_fcn, double a0);

void InitKernelScale(Kernel &kernel, GridData &grid, double scale, function<double(double)> kernel_fcn, double a0);

void InitKernelCode(Kernel &kernel, GridData &grid, double rad, int code);

// calculate are using phi
double ComputeAreaGrid(double ***phi, Kernel &kernel, GridData &grid);

// calculate area using bsurface
double ComputeAreaFcn(const Kernel &kernel, const GridData &grid, const std::function<bool(va::Vector)>& bsurface);

double ComputeSurfaceIntegral(const Kernel &kernel, const GridData &grid, const function<bool(va::Vector)>& bsurface, const function<double(va::Vector)>& fcn);

// compute integrate a0(theta) from 0 to 1 using midpoint rule with n number of grid 
double ComputeA0Midpoint(int n, const std::function<double(double)>& kernel_fcn);

// compute denominator of kernel using midpoint rule
double ComputeDenomMidpoint(int n, double kernel_real_radius, const GridData& grid, std::function<double(double)> kernel_fcn);

// linear convolution with zero padding
double ConvolveLinear(double ***phi, Kernel &kernel);


#endif