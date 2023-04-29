#ifndef KERNEL_H
#define KERNEL_H

#include "util.h"

struct KernelStruct
{
   double ***K; //3d array of kernel
   int rad;
   int start;

   int **index;// index of nonzero element in kernel
   int N; //N+1 is the number of non zero kernel
   KernelStruct();
   KernelStruct(GridData &grid, int radius);//constructor
   KernelStruct(GridData &grid, double radius);//constructor
   KernelStruct(const KernelStruct &k);//copy constructor
   ~KernelStruct();//destructor
};


double getkernel(int *index, int rad, int dim);
double getkernel(double r);

// integrate a_0(theta) from 0 to 1 
double getdenom(int N);

// simpson's rule to get denominator
double getdenom(double rad, int N);

// initialize kernel
void getinitkernel2(KernelStruct &kernel, GridData &grid, int rad);

void getinitkerneldouble(KernelStruct &kernel, GridData &grid, double radius);

double getlengthlist(double ***phi, KernelStruct &kernel, GridData &grid);

#endif