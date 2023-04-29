#include <fftw3.h>
#include "util.h"
#include "surf.h"
#include "kernel2.h"
#include "Solute.h"
#include "vism.h"

extern bool globpow2;

struct FFTWStruct
{
   static const int dim = 3;
//   fftw_complex *datain;
//   fftw_complex *dataout, *kernelout;
   fftw_complex *data;
   fftw_complex *kernel;
   fftw_plan planfor;
   fftw_plan planback;
   int N;
   int nx[dim];
};

double getlengthfftw(double ***phi, FFTWStruct &fftw, GridData &grid);
void getinitfftw(FFTWStruct &fftw, KernelStruct &kernel, GridData &grid);
void binaryflowfftw(SurfaceStruct &surf, SoluteStruct &sol, KernelStruct &kernel, FFTWStruct &fftw, GridData &grid);
