#include "fftconv.h"
#include "util.h"
#include <spdlog/spdlog.h>

bool globpow2 = false;
unsigned globfftwplan = FFTW_EXHAUSTIVE;

int sub2ind(int *index, int *nn, int ndim)
{
   int i, r, factor;

   r = index[ndim-1];
   factor = nn[ndim-1]+1;
   for (i = ndim-2; i >= 0; i--)
   {
      r += index[i]*factor;
      factor *= nn[i]+1;
   }

   return r;
}

double getlengthfftw(double ***phi, FFTWStruct &fftw, GridData &grid)
{
   int i, r, tindex[grid.dim];
   double value, Kval;
   double real, imag;
   double dV = 1.0;
   for (r = 0; r < grid.dim; r++)
      dV *= grid.dx[r];

   for (r = 0; r < fftw.N; r++)
   {
//      fftw.datain[r][0] = 1.0;
//      fftw.datain[r][1] = 0.0;
      fftw.data[r][0] = 1.0;
      fftw.data[r][1] = 0.0;
   }

   // set fft.data to be phi
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      if (evalarray(phi,tindex) <= 0.0)
      {
         r = sub2ind(tindex,fftw.nx,fftw.dim);  
//         fftw.datain[r][0] = 0.0;
         fftw.data[r][0] = 0.0;
      }

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   // fft of phi
   fftw_execute(fftw.planfor);

   // product of convolution
   for (r = 0; r < fftw.N; r++)
   {
//      fftw.datain[r][0] = fftw.dataout[r][0]*fftw.kernelout[r][0]-
//                          fftw.dataout[r][1]*fftw.kernelout[r][1];
//      fftw.datain[r][1] = fftw.dataout[r][0]*fftw.kernelout[r][1]+
//                          fftw.dataout[r][1]*fftw.kernelout[r][0];
      real = fftw.data[r][0];
      imag = fftw.data[r][1];
      fftw.data[r][0] = real*fftw.kernel[r][0]-imag*fftw.kernel[r][1];
      fftw.data[r][1] = real*fftw.kernel[r][1]+imag*fftw.kernel[r][0];
   }

   // ifft of product
   fftw_execute(fftw.planback);

   //  compute area by integrating inside
   value = 0.0;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      if (evalarray(phi,tindex) <= 0.0)
      {
         r = sub2ind(tindex,fftw.nx,fftw.dim);  
//         Kval = fftw.dataout[r][0]/fftw.N;
         Kval = fftw.data[r][0]/fftw.N;
         if (Kval != 0.0)
            value += Kval*dV*dV;
      }

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   return value;
}

// fftw.nx is upper bound power of two of kernel size

void getinitfftw(FFTWStruct &fftw, KernelStruct &kernel, GridData &grid)
{
   int i, tindex[grid.dim], sindex[grid.dim];

   if (globpow2 == true){
      for (i = 0; i < grid.dim; i++)
         fftw.nx[i] = static_cast<int>(round(exp(ceil(log2(grid.nx[i]+kernel.rad+1))*log(2))))-1;
   } else {
      for (i = 0; i < grid.dim; i++)
         fftw.nx[i] =grid.nx[i]+static_cast<int>(kernel.rad)-1;
   }
   
   fftw.N = 1;
   for (i = 0; i < grid.dim; i++)
      fftw.N *= fftw.nx[i]+1;
//   fftw.datain = fftw_alloc_complex(fftw.N);
//   fftw.dataout = fftw_alloc_complex(fftw.N);
//   fftw.kernelout = fftw_alloc_complex(fftw.N);
   fftw.data = fftw_alloc_complex(fftw.N);
   fftw.kernel = fftw_alloc_complex(fftw.N);

   fftw_import_wisdom_from_filename(PROJECTDIR "/binaryflowfftwisdom.dat");
   cout << "starting forward plan" << endl;
//   fftw.planfor = fftw_plan_dft_3d(fftw.nx[0]+1,fftw.nx[1]+1,fftw.nx[2]+1,fftw.datain,
//                                   fftw.dataout,FFTW_FORWARD,globfftwplan);
   fftw.planfor = fftw_plan_dft_3d(fftw.nx[0]+1,fftw.nx[1]+1,fftw.nx[2]+1,fftw.data,
                                   fftw.data,FFTW_FORWARD,globfftwplan);
   cout << "starting backward plan" << endl;
//   fftw.planback = fftw_plan_dft_3d(fftw.nx[0]+1,fftw.nx[1]+1,fftw.nx[2]+1,fftw.datain,
//                                    fftw.dataout,FFTW_BACKWARD,globfftwplan);
   fftw.planback = fftw_plan_dft_3d(fftw.nx[0]+1,fftw.nx[1]+1,fftw.nx[2]+1,fftw.data,
                                    fftw.data,FFTW_BACKWARD,globfftwplan);
   fftw_export_wisdom_to_filename(PROJECTDIR "/binaryflowfftwisdom.dat");

   cout << "FFT info " << fftw.nx[0] << " " << fftw.nx[1] << " " << fftw.nx[2] 
        << " and " << fftw.N << endl;

   // initialize as 0
   for (i = 0; i < fftw.N; i++)
   {
//      fftw.datain[i][0] = 0.0;
//      fftw.datain[i][1] = 0.0;
      fftw.data[i][0] = 0.0;
      fftw.data[i][1] = 0.0;
   }


   // copy kernel to 1d fftw.data
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= 2*kernel.rad)
   {
      for (i = 0; i < grid.dim; i++)
      {
         sindex[i] = tindex[i]-kernel.rad;
         if (sindex[i] < 0)
            sindex[i] += fftw.nx[i]+1;
      }
//      fftw.datain[sub2ind(sindex,fftw.nx,fftw.dim)][0] = evalarray(kernel.K,tindex);
      fftw.data[sub2ind(sindex,fftw.nx,fftw.dim)][0] = evalarray(kernel.K,tindex);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > 2*kernel.rad; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

/*
   int r = 0;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= fftw.nx[0])
   {
      chi2d << fftw.datain[r][0] << " " << endl;
      if (r != sub2ind(tindex,fftw.nx,fftw.dim))
      {
         cout << "bad at " << r << " " << sub2ind(tindex,fftw.nx,fftw.dim) << endl;
         exit(1);
      }
      r++;

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > fftw.nx[i]; i--)
      {
         if (i == grid.dim-1)
            chi2d << endl;
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
*/

   cout << "execute for kernel" << endl;
   fftw_execute(fftw.planfor);
   cout << "done execute for kernel" << endl;

   for (i = 0; i < fftw.N; i++)
   {
//      fftw.kernelout[i][0] = fftw.dataout[i][0];
//      fftw.kernelout[i][1] = fftw.dataout[i][1];
      fftw.kernel[i][0] = fftw.data[i][0];
      fftw.kernel[i][1] = fftw.data[i][1];
   }
}


void binaryflowfftw(SurfaceStruct &surf, SoluteStruct &sol, KernelStruct &kernel, FFTWStruct &fftw, GridData &grid)
{
   int i, step, numchange, minchange = 1, maxstep = 500;
   int tindex[grid.dim], sindex[grid.dim]; 
   double echange, phivalue, theenergy, real1, imag1, real2, imag2;

   int r, rindex[grid.dim];
   double Kval, maxerr = 0.0;
   double ***temp = matrix<double>(grid.nx[0],grid.nx[1],grid.nx[2]);
   clock_t time1, time2, steptime1, steptime2, fntime1, fntime2;
   
   char ***tube = matrix<char>(grid.nx[0],grid.nx[1],grid.nx[2]);
   int count; // count number of flipping in each step
   int totalcount = 0;

   fntime1 = clock();
   // time1 = clock();
   // binaryenergyfftw(theenergy,surf,sol,fftw,grid);
   // time2 = clock();
   // cout << "binaryenergyfftw took " << static_cast<double>(time2-time1)/CLOCKS_PER_SEC
   //      << " SECONDS" << endl;

   // mark near interface
   for (numchange = minchange,step = 1; numchange >= minchange && step <= maxstep; 
        step++)
   {
      cout << "step = " << step << endl;

      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         if (nexttointerface(surf.phi,tindex,grid))
            setvalarray<char>(tube,tindex,(char)1);
         else
            setvalarray<char>(tube,tindex,(char)0);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }

      steptime1 = clock();
      numchange = 0;

      time1 = clock();
      

      // set phi to ffw.data, inside 
      r = 0;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= fftw.nx[0])
      {
         for (i = 0; i < grid.dim && tindex[i] <= grid.nx[i]; i++);
         if (i >= grid.dim) 
         {
            // tindex inside original grid
            phivalue = evalarray(surf.phi,tindex);
            fftw.data[r][0] = phivalue;
         }
         else
         {
            // tindex outside original grid due to padding, set 1
            fftw.data[r][0] = 1.0;
         }
         fftw.data[r][1] = 0.0; //complex part always 0
         r++;
   
         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > fftw.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      time2 = clock();
      cout << "   init fftw took " << static_cast<double>(time2-time1)/CLOCKS_PER_SEC
           << " SECONDS" << endl;
      time1 = time2;

      fftw_execute(fftw.planfor);
      time2 = clock();
      cout << "   first fftw took " << static_cast<double>(time2-time1)/CLOCKS_PER_SEC
           << " SECONDS" << endl;
      time1 = time2;

      for (r = 0; r < fftw.N; r++)
      {
/*
         real1 = fftw.dataout[r][0];
         imag1 = fftw.dataout[r][1];
         real2 = fftw.kernelout[r][0];
         imag2 = fftw.kernelout[r][1];
         fftw.datain[r][0] = real1*real2-imag1*imag2;
         fftw.datain[r][1] = real1*imag2+imag1*real2;
*/
         real1 = fftw.data[r][0];
         imag1 = fftw.data[r][1];
         fftw.data[r][0] = real1*fftw.kernel[r][0]-imag1*fftw.kernel[r][1];
         fftw.data[r][1] = real1*fftw.kernel[r][1]+imag1*fftw.kernel[r][0];
      }
      time2 = clock();
      cout << "   convolution took " << static_cast<double>(time2-time1)/CLOCKS_PER_SEC
           << " SECONDS" << endl;
      time1 = time2;
     
      fftw_execute(fftw.planback);
      for (r = 0; r < fftw.N; r++)
      {
//         fftw.dataout[r][0] /= fftw.N;
         fftw.data[r][0] /= fftw.N;
      }
      time2 = clock();
      cout << "   second fftw took " << static_cast<double>(time2-time1)/CLOCKS_PER_SEC
           << " SECONDS" << endl;
      time1 = time2;

      count = 0;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
//         echange = surf.gamma0*fftw.dataout[sub2ind(tindex,fftw.nx,fftw.dim)][0]*
//                   grid.dV-sol.rho0*evalarray(sol.LJ,tindex);
         echange = ::gamma0*fftw.data[sub2ind(tindex,fftw.nx,fftw.dim)][0]*grid.dV - sol.rho0*evalarray(sol.LJ,tindex) - evalarray(sol.GB,tindex);
         phivalue = evalarray(surf.phi,tindex);
         if (phivalue < 0.0)
            echange = -echange;
         if (echange < 0.0)
         {
            setvalarray(surf.phi,tindex,-phivalue);
            numchange++;
            if (evalarray(tube,tindex) == (char)0)
            {
//               cout << "changing point away from interface " << tindex[0] << " " 
//                    << tindex[1] << " " << tindex[2] << endl;
               count++;
            }
         }

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      time2 = clock();
      cout << "   flipping took " << static_cast<double>(time2-time1)/CLOCKS_PER_SEC
           << " SECONDS" << endl;

      cout << "   changed " << numchange << " and " << count << " away from interface" 
           << endl;

      steptime2 = clock();
      cout << "   whole step took " 
           << static_cast<double>(steptime2-steptime1)/CLOCKS_PER_SEC 
           << " SECONDS" << endl;
      steptime1 = steptime2;

      totalcount += numchange;
      // Energy e = VismEnergyOnly(sol, sol.LJ, sol.GB, surf, kernel, grid);
      // spdlog::info( "step {}, fft energy surf={}, ljswbox={}, total={}\n", step,  e.surf, e.ljswbox,  e.total());
   }

   // time1 = clock();
   // binaryenergyfftw(theenergy,surf,sol,fftw,grid);
   // time2 = clock();
   // cout << "binaryenergyfftw took " << static_cast<double>(time2-time1)/CLOCKS_PER_SEC
   //      << " SECONDS" << endl;

   fntime2 = clock();
   spdlog::info("function took {} SECONDS\n",static_cast<double>(fntime2-fntime1)/CLOCKS_PER_SEC);
   spdlog::info("fft total flip {} \n", totalcount);
}
