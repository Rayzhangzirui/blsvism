//this file supplement cfangpu.cpp
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <queue> //in bfs
#include <array>
#include <algorithm>
//--

#include "vism_addon.h"
#include "globals.h"
#include "cfa.h"

using namespace std;
//--

#ifdef OPENMP
#include <omp.h>
#endif

#ifdef OPENCL
#include "vismgpu.h"
#endif


double LJwater(int *index, SoluteStruct &sol, GridData &grid)
{
	double val = 0.0;
	double epsilon, sigma, x[grid.dim];
	

	sub2coord(x,index,grid);

	for (int r = 0; r < sol.dnum; r++)
	{
		sigma = (sol.sigma[r] + WATERSIG)/2;
		epsilon = sqrt(sol.epsilon[r] * WATEREPS);
		
		double temp = 0.0;
		for (int s = 0; s < grid.dim; s++)
			temp += (x[s]-sol.d[r][s])*(x[s]-sol.d[r][s]);
		if (temp <= TOL)
			temp = TOL;
		// val += 4.0*epsilon*
		//        (exp(12.0*log(sigma/temp))-exp(6.0*log(sigma/temp)));
		val += 4.0*epsilon*(pow(sigma*sigma/temp,6.0)-pow(sigma*sigma/temp,3.0));//pow is faster than exp log, 60x60x60 grid, p53/mdm2, time is 2s vs. 13s.
	}

	return val;
}


void seq_getinitsurf(SurfaceStruct &surf, SoluteStruct &sol, GridData &grid,bool tightfit){
	int tindex[grid.dim];
	double x[grid.dim];

	for (int i=0;i<=grid.nx[0];i++){
		for (int j=0;j<=grid.nx[1];j++){
			for (int k=0;k<=grid.nx[2];k++){
				tindex[0]=i;
				tindex[1]=j;
				tindex[2]=k;			
				sub2coord(x,tindex,grid);
				// setvalarray(surf.phi,tindex,1.0);
				if (tightfit){
					//if tightfit, set -1 if close to any atom
					for(int n=0; n<sol.dnum;n++){
						if (sqrt((x[0]-sol.d[n][0])*(x[0]-sol.d[n][0])+
								(x[1]-sol.d[n][1])*(x[1]-sol.d[n][1])+
								(x[2]-sol.d[n][2])*(x[2]-sol.d[n][2])) <= (sol.sigma[n]+WATERSIG)/2) {
							setvalarray(surf.phi,tindex,-1.0);
						}
					}
				}else{
					//if loosefit, set everything as inside
					setvalarray(surf.phi,tindex,-1.0);
				}
			}
		}
	}

}

//sequential initialize surface
void getinitsurf(SurfaceStruct &surf, SoluteStruct &sol, GridData &grid, bool tightfit){
	
	#ifdef OPENCL
		ocl_initsurf(surf, sol, grid, tightfit);
	#else
		seq_getinitsurf(surf, sol, grid, tightfit);
	#endif
	
	
}


void seqLJ(SoluteStruct &sol, double*** LJ, GridData &grid){
	#pragma omp parallel for collapse(3)
		for(int i=0; i<=grid.nx[0]; i++){
			for(int j=0; j<=grid.nx[1]; j++){
				for(int k=0; k<=grid.nx[2]; k++){
					int tindex[3] = {i,j,k};
					setvalarray(LJ,tindex,LJwater(tindex,sol,grid));
				}
			}
		}
}

//
void cal_LJ(SoluteStruct &sol, double*** LJ, GridData &grid){
	#ifdef OPENCL
	parLJm2(sol, LJ, grid);
	#else
	seqLJ(sol, LJ,grid);
	#endif
}



// compute total LJ solute water in box
double OutsideIntegral(SurfaceStruct &surf, double ***f, GridData &grid, double scale){
    double total = 0.0;
    #pragma omp parallel for collapse(3)
    for (int i = 0;i <= grid.nx[0]; i++){
        for (int j = 0;j <= grid.nx[1]; j++){
            for (int k = 0;k <= grid.nx[2]; k++){
                if (surf.phi[i][j][k] > 0.0){
                    #pragma omp atomic update
                    total += f[i][j][k]*grid.dV;
                }
            }
        }
    }

    return total*scale;
}

// given surface, and cut off, compute contribution LJ of a single atom
double GLjSolWatAtom(SoluteStruct &sol, int ind, SurfaceStruct &surf, double cutoff, GridData &grid){
    double LJpart = 0.0;
    double sigma = (sol.sigma[ind] + WATERSIG)/2;
	double epsilon = sqrt(sol.epsilon[ind] * WATEREPS);

    Coord r = {sol.d[ind][0], sol.d[ind][1], sol.d[ind][2]};
    #pragma omp parallel for collapse(3)
    for (int i = 0;i <= grid.nx[0]; i++){
        for (int j = 0;j <= grid.nx[1]; j++){
            for (int k = 0;k <= grid.nx[2]; k++){

            	Coord x = sub2coord(Sub {{i,j,k}}, grid);
            	double d = dist2(x,r);
            	if(d < cutoff & surf.phi[i][j][k] > 0.0){
            		#pragma omp atomic update
            		LJpart += 4.0*epsilon*(pow(sigma/d,12.0)-pow(sigma/d,6.0));
            	}
            }
        }
    }
    
    return LJpart;
}





double outsidebox1D(int r, const SoluteStruct& sol, GridData grid)
{
// assumes gridbox is cube
	int ntheta = 100*grid.nx[0];
	double SS = (grid.b[0]-grid.a[0])/2.0+0.5*grid.dx[0];
//   double SS = (grid.b[0]-grid.a[0])/2.0;
	double theta, dtheta;
	double sigma, epsilon;
	double x0, y0, z0, theta1, theta2, C, value, RR;
	double x, y, z, xx, yy, zz, temp, temp1, value2 = 0.0;
	int i, j, k;

	// epsilon = sol.epsilon[r];
	// sigma = sol.sigma[r];
	// average with water

	sigma = (sol.sigma[r] + WATERSIG)/2;
	epsilon = sqrt(sol.epsilon[r] * WATEREPS);

	x0 = sol.d[r][0]-(grid.a[0]+grid.b[0])/2.0;
	y0 = sol.d[r][1]-(grid.a[1]+grid.b[1])/2.0;
	z0 = sol.d[r][2]-(grid.a[2]+grid.b[2])/2.0;

	theta1 = atan2(-SS-y0,SS-x0);
	theta2 = atan2(SS-y0,SS-x0);
	C = fabs(SS-x0);
	value = exp(12.0*log(sigma))/(9.0*exp(9.0*log(C)))*63.0/256.0*M_PI*
			(sin(theta2)*(exp(8.0*log(cos(theta2)))/9.0+
						 8.0/9.0*(exp(6.0*log(cos(theta2)))/7.0+
									6.0/7.0*(exp(4.0*log(cos(theta2)))/5.0+
											4.0/5.0*(cos(theta2)*cos(theta2)/3.0
													+2.0/3.0))))-
			sin(theta1)*(exp(8.0*log(cos(theta1)))/9.0+
						 8.0/9.0*(exp(6.0*log(cos(theta1)))/7.0+
									6.0/7.0*(exp(4.0*log(cos(theta1)))/5.0+
											4.0/5.0*(cos(theta1)*cos(theta1)/3.0
													+2.0/3.0)))))-
			exp(6.0*log(sigma))/(3.0*exp(3.0*log(C)))*3.0/8.0*M_PI*
			(sin(theta2)*(cos(theta2)*cos(theta2)/3.0+2.0/3.0)-
			sin(theta1)*(cos(theta1)*cos(theta1)/3.0+2.0/3.0));
	dtheta = (theta2-theta1)/ntheta;
	for (i = 0; i <= ntheta; i++)
	{
		theta = theta1+i*dtheta;
		RR = C/sigma/cos(theta);
		z = SS;
		zz = (z-z0)/sigma;
		temp1 = RR*RR+zz*zz;
		temp = -sigma*sigma*sigma/10.0*
			 (-zz/(8.0*exp(log(RR*RR)+4.0*log(temp1)))-
				7.0*zz/(48.0*exp(2.0*log(RR*RR)+3.0*log(temp1)))-
				35.0*zz/(192.0*exp(3.0*log(RR*RR)+2.0*log(temp1)))-
				35.0*zz/(128.0*exp(4.0*log(RR*RR)+log(temp1)))+
				35.0/128.0*(M_PI/2.0-atan2(zz,RR))/exp(9.0*log(RR))-
				signum(zz)/9.0/exp(9.0*log(fabs(zz))));
		temp -= -sigma*sigma*sigma/4.0*
				(-zz/(2.0*RR*RR*temp1)+1.0/2.0*(M_PI/2.0-atan2(zz,RR))/(RR*RR*RR)-
				1.0/3.0/(zz*zz*zz));
		z = -SS;
		zz = (z-z0)/sigma;
		temp1 = RR*RR+zz*zz;
		temp += -sigma*sigma*sigma/10.0*
				(zz/(8.0*exp(log(RR*RR)+4.0*log(temp1)))+
				7.0*zz/(48.0*exp(2.0*log(RR*RR)+3.0*log(temp1)))+
				35.0*zz/(192.0*exp(3.0*log(RR*RR)+2.0*log(temp1)))+
				35.0*zz/(128.0*exp(4.0*log(RR*RR)+log(temp1)))+
				35.0/128.0*(atan2(zz,RR)+M_PI/2.0)/exp(9.0*log(RR))+
				signum(zz)/9.0/exp(9.0*log(fabs(zz))));
		temp -= -sigma*sigma*sigma/4.0*
				(zz/(2.0*RR*RR*temp1)+1.0/2.0*(atan2(zz,RR)+M_PI/2.0)/(RR*RR*RR)+
				1.0/3.0/(zz*zz*zz));
//      if (i > 0 && i < ntheta)
//         value2 += temp*dtheta;
//      else
//         value2 += 0.5*temp*dtheta;
		if (i > 0 && i < ntheta && i%2 == 1)
		 value2 += 4.0*temp*dtheta/3.0;
		else if (i > 0 && i < ntheta)
		 value2 += 2.0*temp*dtheta/3.0;
		else
		 value2 += temp*dtheta/3.0;
	}

	theta1 = atan2(-SS+x0,SS-y0);
	theta2 = atan2(SS+x0,SS-y0);
	C = fabs(SS-y0);
	value += exp(12.0*log(sigma))/(9.0*exp(9.0*log(C)))*63.0/256.0*M_PI*
			(sin(theta2)*(exp(8.0*log(cos(theta2)))/9.0+
							8.0/9.0*(exp(6.0*log(cos(theta2)))/7.0+
									6.0/7.0*(exp(4.0*log(cos(theta2)))/5.0+
											4.0/5.0*(cos(theta2)*cos(theta2)/3.0
													 +2.0/3.0))))-
			 sin(theta1)*(exp(8.0*log(cos(theta1)))/9.0+
							8.0/9.0*(exp(6.0*log(cos(theta1)))/7.0+
									6.0/7.0*(exp(4.0*log(cos(theta1)))/5.0+
											4.0/5.0*(cos(theta1)*cos(theta1)/3.0
													 +2.0/3.0)))))-
			exp(6.0*log(sigma))/(3.0*exp(3.0*log(C)))*3.0/8.0*M_PI*
			(sin(theta2)*(cos(theta2)*cos(theta2)/3.0+2.0/3.0)-
			 sin(theta1)*(cos(theta1)*cos(theta1)/3.0+2.0/3.0));
	dtheta = (theta2-theta1)/ntheta;
	for (i = 0; i <= ntheta; i++)
	{
		theta = theta1+i*dtheta;
		RR = C/sigma/cos(theta);
		z = SS;
		zz = (z-z0)/sigma;
		temp1 = RR*RR+zz*zz;
		temp = -sigma*sigma*sigma/10.0*
			 (-zz/(8.0*exp(log(RR*RR)+4.0*log(temp1)))-
				7.0*zz/(48.0*exp(2.0*log(RR*RR)+3.0*log(temp1)))-
				35.0*zz/(192.0*exp(3.0*log(RR*RR)+2.0*log(temp1)))-
				35.0*zz/(128.0*exp(4.0*log(RR*RR)+log(temp1)))+
				35.0/128.0*(M_PI/2.0-atan2(zz,RR))/exp(9.0*log(RR))-
				signum(zz)/9.0/exp(9.0*log(fabs(zz))));
		temp -= -sigma*sigma*sigma/4.0*
				(-zz/(2.0*RR*RR*temp1)+1.0/2.0*(M_PI/2.0-atan2(zz,RR))/(RR*RR*RR)-
				1.0/3.0/(zz*zz*zz));
		z = -SS;
		zz = (z-z0)/sigma;
		temp1 = RR*RR+zz*zz;
		temp += -sigma*sigma*sigma/10.0*
				(zz/(8.0*exp(log(RR*RR)+4.0*log(temp1)))+
				7.0*zz/(48.0*exp(2.0*log(RR*RR)+3.0*log(temp1)))+
				35.0*zz/(192.0*exp(3.0*log(RR*RR)+2.0*log(temp1)))+
				35.0*zz/(128.0*exp(4.0*log(RR*RR)+log(temp1)))+
				35.0/128.0*(atan2(zz,RR)+M_PI/2.0)/exp(9.0*log(RR))+
				signum(zz)/9.0/exp(9.0*log(fabs(zz))));
		temp -= -sigma*sigma*sigma/4.0*
				(zz/(2.0*RR*RR*temp1)+1.0/2.0*(atan2(zz,RR)+M_PI/2.0)/(RR*RR*RR)+
				1.0/3.0/(zz*zz*zz));
//      if (i > 0 && i < ntheta)
//         value2 += temp*dtheta;
//      else
//         value2 += 0.5*temp*dtheta;
		if (i > 0 && i < ntheta && i%2 == 1)
		 value2 += 4.0*temp*dtheta/3.0;
		else if (i > 0 && i < ntheta)
		 value2 += 2.0*temp*dtheta/3.0;
		else
		 value2 += temp*dtheta/3.0;
	}


	theta1 = atan2(-SS+y0,SS+x0);
	theta2 = atan2(SS+y0,SS+x0);
	C = fabs(SS+x0);
	value += exp(12.0*log(sigma))/(9.0*exp(9.0*log(C)))*63.0/256.0*M_PI*
			(sin(theta2)*(exp(8.0*log(cos(theta2)))/9.0+
							8.0/9.0*(exp(6.0*log(cos(theta2)))/7.0+
									6.0/7.0*(exp(4.0*log(cos(theta2)))/5.0+
											4.0/5.0*(cos(theta2)*cos(theta2)/3.0
													 +2.0/3.0))))-
			 sin(theta1)*(exp(8.0*log(cos(theta1)))/9.0+
							8.0/9.0*(exp(6.0*log(cos(theta1)))/7.0+
									6.0/7.0*(exp(4.0*log(cos(theta1)))/5.0+
											4.0/5.0*(cos(theta1)*cos(theta1)/3.0
													 +2.0/3.0)))))-
			exp(6.0*log(sigma))/(3.0*exp(3.0*log(C)))*3.0/8.0*M_PI*
			(sin(theta2)*(cos(theta2)*cos(theta2)/3.0+2.0/3.0)-
			 sin(theta1)*(cos(theta1)*cos(theta1)/3.0+2.0/3.0));
	dtheta = (theta2-theta1)/ntheta;
	for (i = 0; i <= ntheta; i++)
	{
		theta = theta1+i*dtheta;
		RR = C/sigma/cos(theta);
		z = SS;
		zz = (z-z0)/sigma;
		temp1 = RR*RR+zz*zz;
		temp = -sigma*sigma*sigma/10.0*
			 (-zz/(8.0*exp(log(RR*RR)+4.0*log(temp1)))-
				7.0*zz/(48.0*exp(2.0*log(RR*RR)+3.0*log(temp1)))-
				35.0*zz/(192.0*exp(3.0*log(RR*RR)+2.0*log(temp1)))-
				35.0*zz/(128.0*exp(4.0*log(RR*RR)+log(temp1)))+
				35.0/128.0*(M_PI/2.0-atan2(zz,RR))/exp(9.0*log(RR))-
				signum(zz)/9.0/exp(9.0*log(fabs(zz))));
		temp -= -sigma*sigma*sigma/4.0*
				(-zz/(2.0*RR*RR*temp1)+1.0/2.0*(M_PI/2.0-atan2(zz,RR))/(RR*RR*RR)-
				1.0/3.0/(zz*zz*zz));
		z = -SS;
		zz = (z-z0)/sigma;
		temp1 = RR*RR+zz*zz;
		temp += -sigma*sigma*sigma/10.0*
				(zz/(8.0*exp(log(RR*RR)+4.0*log(temp1)))+
				7.0*zz/(48.0*exp(2.0*log(RR*RR)+3.0*log(temp1)))+
				35.0*zz/(192.0*exp(3.0*log(RR*RR)+2.0*log(temp1)))+
				35.0*zz/(128.0*exp(4.0*log(RR*RR)+log(temp1)))+
				35.0/128.0*(atan2(zz,RR)+M_PI/2.0)/exp(9.0*log(RR))+
				signum(zz)/9.0/exp(9.0*log(fabs(zz))));
		temp -= -sigma*sigma*sigma/4.0*
				(zz/(2.0*RR*RR*temp1)+1.0/2.0*(atan2(zz,RR)+M_PI/2.0)/(RR*RR*RR)+
				1.0/3.0/(zz*zz*zz));
//      if (i > 0 && i < ntheta)
//         value2 += temp*dtheta;
//      else
//         value2 += 0.5*temp*dtheta;
		if (i > 0 && i < ntheta && i%2 == 1)
		 value2 += 4.0*temp*dtheta/3.0;
		else if (i > 0 && i < ntheta)
		 value2 += 2.0*temp*dtheta/3.0;
		else
		 value2 += temp*dtheta/3.0;
	}


	theta1 = atan2(-SS-x0,SS+y0);
	theta2 = atan2(SS-x0,SS+y0);
	C = fabs(SS+y0);
	value += exp(12.0*log(sigma))/(9.0*exp(9.0*log(C)))*63.0/256.0*M_PI*
			(sin(theta2)*(exp(8.0*log(cos(theta2)))/9.0+
							8.0/9.0*(exp(6.0*log(cos(theta2)))/7.0+
									6.0/7.0*(exp(4.0*log(cos(theta2)))/5.0+
											4.0/5.0*(cos(theta2)*cos(theta2)/3.0
													 +2.0/3.0))))-
			 sin(theta1)*(exp(8.0*log(cos(theta1)))/9.0+
							8.0/9.0*(exp(6.0*log(cos(theta1)))/7.0+
									6.0/7.0*(exp(4.0*log(cos(theta1)))/5.0+
											4.0/5.0*(cos(theta1)*cos(theta1)/3.0
													 +2.0/3.0)))))-
			exp(6.0*log(sigma))/(3.0*exp(3.0*log(C)))*3.0/8.0*M_PI*
			(sin(theta2)*(cos(theta2)*cos(theta2)/3.0+2.0/3.0)-
			 sin(theta1)*(cos(theta1)*cos(theta1)/3.0+2.0/3.0));
	dtheta = (theta2-theta1)/ntheta;
	for (i = 0; i <= ntheta; i++)
	{
		theta = theta1+i*dtheta;
		RR = C/sigma/cos(theta);
		z = SS;
		zz = (z-z0)/sigma;
		temp1 = RR*RR+zz*zz;
		temp = -sigma*sigma*sigma/10.0*
			 (-zz/(8.0*exp(log(RR*RR)+4.0*log(temp1)))-
				7.0*zz/(48.0*exp(2.0*log(RR*RR)+3.0*log(temp1)))-
				35.0*zz/(192.0*exp(3.0*log(RR*RR)+2.0*log(temp1)))-
				35.0*zz/(128.0*exp(4.0*log(RR*RR)+log(temp1)))+
				35.0/128.0*(M_PI/2.0-atan2(zz,RR))/exp(9.0*log(RR))-
				signum(zz)/9.0/exp(9.0*log(fabs(zz))));
		temp -= -sigma*sigma*sigma/4.0*
				(-zz/(2.0*RR*RR*temp1)+1.0/2.0*(M_PI/2.0-atan2(zz,RR))/(RR*RR*RR)-
				1.0/3.0/(zz*zz*zz));
		z = -SS;
		zz = (z-z0)/sigma;
		temp1 = RR*RR+zz*zz;
		temp += -sigma*sigma*sigma/10.0*
				(zz/(8.0*exp(log(RR*RR)+4.0*log(temp1)))+
				7.0*zz/(48.0*exp(2.0*log(RR*RR)+3.0*log(temp1)))+
				35.0*zz/(192.0*exp(3.0*log(RR*RR)+2.0*log(temp1)))+
				35.0*zz/(128.0*exp(4.0*log(RR*RR)+log(temp1)))+
				35.0/128.0*(atan2(zz,RR)+M_PI/2.0)/exp(9.0*log(RR))+
				signum(zz)/9.0/exp(9.0*log(fabs(zz))));
		temp -= -sigma*sigma*sigma/4.0*
				(zz/(2.0*RR*RR*temp1)+1.0/2.0*(atan2(zz,RR)+M_PI/2.0)/(RR*RR*RR)+
				1.0/3.0/(zz*zz*zz));
//      if (i > 0 && i < ntheta)
//         value2 += temp*dtheta;
//      else
//         value2 += 0.5*temp*dtheta;
		if (i > 0 && i < ntheta && i%2 == 1)
		 value2 += 4.0*temp*dtheta/3.0;
		else if (i > 0 && i < ntheta)
		 value2 += 2.0*temp*dtheta/3.0;
		else
		 value2 += temp*dtheta/3.0;
	}
//   cout << "value1 = " << sol.rho0*4.0*epsilon*value << endl;
//   cout << "value2 = " << sol.rho0*4.0*epsilon*value2 << endl;
	value += value2;

	value *= sol.rho0*4.0*epsilon;
//   cout << "value4 = " << value << endl;

	return value;
}


// compute outside box LJ
double LJoutsidebox3D(SoluteStruct& sol, const GridData &grid){
	double e = 0.0;
	for (int i = 0; i < sol.dnum; i++){
	      e += outsidebox1D(i, sol, grid);
	   }
	return e;
}
