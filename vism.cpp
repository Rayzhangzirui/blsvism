#define CL_SILENCE_DEPRECATION //opencl deprecated in macOS 10.14 Mojave

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <sstream>
#include <string>

#include "vism.h"
#include "vism_addon.h"
#include "globals.h"

#include <spdlog/spdlog.h>
#include <fmt/core.h>

using namespace std;



double EnergyChangeField(int *index, SurfaceStruct &surf, double ***field, GridData &grid){
	double dV = grid.dx[0] * grid.dx[1] * grid.dx[2];
	double part = dV * evalarray(field,index);
	part = (evalarray(surf.phi,index) > 0.0)? -part:part; 	
	return part;

}

// calculate change of surface energy if flipped index, it depends on nearby points used by kernel
double EnergyChangeSurf(int *index, SurfaceStruct &surf, SoluteStruct &sol, KernelStruct &kernel, GridData &grid)
{
	int r, s, m, n, t, tindex[grid.dim], rindex[grid.dim], sindex[grid.dim];
	double temp, grad, surfpart = 0.0, Kval, phival;
	int i, L, rad = kernel.rad;

	double dV = 1.0;
	for (int r = 0; r < grid.dim; r++)
		dV *= grid.dx[r];
	// loop through nonzero kernel element, surfpart is sum of kernel inside(phi<0) - sum of kernel outside
	// suppose index is inside (phi<0), then change in surface is inside kernel - outside kernel
	for (int r = 0; r < kernel.N; r++)
	{
		for (i = 0; i < grid.dim; i++)
		 sindex[i] = kernel.index[r][i];//sindex =  non-zero local index in kernel
		for (i = 0; i < grid.dim; i++)
		 tindex[i] = index[i]+sindex[i]-kernel.rad; // tindex = non-zero global index
		Kval = evalarray(kernel.K,sindex);
		if (Kval != 0.0)
		{
		 for (i = 0; i < grid.dim && tindex[i] >= 0 && tindex[i] <= grid.nx[i]; i++);// if tindex is in box, i = 3, otherwise i <= 2
		 if (i >= grid.dim)
			phival = evalarray(surf.phi,tindex);
		 if (i >= grid.dim && phival <= 0.0) // sum inside
			surfpart += Kval;
		 else if (i < grid.dim || phival > 0.0) // subtracted by sum outside, out of grid always counted as outside
			surfpart -= Kval; 
		}
	}
	//if current phi>0, if flipped, become inside, LJ GB reduce,
	if (evalarray(surf.phi,index) > 0.0)
	{
		surfpart = -surfpart;
	}

	surfpart *= gamma0*dV*dV;

	return surfpart;
}


// calculate change of total energy if flipped index, it depends on nearby points used by kernel
double energychangelist(int *index, SurfaceStruct &surf, SoluteStruct &sol, double ***LJ, double ***GB,
						KernelStruct &kernel, GridData &grid)
{
	int r, s, m, n, t, tindex[grid.dim], rindex[grid.dim], sindex[grid.dim];
	double temp, grad, surfpart = 0.0, LJpart = 0.0, GBpart = 0.0, dV, Kval, phival;
	int i, L, rad = kernel.rad;

	dV = 1.0;
	for (r = 0; r < grid.dim; r++)
		dV *= grid.dx[r];
	// loop through nonzero kernel element, surfpart is sum of kernel inside(phi<0) - sum of kernel outside
	for (r = 0; r < kernel.N; r++)
	{
		for (i = 0; i < grid.dim; i++)
		 sindex[i] = kernel.index[r][i];//sindex =  non-zero local index in kernel
		for (i = 0; i < grid.dim; i++)
		 tindex[i] = index[i]+sindex[i]-kernel.rad; // tindex = non-zero global index
		Kval = evalarray(kernel.K,sindex);
		if (Kval != 0.0)
		{
		 for (i = 0; i < grid.dim && tindex[i] >= 0 && tindex[i] <= grid.nx[i]; i++);
		 if (i >= grid.dim)
			phival = evalarray(surf.phi,tindex);
		 if (i >= grid.dim && phival <= 0.0) // phi<0 and inbound is inside
			surfpart += Kval;
		 else if (i < grid.dim || phival > 0.0) // phi>0 or out of bound is outside
			surfpart -= Kval; 
		}
	}

//   LJpart = LJwater(index,sol,grid);
	LJpart = evalarray(LJ,index);
	GBpart = evalarray(GB,index);

	//:if current phi>0, if flipped, become inside, LJ GB reduce,
	if (evalarray(surf.phi,index) > 0.0)
	{
		surfpart = -surfpart;
		LJpart = -LJpart;
		GBpart = -GBpart;
	
	}

	surfpart *= gamma0*dV*dV;
	LJpart *= rho0*dV;
	GBpart *= dV;

	return surfpart+LJpart+GBpart;
}

// suppose pt at tindex has changed sign recently.  want to return the difference
// in energy change if pt at index changes sign before and after the change at tindex
double energychangept2(int *index, int *tindex, SurfaceStruct &surf, SoluteStruct &sol,
						KernelStruct &kernel, GridData &grid)
{

	int r, s, m, n, t, rindex[grid.dim], sindex[grid.dim];
	double temp, grad, surfpart = 0.0, LJpart = 0.0, dV, Kval;
	int i, L, rad = kernel.rad;

	dV = 1.0;
	for (r = 0; r < grid.dim; r++)
		dV *= grid.dx[r];

	for (i = 0; i < grid.dim; i++)
		sindex[i] = index[i]-tindex[i]+rad;//sindex is index of [point at *index] in kernel of tindex
	Kval = evalarray(kernel.K,sindex);

	if (Kval != 0.0)
	{
		for (r = 0; r < grid.dim && tindex[r] >= 0 && tindex[r] <= grid.nx[r]; r++);//check out grid bound
		if ((evalarray(surf.phi,index) > 0.0 && //index on outside and tindex on inside
			(r >= grid.dim && evalarray(surf.phi,tindex) <= 0.0)) || 
			(evalarray(surf.phi,index) <= 0.0 && //index on inside and (tindex on ouside or outofbound)
			(r < grid.dim || evalarray(surf.phi,tindex) > 0.0)))
		 surfpart -= Kval*dV*dV;
		else if ((evalarray(surf.phi,index) <= 0.0 && // index and tindex both inside
				(r >= grid.dim && evalarray(surf.phi,tindex) <= 0.0)) ||
				(evalarray(surf.phi,index) > 0.0 && // index and tindex both outside (or tindex outofbound)
				(r < grid.dim || evalarray(surf.phi,tindex) > 0.0)))
		 surfpart += Kval*dV*dV;
	}

	surfpart *= 2.0;

	surfpart *= gamma0;

	return surfpart;
}


// solute water energy inside box
double binaryenergylist(double &theenergy, SurfaceStruct &surf, SoluteStruct &sol, double ***LJ, KernelStruct &kernel, GridData &grid)
{
	int i, r, s, rad = kernel.rad;
	int tindex[grid.dim], rindex[grid.dim];
	double temp, grad, surfpart = 0.0, LJpart = 0.0, GBpart = 0.0, dV;
	clock_t time1, time2;

	dV = 1.0;
	for (r = 0; r < grid.dim; r++)
		dV *= grid.dx[r];


	surfpart = getlengthlist(surf.phi,kernel,grid);


	for (i = 0; i < grid.dim; i++)
		tindex[i] = 0;
	while (tindex[0] <= grid.nx[0])
	{
		if (evalarray(surf.phi,tindex) > 0.0)
		{
		 LJpart += evalarray(LJ,tindex)*dV;
		}

		(tindex[grid.dim-1])++;
		for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
		{
		 tindex[i] = 0;
		 (tindex[i-1])++;
		}
	}

	surfpart *= gamma0;
	LJpart *= rho0;
	
	cout << "Total energy = " << surfpart+LJpart+GBpart << endl;
	cout << "   Surf = " << surfpart << endl;
	cout << "   LJ = " << LJpart << endl;
	
	theenergy = surfpart+LJpart;

	return theenergy;
}


//look at the nearby gridpoint, check if sign different, tindex is global index of grid
char nexttointerface(double ***phi, int *tindex, GridData &grid)
{
	int i, rindex[grid.dim];;
	char away = 1;
	for (i = 0; i < grid.dim; i++)
		rindex[i] = tindex[i]-1; //initialize at a corner, increase to opposite corner
	while (away && rindex[0] <= tindex[0]+1)
	{
		for (i = 0; i < grid.dim && rindex[i] >= 0 && rindex[i] <= grid.nx[i]; i++);
		if (i >= grid.dim)//if rindex not out of bound
		 if ((evalarray(phi,tindex) < 0.0)+(evalarray(phi,rindex) < 0.0) == 1) // phi(tindex) and phi(rindex) has different sign
			away = 0;//is nexttointerface
		 else;
		else if (evalarray(phi,tindex) < 0.0)//if rindex out of bound, then tindex is on boundary
		 away = 0; //is nexttointerface

		(rindex[grid.dim-1])++;
		for (i = grid.dim-1; i > 0 && rindex[i] > min(tindex[i]+1,grid.nx[i]); i--)
		{
		 rindex[i] = max(tindex[i]-1,0);
		 (rindex[i-1])++;
		}
	}

	return (!away);
}

/*

heap.loc, 3d, going from grid to element in the heap
heap.value, 3d, store du
heap.tube, 3d, keep track of interface point
*/

FlipInfo binaryflowheaplistinterfaceonly(SurfaceStruct &surf, SoluteStruct &sol, double ***LJ, double ***GB,
	KernelStruct &kernel, GridData &grid)
{
	int i, r, step, numchange, minchange = 1, numcheck[2], maxchange = 500;
	int tindex[grid.dim], rindex[grid.dim], sindex[grid.dim]; 
	int rad = kernel.rad, L, themax;
	double theenergy, energydiff, Kval;
	HeapStruct heap;
	heap.head = NULL;
	heap.tail = NULL;
	heap.loc = heapmatrix(grid.nx[0],grid.nx[1],grid.nx[2]);
	heap.value = matrix<double>(grid.nx[0],grid.nx[1],grid.nx[2]);
	heap.tube = matrix<char>(grid.nx[0],grid.nx[1],grid.nx[2]);
	HeapElt *current;
	clock_t time1, time2;

	for (i = 0; i < 2; i++)
		numcheck[i] = 0;
	
	time1 = clock();
	// binaryenergylist(theenergy,surf,sol,kernel,grid);
	// time2 = clock();
	// cout << "binaryenergy took " << static_cast<double>(time2-time1)/CLOCKS_PER_SEC
	// 	<< " SECONDS" << endl;
	// time1 = time2;
	

	for (i = 0; i < grid.dim; i++)
		tindex[i] = 0;
	while (tindex[0] <= grid.nx[0])
	{
		if (nexttointerface(surf.phi,tindex,grid))
		{
		 energydiff = energychangelist(tindex,surf,sol,LJ,GB,kernel,grid);
		 (numcheck[0])++;
		 if (energydiff < 0.0)
			addtoheap(heap,energydiff,tindex);
		 else
			setvalarray(heap.value,tindex,energydiff);
		 setvalarray(heap.tube,tindex,(char)1);
		}
		else
		 setvalarray(heap.tube,tindex,(char)0);

		(tindex[grid.dim-1])++;
		for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
		{
		 tindex[i] = 0;
		 (tindex[i-1])++;
		}
	}
	// heaptotals(heap);
	time2 = clock();
	spdlog::info("heap init took {}s.",static_cast<double>(time2-time1)/CLOCKS_PER_SEC);
	time1 = time2;

	numchange = 0;
	while (heap.head != NULL && (*(heap.head)).index[0] >= 0)
	{
		numchange++;
		current = heap.head;
		for (i = 0; i < grid.dim; i++)
		 tindex[i] = (*current).index[i]; // tindex is current index
		energydiff = evalarray(heap.value,tindex);
		// theenergy += energydiff; 
		setvalarray(surf.phi,tindex,-evalarray(surf.phi,tindex));//flip
		setvalarray(heap.value,tindex,-energydiff);//reserse energydiff
		fixheapeltdelete(heap,current);//pop current
		for (r = 0; r < kernel.N; r++) //after flipping, interface changes, energydiff of kernel nbr changes
		{
		 for (i = 0; i < grid.dim; i++)
			sindex[i] = kernel.index[r][i];
		 Kval = evalarray(kernel.K,sindex);
		 if (Kval != 0.0) // look at nonzero kernel elements
		 {
			for (i = 0; i < grid.dim; i++)
				rindex[i] = tindex[i]+sindex[i]-rad;//rindex global index of pt in kernel
			for (i = 0; i < grid.dim && rindex[i] >= 0 && rindex[i] <= grid.nx[i]; 
				 i++); //i = 3 if not out of bound
			if (i >= grid.dim)
			{
				current = evalarray(heap.loc,rindex);
				themax = 0;
				for (i = 0; i < grid.dim; i++)
					if (abs(rindex[i]-tindex[i]) > themax)
					 themax = abs(rindex[i]-tindex[i]); // the max is Linf distance between rindex and tindex
				if ((themax > 1 && evalarray(heap.tube,rindex)) || //rindex is not direct neibor, but in tube
					(themax <= 1 && nexttointerface(surf.phi,rindex,grid)))//rindex is direct neighbor
				{
					if (evalarray(heap.tube,rindex)) //if already in tube, update energydiff
					{
					 energydiff = evalarray(heap.value,rindex)+
									energychangept2(rindex,tindex,surf,sol,kernel,grid);
					 (numcheck[1])++;
					}
					else // if not in tube, calculate energy change and add to tube
					{
					 energydiff = energychangelist(rindex,surf,sol,LJ,GB,kernel,grid);
					 setvalarray(heap.tube,rindex,(char)1);
					 (numcheck[0])++;
					}


					if (energydiff < 0.0) // if energy diff is neg,
					{
					 if (current != NULL) // if rindex in heap, change value, fix heap
						fixheapeltreplace(heap,energydiff,rindex);
					 else //if rindex not in heap, add to heap
						addtoheap(heap,energydiff,rindex);
					}
					else //if energy diff is pos
					{
					 if (current != NULL) // if rindex in heap, remove from heap
						fixheapeltdelete(heap,current);
					 setvalarray(heap.value,rindex,energydiff);
					}
				}
				else // rindex not in tube, not direct neighbor
				{
					if (current != NULL)
					 fixheapeltdelete(heap,current);
					setvalarray(heap.tube,rindex,(char)0);
				}
			}
		 }
		}
	}
	// cout << "   changed " << numchange << endl;

	// cout << "   checked " << numcheck[0] << " full and " << numcheck[1] << " partial" 
	// 	<< endl;
	// cout << "   energy = " << theenergy << endl;
	time2 = clock();
	double flowtime = static_cast<double>(time2-time1)/CLOCKS_PER_SEC;
	spdlog::info("heap flow took {}s.",flowtime);
	spdlog::info("heap changed {}.",numchange);
	
	// time1 = time2;
	// binaryenergylist(theenergy,surf,sol,kernel,grid);
	// time2 = clock();
	// cout << "binaryenergy took " << static_cast<double>(time2-time1)/CLOCKS_PER_SEC
	// 	<< " SECONDS" << endl;
	// time1 = time2;

	free_heapmatrix(heap.loc, grid.nx[0],grid.nx[1],grid.nx[2]);
	free_matrix<double>(heap.value, grid.nx[0],grid.nx[1],grid.nx[2]);
	free_matrix<char>(heap.tube, grid.nx[0],grid.nx[1],grid.nx[2]);
	FlipInfo info = {numchange, flowtime};
	return info;
}



double GSurf(SurfaceStruct &surf, KernelStruct &kernel, GridData &grid){
	return gamma0 * getlengthlist(surf.phi, kernel, grid);
}


// Initialie surf, lj, cfa. Compute energy
Energy VismEnergyInit(SoluteStruct &sol, double*** lj, double ***cfa, SurfaceStruct& surf, KernelStruct& kernel, GridData& grid, bool tightfit){
	getinitsurf(surf, sol, grid, tightfit);
    cal_LJ(sol, lj, grid);
    if(cfa){
    	seqCFA(sol,cfa,grid);	
    }
    Energy e = VismEnergyOnly(sol, lj, cfa, surf, kernel, grid);
    return e;
}

// Initialize lj, cfa, surf. Then flow.
std::tuple<Energy, FlipInfo> VismEnergyInitFlow(SoluteStruct &sol, double*** lj, double ***cfa, SurfaceStruct& surf, KernelStruct& kernel, GridData& grid, bool tightfit){
	VismEnergyInit(sol, lj, cfa, surf, kernel, grid, tightfit);
    
    // flow
    FlipInfo info = binaryflowheaplistinterfaceonly(surf, sol, lj, cfa, kernel, grid);

    Energy e = VismEnergyOnly(sol, lj, cfa, surf, kernel, grid);
    return {e, info};
}

// given lj, cfa, surf, kernel, only compute the energy
Energy VismEnergyOnly(SoluteStruct &sol, double*** lj, double ***cfa, SurfaceStruct& surf, KernelStruct& kernel, GridData& grid){
	Energy e;
    e.surf = GSurf(surf,kernel,grid);
    
    if(lj){
	    e.ljswbox = OutsideIntegral(surf, lj, grid, rho0);
	    e.ljswout = LJoutsidebox3D(sol, grid); // skip outside box ljsolwat, does not change much, small (<1) compared with others 	
    }
    
    if(cfa){
    	e.eswbox = OutsideIntegral(surf, cfa, grid, 1.0);
    	e.eswout = GB2outsidebox3D(sol, grid); // skip outside box ljsolwat, does not change much, small (<1) compared with others 	
    }
    return e;
}


