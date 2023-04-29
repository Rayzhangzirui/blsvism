#include "kernel.h"
#include <cmath>
#include "globals.h"


//kenel function, grid point at the center of the face has r=1
double getkernel(double r)
{
	if (r <= 1.0)
     // return 1.0;
		return sin(M_PI*r)*sin(M_PI*r);
//      return 1.0-r;
     	// return cos(M_PI*r)+1.0;
//      return log(4.0*exp(r))+log(2.0*exp(r-1));
	else
		return 0.0;
}

// used in getinitkernel
double getkernel(int *index, int rad, int dim)
{
	int i, r;

	r = 0;
	for (i = 0; i < dim; i++)
		r += (index[i]-rad)*(index[i]-rad);
	
	return getkernel(sqrt(static_cast<double>(r))/rad);
}

// integrate a_0(theta) from 0 to 1 
double getdenom(int N)
{
	int i;
	double thefactor, value = 0.0;
	double r, dr;

	dr = 1.0/N;
	for (i = 0; i <= N; i++)
	{
		r = i*dr;
		if (i == 0 || i == N)
		 thefactor = dr/3.0;
		else if (i%2 == 1)
		 thefactor = 4.0*dr/3.0;
		else
		 thefactor = 2.0*dr/3.0;

		value += getkernel(r)*r*r*r*thefactor;//why r^3
	}
	value *= M_PI;

	return value;
}
// simpson's rule to get denominator
double getdenom(double rad, int N)
{
	return getdenom(N)*exp(4.0*log(rad));
}


//defalut constructor
KernelStruct::KernelStruct(){};

//constructor
KernelStruct::KernelStruct(GridData &grid, int radius){
	getinitkernel2(*this, grid, radius);
}

KernelStruct::KernelStruct(GridData &grid, double radius){
	getinitkerneldouble(*this, grid, radius);
}

//destructor 
KernelStruct::~KernelStruct(){
	free_matrix<double>(this->K, 2*rad, 2*rad, 2*rad);
	free_matrix<int>(this->index,(2*rad+1)*(2*rad+1)*(2*rad+1)-1, 3-1);
}

// copy constructor 
KernelStruct::KernelStruct(const KernelStruct &rhs){
	rad = rhs.rad;
	start = rhs.start;
	N = rhs.N;
	K = matrix<double>(2*rad,2*rad,2*rad);
	for (int i = 0; i <= 2*rad; i++){
		for (int j = 0; j <= 2*rad; j++){
			for (int k = 0; k <= 2*rad; k++){
				K[i][j][k] = rhs.K[i][j][k];
			}
		}
	}
	int n = pow(2*rad+1,3);
	index = matrix<int>(n-1,2-1);
	for (int i = 0; i < n; i++){
		index[i][0] = rhs.index[i][0];
		index[i][1] = rhs.index[i][1];
		index[i][2] = rhs.index[i][2];
	}
}



void getinitkerneldouble(KernelStruct &kernel, GridData &grid, double drad){
	
	int intrad = ceil(drad/grid.dx[0]);
	kernel.rad = intrad;
	kernel.K = matrix<double>(2*intrad,2*intrad,2*intrad); //matrix of size 2*rad+1
	kernel.index = matrix<int>(pow(2*intrad+1,3),grid.dim-1); // may not use all the rows
	double thedenom = getdenom( drad, 1000);

   	kernel.N = 0;
   	int tindex[3];
	for(int i = 0; i <= 2*intrad; i++){
		for(int j = 0; j <= 2*intrad; j++){
			for(int k = 0; k <= 2*intrad; k++){
				tindex[0] = i;
				tindex[1] = j;
				tindex[2] = k;

				double r = sqrt( pow(i-intrad,2)+pow(j-intrad,2) +pow(k-intrad,2)) * grid.dx[0]; // distance to center
				
				double value = getkernel(r/drad);
				if (fabs(value) > TOL){
					kernel.index[kernel.N][0] = i;
					kernel.index[kernel.N][1] = j;
					kernel.index[kernel.N][2] = k;
					kernel.N++;
				}
				value /= thedenom;
				setvalarray(kernel.K,tindex,value);		
			}
		}
	}

}


// modify getinitkernel
void getinitkernel2(KernelStruct &kernel, GridData &grid, int rad){
     
	if(rad==0){
		// kernelrad is a multiple of h, say k, kernel raidus should scale with sqrt(h)
		// kh = C sqrt(h), therefore k = C sqrt(1/h)
		//if rad not specified, for h = 1, r = 3. r scales with sqrt(h)
		kernel.rad = (int) ceil(3*sqrt(1.0/grid.dx[0]));
	}
	else{
		kernel.rad = rad;
	}
	
	// cout<<"kernel.rad = "<<kernel.rad<<endl;

	kernel.K = matrix<double>(2*kernel.rad,2*kernel.rad,2*kernel.rad); //matrix of size 2*rad+1
	kernel.index = matrix<int>(pow(2*kernel.rad+1,3),grid.dim-1);
	double thedenom = getdenom(kernel.rad*grid.dx[0],1000);
//   thedenom = getdenom(-grid.a[0],kernel.rad,grid);
//   thedenom = getdenomkernel(kernel.rad,grid);
	// cout << "using denom" << thedenom << endl;

   	kernel.N = 0;
   	int tindex[3];
	for(int i = 0; i <= 2*kernel.rad; i++){
		for(int j = 0; j <= 2*kernel.rad; j++){
			for(int k = 0; k <= 2*kernel.rad; k++){
					tindex[0] = i;
					tindex[1] = j;
					tindex[2] = k;
					double value = getkernel(tindex, kernel.rad, grid.dim); //value of kernel at tindex
					if (fabs(value) > 1.0e-15){
						kernel.index[kernel.N][0] = i;
						kernel.index[kernel.N][1] = j;
						kernel.index[kernel.N][2] = k;
						kernel.N++;
					}
					value /= thedenom;
					setvalarray(kernel.K,tindex,value);
			}
		}
	}

  //  cout << "kernel has " << kernel.N << " non-zero, non-center elements out of " 
		// << pow(2*kernel.rad+1,3) << ", or " 
		// << static_cast<double>(kernel.N)/pow(2*kernel.rad+1,3)*100.0 << " percent" 
		// << endl;
}



// get surface area of phi
double getlengthlist(double ***phi, KernelStruct &kernel, GridData &grid)
{
	int i, r, L, tindex[grid.dim], rindex[grid.dim], sindex[grid.dim];
	double value = 0.0, Kval, rad = kernel.rad;
	double dV = 1.0;
	for (r = 0; r < grid.dim; r++)
		dV *= grid.dx[r];

	for (i = 0; i < grid.dim; i++) //tindex = current grid index
		tindex[i] = 0;
	while (tindex[0] <= grid.nx[0])
	{
		if (evalarray(phi,tindex) <= 0.0) //if inside
		{
		 for (r = 0; r < kernel.N; r++) //loop through nonzero kernel index
		 {
			for (i = 0; i < grid.dim; i++)
				sindex[i] = kernel.index[r][i]; //sindex = kernel point local index
			for (i = 0; i < grid.dim; i++)
				rindex[i] = tindex[i]+sindex[i]-kernel.rad;  //rindex = kernel point global index
			Kval = evalarray(kernel.K,sindex);
			for (i = 0; i < grid.dim && rindex[i] >= 0 && rindex[i] <= grid.nx[i]; i++); //if valid in bound, i = 2
			if (Kval != 0.0 && (i < grid.dim || evalarray(phi,rindex) > 0.0))
				value += Kval*dV*dV;
		 }
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