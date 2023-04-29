#include "kernel2.h"

#ifdef OPENMP
#include <omp.h>
#endif


const double A0_SIN2 = (M_PI * M_PI - 3) / (8 * M_PI); // a0 for sin(pi r)^2
const function<double(double)> KERNEL_SIN2 = [](double r)->double {return (r<=1.0)? sin(M_PI*r)*sin(M_PI*r): 0.0;};

const double A0_IND = M_PI/4; // a0 for indicator B(0)
const function<double(double)> KERNEL_IND = [](double r)->double {return (r<=1.0)? 1.0: 0.0;};


const double A0_COS = 12/pow(M_PI,3) - 3/M_PI + M_PI/4; // a0 for indicator B(0)
const function<double(double)> KERNEL_COS = [](double r)->double {return (r<=1.0)? cos(M_PI*r)+1: 0.0;};

// vanishing second moment
// see mathematica
const double A0_EXP2 = M_PI * 0.0925289228938235; // a0 for indicator B(0)
const function<double(double)> KERNEL_EXP2 = [](double r)->double {return (r<=1.0)? exp(2.0/(pow(2.0*r-1,2)-1.0))*(3196.1015220946833*r*r - 3457.6211113812255*r + 852.9832518883903): 0.0;};

const double A0_EXP1 = M_PI * (-0.2393733599188933); // a0 for indicator B(0)
const function<double(double)> KERNEL_EXP1 = [](double r)->double {return (r<=1.0)? exp(2.0/(pow(2.0*r-1,2)-1.0))*(-261.5195892865372*r+145.7876577089403): 0.0;};



va::Vector Sub2Vector(Sub index, const GridData& grid){
   return va::Vector(grid.a[0]+index[0]*grid.dx[0],
                     grid.a[1]+index[1]*grid.dx[1],
                     grid.a[2]+index[2]*grid.dx[2]);
}

// initialize as multiples of dx
Kernel::Kernel(GridData &grid, int multiple_dx, int code){
	InitKernelCode(*this, grid, multiple_dx * grid.dx[0]+1e-12, code);
};

// initialize as multiples of dx
Kernel::Kernel(GridData &grid, double rad, int code){
	InitKernelCode(*this, grid, rad, code);
};

// set kernel to zero, for debugging
void Kernel::set_zero(){
	index.clear();
	value.clear();
	index.push_back({0,0,0});
	value.push_back(0);
};


Kernel::~Kernel(){
	free_matrix<double>(K, N, N, N);
}

// initialize kernel, rad is absolute radiuf of the kernel
void InitKernelDoubleRadius(Kernel &kernel, GridData &grid, double rad, function<double(double)> kernel_fcn, double a0){


	kernel.index.clear();
	kernel.value.clear();
	
	kernel.radius = rad;	
	int range = ceil(rad/grid.dx[0]); // integer of half the box side
	kernel.rad = range;

	kernel.N = 2 * range;
	kernel.K = matrix<double>(kernel.N, kernel.N, kernel.N);
	
	double thedenom = a0*pow(kernel.radius, grid.dim + 1);
   	
	for(int i = -range; i <= range; i++){
		for(int j = -range; j <= range; j++){
			for(int k = -range; k <= range; k++){
				Sub tindex = {i,j,k};
				double r = sqrt( i*i + j*j + k*k ) * grid.dx[0]; // distance to center
				double sr = r / rad; // scaled radius
				
				if(sr <= 1.0){
					double value = kernel_fcn(sr)/ thedenom;
					kernel.index.push_back(tindex);
					kernel.value.push_back(value);
					kernel.K[i+range][j+range][k+range] = value;
				}
			}
		}
	}
	
}


void InitKernelCode(Kernel &kernel, GridData &grid, double rad, int code){
/* 
code = 0 indicator
code = 1 sin squred
code = 2 cos 
code = 3 gauss
*/ 
	function<double(double)> kernel_fcn;
	double a0;
	switch(code){
		case 0 :
			cout<<"indicator kernel" <<endl;
			kernel_fcn = KERNEL_IND;
			a0 = A0_IND;
			break;
		case 1 :
			cout<<"sin squred kernel" <<endl;
			kernel_fcn = KERNEL_SIN2;
			a0 = A0_SIN2;
			break;
		case 2 :
			cout<<"cos kernel" <<endl;
			kernel_fcn = KERNEL_COS;
			a0 = A0_COS;
			break;
		case 3 :
			cout<<"exp2 kernel" <<endl;
			kernel_fcn = KERNEL_EXP2;
			a0 = A0_EXP2;
			break;
		case 4 :
			cout<<"exp1 kernel" <<endl;
			kernel_fcn = KERNEL_EXP1;
			a0 = A0_EXP1;
			break;

		default:
			cerr<<"unknown kernel"<<endl;
			exit(1);
	}

	InitKernelDoubleRadius(kernel, grid, rad, kernel_fcn, a0);
}




double ComputeAreaGrid(double ***phi, Kernel &kernel, GridData &grid)
{
	double value = 0.0;
	double dV = 1.0;
	int debugcount = 0; // counter for debuging message
	
	for (int r = 0; r < grid.dim; r++)
		dV *= grid.dx[r];

	#pragma omp parallel for collapse(3)
	for(int i = 0; i <= grid.nx[0]; i++){
		for(int j = 0; j <= grid.nx[1]; j++){
			for(int k = 0; k <= grid.nx[2]; k++){
				Sub tindex = {i, j, k};
				

				if (evalarray(phi,tindex.data()) <= 0.0){ //if inside
					for(int r = 0; r < kernel.index.size(); r ++){
						Sub sindex = Add(kernel.index[r],tindex);
						if( outofbound(sindex,grid) || evalarray(phi,sindex.data()) > 0.0){
							// okay if sindex is out of bound
							#pragma omp atomic update
							value += kernel.value[r] * dV * dV;
							
							// #ifdef DEBUG
							// if (debugcount < 5){
							// 	printf("%i, %i, %i, %f\n",i,j,k, kernel.value[r] );	
							// 	debugcount ++;
							// }
							// #endif
						}
					}
				}

			}
		}
	}
	return value;
}



// calculate area of large grid number,
double ComputeAreaFcn(const Kernel &kernel, const GridData &grid,const function<bool(va::Vector)>& bsurface)
{
	double value = 0.0;
	double dV = 1.0;
	int debugcount = 0; // counter for debuging message
	
	for (int r = 0; r < grid.dim; r++)
		dV *= grid.dx[r];

	#pragma omp parallel for collapse(3)
	for(int i = 0; i <= grid.nx[0]; i++){
		for(int j = 0; j <= grid.nx[1]; j++){
			for(int k = 0; k <= grid.nx[2]; k++){
				Sub tindex = {i, j, k};

				va::Vector x = Sub2Vector(tindex, grid);

				if (!bsurface(x)){ 
					//if inside
					for(int r = 0; r < kernel.index.size(); r ++){
						Sub sindex = Add(kernel.index[r],tindex);
						va::Vector y = Sub2Vector(sindex, grid);

						if( bsurface(y)){
						// if outside
							#pragma omp atomic update
							value += kernel.value[r] * dV * dV;
							
							#ifdef DEBUG
							if (debugcount < 5){
								printf("%i, %i, %i, %f\n",i,j,k, kernel.value[r] );	
								debugcount ++;
							}
							#endif
						}
					}
				}

			}
		}
	}
	return value;
}


// calculate area of large grid number,
double ComputeSurfaceIntegral(const Kernel &kernel, const GridData &grid, const function<bool(va::Vector)>& bsurface, const function<double(va::Vector)>& fcn)
{
	double value = 0.0;
	double dV = 1.0;
	int debugcount = 0; // counter for debuging message
	
	for (int r = 0; r < grid.dim; r++)
		dV *= grid.dx[r];

	#pragma omp parallel for collapse(3)
	for(int i = 0; i <= grid.nx[0]; i++){
		for(int j = 0; j <= grid.nx[1]; j++){
			for(int k = 0; k <= grid.nx[2]; k++){
				Sub tindex = {i,j,k};

				va::Vector x = Sub2Vector(tindex, grid);

				if (!bsurface(x)){ 
					//if inside
					for(int r = 0; r < kernel.index.size(); r ++){
						Sub sindex = Add(kernel.index[r],tindex); // local sub 
						va::Vector y = Sub2Vector(sindex, grid); // global coord

						if( bsurface(y)){
						// if outside
							#pragma omp atomic update
							value += kernel.value[r] * (fcn(x)+fcn(y))/2 * dV * dV;

							#ifdef DEBUG
							if (debugcount < 3){
								printf("%i, %i, %i, %f, %f\n", i,j,k, kernel.value[r],value);	
								debugcount ++;
							}
							#endif
						}
					}
				}

			}
		}
	}
	return value;
}

// compute integrate a0(theta) from 0 to 1 using midpoint rul with n number of grid 
double ComputeA0Midpoint(int n, const function<double(double)>& kernel_fcn){
	double h = 1.0/(double)n;
	double dv = pow(h,3);
	double value = 0;
	va::Vector o(0,0,0);
	// #pragma omp parallel for collapse(4)
	for(double theta = h/2; theta < 1.0; theta += h){
		for(double i = theta+h/2; i < 1.0; i += h){
			for(double j = -1.0+h/2; j < 1.0; j += h){
				for(double k = -1.0+h/2; k < 1.0; k += h){
					va::Vector x(i,j,k);
					double r = norm( x - o);
					#pragma omp atomic update
					value += kernel_fcn(r);
				}
			}
		}
	}
	value *= dv*h;
	return value;
}


// compute denominator of kernel
double ComputeDenomMidpoint(int n, double kernel_real_radius, const GridData& grid, function<double(double)> kernel_fcn){
	double a0 = ComputeA0Midpoint(n,kernel_fcn);
	return a0*pow(kernel_real_radius, grid.dim + 1);
}




// // Compute area by conv kernel, 
// double ComputeAreaConv(double ***phi, Kernel &kernel, GridData &grid){
// 	double value = 0.0;

// 	#pragma omp parallel for collapse(3)
// 	for(int i = 0; i <= grid.nx[0]; i++){
// 		for(int j = 0; j <= grid.nx[1]; j++){
// 			for(int k = 0; k <= grid.nx[2]; k++){
// 				Sub tindex = {i, j, k};
				
// 				if (evalarray(phi,tindex.data()) <= 0.0){ //if inside
// 					for(int r = 0; r < kernel.index.size(); r ++){
// 						Sub sindex = Add(kernel.index[r],tindex);
// 						if( outofbound(sindex,grid) || evalarray(phi,sindex.data()) > 0.0){
// 							// okay if sindex is out of bound
// 							#pragma omp atomic update
// 							value += kernel.value[r] * grid.dV * grid.dV;
							
// 							// #ifdef DEBUG
// 							// if (debugcount < 5){
// 							// 	printf("%i, %i, %i, %f\n",i,j,k, kernel.value[r] );	
// 							// 	debugcount ++;
// 							// }
// 							// #endif
// 						}
// 					}
// 				}

// 			}
// 		}
// 	}
// 	return value;

// }