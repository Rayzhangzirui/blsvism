#include <omp.h>
#include "ompvism.h"


//sequential initialize GB
void ompCFA(SoluteStruct &sol, double*** GB, GridData &grid){
	#pragma omp parallel for collapse(3)
	for (int i=0;i<=grid.nx[0];i++){
		for (int j=0;j<=grid.nx[1];j++){
			for (int k=0;k<=grid.nx[2];k++){
				int tindex[3]={i,j,k};
				setvalarray(GB,tindex,GBwater(tindex,sol,grid));
			}
		}
	}
}

void ompCFAcomp(SoluteStruct &sol, double**** CFAcomp,  GridData &grid){
	#pragma omp parallel for collapse(3)
	for (int i = 0; i <= grid.nx[0]; i++){
		for (int j = 0; j <= grid.nx[1]; j++){
			for (int k = 0; k <= grid.nx[2]; k++){
				Sub tindex = {i,j,k};
				Coord x = sub2coord(tindex,grid);
				double w[3];
				double comp0 = 0;
				double comp1 = 0;
				double comp2 = 0;
				for (int r = 0; r<sol.dnum; r++){
					w[0] = x[0]-sol.d[r][0]; //w(s) = x(s)-x_i(s), s= 1 2 3 dimension		
					w[1] = x[1]-sol.d[r][1];
					w[2] = x[2]-sol.d[r][2];
					double gradw3 = pow(w[0]*w[0]+w[1]*w[1]+w[2]*w[2],1.5);//gradw3 =|x-xi|^2
					if (gradw3 < 1.0e-10){
						gradw3 = 1.0e-10;
					}
					
					comp0 += sol.Q[r]*w[0]/gradw3;
					comp1 += sol.Q[r]*w[1]/gradw3;
					comp2 += sol.Q[r]*w[2]/gradw3;
					
				}

				CFAcomp[i][j][k][0] = comp0;
				CFAcomp[i][j][k][1] = comp1;
				CFAcomp[i][j][k][2] = comp2;
			}
		}
	}
}
// even slower
// void ompCFAcomp(SoluteStruct &sol, double**** CFAcomp,  GridData &grid){
// 	#pragma omp parallel for collapse(4)
// 	for (int i = 0; i <= grid.nx[0]; i++){
// 		for (int j = 0; j <= grid.nx[1]; j++){
// 			for (int k = 0; k <= grid.nx[2]; k++){
// 				for (int r = 0; r<sol.dnum; r++){
// 					Sub tindex = {i,j,k};
// 					Coord x = sub2coord(tindex,grid);
// 					double w[3];
// 					w[0] = x[0]-sol.d[r][0]; //w(s) = x(s)-x_i(s), s= 1 2 3 dimension		
// 					w[1] = x[1]-sol.d[r][1];
// 					w[2] = x[2]-sol.d[r][2];
// 					double gradw3 = pow(w[0]*w[0]+w[1]*w[1]+w[2]*w[2],1.5);//gradw3 =|x-xi|^2
// 					if (gradw3 < 1.0e-10){
// 						gradw3 = 1.0e-10;
// 					}
// 					#pragma omp atomic update
// 					CFAcomp[i][j][k][0] += sol.Q[r]*w[0]/gradw3;
// 					#pragma omp atomic update
// 					CFAcomp[i][j][k][1] += sol.Q[r]*w[1]/gradw3;
// 					#pragma omp atomic update
// 					CFAcomp[i][j][k][2] += sol.Q[r]*w[2]/gradw3;
					
// 				}
// 			}
// 		}
// 	}
// }

//sequential initialize LJ
void ompLJ(SoluteStruct &sol,double*** LJ, GridData &grid){
	#pragma omp parallel for collapse(3)
	for(int i=0; i<=grid.nx[0]; i++){
		for(int j=0; j<=grid.nx[1]; j++){
			for(int k=0; k<=grid.nx[2]; k++){
				int tindex[3]={i,j,k};
				setvalarray(LJ,tindex,LJwater(tindex,sol,grid));
			}
		}
	}
}

//sequential initialize LJ1d
void ompLJ1d(SoluteStruct &sol,double* LJ1d, GridData &grid){
	#pragma omp parallel for
	for(int i=0; i<grid.N; i++){
		Sub tindex = ind2sub(grid,i);
		LJ1d[i] = LJwater(tindex.data(),sol,grid);
	}
}


// sequential initialize surface
void omp_initsurf(SurfaceStruct &surf, SoluteStruct &sol, GridData &grid, bool tightfit){
	#pragma omp parallel for collapse(3)
	for(int i=0; i<=grid.nx[0]; i++){
		for(int j=0; j<=grid.nx[1]; j++){
			for(int k=0; k<=grid.nx[2]; k++){
				Sub tindex = {i,j,k};
				Coord x = sub2coord(tindex,grid);
				if (tightfit){
					//if tightfit, set -1 if close to any atom
					surf.phi[i][j][k] = 1.0;
					for(int n = 0; n <sol.dnum; n++){
						if (sqrt((x[0]-sol.d[n][0])*(x[0]-sol.d[n][0])+
								 (x[1]-sol.d[n][1])*(x[1]-sol.d[n][1])+
								 (x[2]-sol.d[n][2])*(x[2]-sol.d[n][2])) <= (sol.sigma[n]+WATERSIG)/2) {
								surf.phi[i][j][k] = -1.0;
								break;
							}
						}
				}else{
					//if loosefit, set everything as inside
					surf.phi[i][j][k] = -1.0;
				}
			}
		}
	}
}


double omp_LJsolsol(SoluteStruct &Asol,SoluteStruct &Bsol){
	double temp = 0.0;
	#pragma omp parallel for collapse(2)
	for (int i=0; i<Asol.dnum; i++){
		for (int j=0; j<Bsol.dnum; j++){
			double sigma = 0.5*(Asol.sigma[i] + Bsol.sigma[j]);
			double epsilon = sqrt(Asol.epsilon[i] * Bsol.epsilon[j]);
			double r = sqrt(pow(Asol.d[i][0] - Bsol.d[j][0],2)+
							pow(Asol.d[i][1] - Bsol.d[j][1],2)+
							pow(Asol.d[i][2] - Bsol.d[j][2],2));
			#pragma omp atomic update
			temp += 4.0*epsilon*(pow(sigma/r,12.0)-pow(sigma/r,6.0));
		}
	}
	return temp;
}


double omp_outsidebox_cfa(SoluteStruct &sol, GridData &boxgrid, int ntheta, int nphi, int nr){
	//spherical coordinate range
	double atheta = 0;
	double aphi = 0;
	double ar = 0;
	double btheta = 2*M_PI;
	double bphi = M_PI;
	double br = 1.0/min(boxgrid.b[0],min(boxgrid.b[1],boxgrid.b[2]));

	double coeff = 1.0/(32.0*M_PI*M_PI*sol.epsilon0)*(1.0/sol.epsilonex-1.0/sol.epsilonin);

	GridData sphgrid(ntheta, nphi, nr, atheta, aphi, ar, btheta, bphi, br);
	double dV = sphgrid.dx[0]*sphgrid.dx[1]*sphgrid.dx[2];
	double total = 0.0;
	int N = (sphgrid.nx[0])*(sphgrid.nx[1])*(sphgrid.nx[2]);
	#pragma omp parallel for collapse(3)
	for(int i = 0; i < ntheta; i++){
		for(int j = 0; j < nphi; j++){
			for(int k = 0; k < nr; k++){
				double theta = atheta + i * sphgrid.dx[0] + sphgrid.dx[0]/2;// avoid 0 and 2pi
				double phi = aphi + j * sphgrid.dx[1] + sphgrid.dx[1]/2;//
				double r = ar + k * sphgrid.dx[2] + sphgrid.dx[2]/2; // avoid r = 0	

				double x =  sin(phi) * cos(theta)/r;
				double y =  sin(phi) * sin(theta)/r;
				double z =  cos(phi)/r;

				double vpart = 0;
				//if out of box
				if( x < boxgrid.a[0] - boxgrid.dx[0]/2 || x > boxgrid.b[0] + boxgrid.dx[0]/2 ||
					y < boxgrid.a[1] - boxgrid.dx[1]/2 || y > boxgrid.b[1] + boxgrid.dx[1]/2 ||
					z < boxgrid.a[2] - boxgrid.dx[2]/2 || z > boxgrid.b[2] + boxgrid.dx[2]/2){

					double integralfactor = sin(phi)/pow(r,4.0) * dV;
					double v1 = 0;
					double v2 = 0;
					double v3 = 0;
					for (int m = 0; m < sol.dnum; m++){
						double q = sol.Q[m];
						double w1 = x - sol.d[m][0];
						double w2 = y - sol.d[m][1];
						double w3 = z - sol.d[m][2];
						double gradw3 = pow(w1*w1 + w2*w2 + w3*w3, 1.5);
						v1 += q*w1/gradw3;
						v2 += q*w2/gradw3;
						v3 += q*w3/gradw3;
					}
					vpart = (v1*v1 + v2*v2 + v3*v3) * integralfactor;
				}

				#pragma omp atomic update
				total += vpart;
			}
		}
	}
	return coeff*total;
}


double omp_ELECsolsol(SoluteStruct &Asol,SoluteStruct &Bsol,double permittivity){
	const double epsilon0 = 0.00014319;
	double sum = 0.0;
	#pragma omp parallel for collapse(2)
	for (int i=0; i<Asol.dnum; i++){
	  for (int j=0; j<Bsol.dnum; j++){
		double r = sqrt(pow(Asol.d[i][0] - Bsol.d[j][0],2)+
						pow(Asol.d[i][1] - Bsol.d[j][1],2)+
						pow(Asol.d[i][2] - Bsol.d[j][2],2));
		#pragma omp atomic update
		sum += Asol.Q[i]*Bsol.Q[j]/r;
	  }
   }
   return sum/(4*M_PI*epsilon0*permittivity);
}


