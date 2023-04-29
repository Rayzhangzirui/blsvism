#include "cfa.h"


double GBwater(int *index, SoluteStruct &sol, GridData &grid)
{
	double val = 0.0;
	double gradw3, x[grid.dim], v[grid.dim], w[grid.dim];

	sub2coord(x,index,grid);

	for (int s = 0; s < grid.dim; s++)
		v[s] = 0.0;
	for (int r = 0; r < sol.dnum; r++)
	{
		for (int s = 0; s < grid.dim; s++)
			w[s] = x[s]-sol.d[r][s]; //w(s) = x(s)-x_i(s), s= 1 2 3 dimension
		gradw3 = w[0]*w[0]+w[1]*w[1]+w[2]*w[2];//gradw3 =|x-xi|^2
		gradw3 *= sqrt(gradw3);//gradw3 =|x-xi|^3
		if (gradw3 < TOL)
			gradw3 = TOL;
		for (int s = 0; s < grid.dim; s++)
			v[s] += sol.Q[r]*w[s]/gradw3;
	}

//   return -22.113*(1.0/sol.epsilonin-1.0/sol.epsilonex)*
//          (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	return 1.0/(32.0*M_PI*M_PI*epsilon0)*(1.0/epsilonex-1.0/epsilonin)*
			(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}


//sequential initialize GB
void seqCFA(SoluteStruct &sol, double*** GB, GridData &grid){
	#pragma omp parallel for collapse(3)
	for (int i=0;i<=grid.nx[0];i++){
		for (int j=0;j<=grid.nx[1];j++){
			for (int k=0;k<=grid.nx[2];k++){
				int tindex[3] = {i,j,k};
				setvalarray(GB,tindex,GBwater(tindex,sol,grid));
			}
		}
	}
}


double GB2outsidebox3D(const SoluteStruct& sol, GridData &grid)
{
	int nrho = 40;
	int ntheta = 40;
	int nphi = 40;
	int ntheta4 = ntheta/4;// ntheta4=10
	double drho;
	double dtheta = 2.0*M_PI/ntheta;// dtheta = 2pi/40 
	double dphi = M_PI/nphi; // dphi = pi/40
	double SS = (grid.b[0]-grid.a[0])/2.0+0.5*grid.dx[0];// SS = 0.5 * gridsize + box boundary
//   double SS = (grid.b[0]-grid.a[0])/2.0;
	int nx = 40;
	int ny = 40;
	int nu = 40;
	double dx = 2.0*SS/nx;
	double dy = 2.0*SS/ny;
	double du;
	double phi, theta, rho, v[grid.dim], w[grid.dim], gradw3, C, x, y, u, theta1, theta2;
	double value = 0.0, temp, coef;
	double midpt[grid.dim];
	int i, j, k, r, s;

	for (i = 0; i < grid.dim; i++)
		midpt[i] = (grid.a[i]+grid.b[i])/2.0;
	coef = 1.0/(32.0*M_PI*M_PI*sol.epsilon0)*(1.0/sol.epsilonex-1.0/sol.epsilonin);
	theta1 = -M_PI/4.0; //-45 degree
	theta2 = M_PI/4.0; //+45 degree
	C = SS;
	for (i = 0; i <= nphi; i++)
	{
		phi = i*dphi;
		for (j = 0; j <= ntheta4; j++)
		{
		 theta = theta1+j*dtheta;
		 drho = cos(theta)*sin(phi)/C/nrho;
		 for (k = 0; k <= nrho; k++)
		 {
			rho = k*drho;
			v[0] = 0.0;
			v[1] = 0.0;
			v[2] = 0.0;
			for (s = 0; s < sol.dnum; s++)
			{
				w[0] = cos(theta)*sin(phi)-rho*(sol.d[s][0]-midpt[0]);
				w[1] = sin(theta)*sin(phi)-rho*(sol.d[s][1]-midpt[1]);
				w[2] = cos(phi)-rho*(sol.d[s][2]-midpt[2]);
				gradw3 = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2])*
						(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
				v[0] += sol.Q[s]*w[0]/gradw3;
				v[1] += sol.Q[s]*w[1]/gradw3;
				v[2] += sol.Q[s]*w[2]/gradw3;
			}
			temp = (v[0]*v[0]+v[1]*v[1]+v[2]*v[2])*sin(phi);
			if (i > 0 && i < nphi && i%2 == 1)
				temp *= 4.0/3.0;
			else if (i > 0 && i < nphi)
				temp *= 2.0/3.0;
			else
				temp *= 1.0/3.0;
			if (j > 0 && j < ntheta4 && j%2 == 1)
				temp *= 4.0/3.0;
			else if (j > 0 && j < ntheta4)
				temp *= 2.0/3.0;
			else
				temp *= 1.0/3.0;
			if (k > 0 && k < nrho && k%2 == 1)
				temp *= 4.0/3.0;
			else if (k > 0 && k < nrho)
				temp *= 2.0/3.0;
			else
				temp *= 1.0/3.0;
			value += coef*temp*drho*dtheta*dphi;
		 }
		}
	}

	theta1 = M_PI/4.0;
	theta2 = 3.0*M_PI/4.0;
	C = SS;
	for (i = 0; i <= nphi; i++)
	{
		phi = i*dphi;
		for (j = 0; j <= ntheta4; j++)
		{
		 theta = theta1+j*dtheta;
		 drho = sin(theta)*sin(phi)/C/nrho;
		 for (k = 0; k <= nrho; k++)
		 {
			rho = k*drho;
			v[0] = 0.0;
			v[1] = 0.0;
			v[2] = 0.0;
			for (s = 0; s < sol.dnum; s++)
			{
				w[0] = cos(theta)*sin(phi)-rho*(sol.d[s][0]-midpt[0]);
				w[1] = sin(theta)*sin(phi)-rho*(sol.d[s][1]-midpt[1]);
				w[2] = cos(phi)-rho*(sol.d[s][2]-midpt[2]);
				gradw3 = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2])*
						(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
				v[0] += sol.Q[s]*w[0]/gradw3;
				v[1] += sol.Q[s]*w[1]/gradw3;
				v[2] += sol.Q[s]*w[2]/gradw3;
			}
			temp = (v[0]*v[0]+v[1]*v[1]+v[2]*v[2])*sin(phi);
			if (i > 0 && i < nphi && i%2 == 1)
				temp *= 4.0/3.0;
			else if (i > 0 && i < nphi)
				temp *= 2.0/3.0;
			else
				temp *= 1.0/3.0;
			if (j > 0 && j < ntheta4 && j%2 == 1)
				temp *= 4.0/3.0;
			else if (j > 0 && j < ntheta4)
				temp *= 2.0/3.0;
			else
				temp *= 1.0/3.0;
			if (k > 0 && k < nrho && k%2 == 1)
				temp *= 4.0/3.0;
			else if (k > 0 && k < nrho)
				temp *= 2.0/3.0;
			else
				temp *= 1.0/3.0;
			value += coef*temp*drho*dtheta*dphi;
		 }
		}
	}

	theta1 = 3.0*M_PI/4.0;
	theta2 = 5.0*M_PI/4.0;
	C = -SS;
	for (i = 0; i <= nphi; i++)
	{
		phi = i*dphi;
		for (j = 0; j <= ntheta4; j++)
		{
		 theta = theta1+j*dtheta;
		 drho = cos(theta)*sin(phi)/C/nrho;
		 for (k = 0; k <= nrho; k++)
		 {
			rho = k*drho;
			v[0] = 0.0;
			v[1] = 0.0;
			v[2] = 0.0;
			for (s = 0; s < sol.dnum; s++)
			{
				w[0] = cos(theta)*sin(phi)-rho*(sol.d[s][0]-midpt[0]);
				w[1] = sin(theta)*sin(phi)-rho*(sol.d[s][1]-midpt[1]);
				w[2] = cos(phi)-rho*(sol.d[s][2]-midpt[2]);
				gradw3 = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2])*
						(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
				v[0] += sol.Q[s]*w[0]/gradw3;
				v[1] += sol.Q[s]*w[1]/gradw3;
				v[2] += sol.Q[s]*w[2]/gradw3;
			}
			temp = (v[0]*v[0]+v[1]*v[1]+v[2]*v[2])*sin(phi);
			if (i > 0 && i < nphi && i%2 == 1)
				temp *= 4.0/3.0;
			else if (i > 0 && i < nphi)
				temp *= 2.0/3.0;
			else
				temp *= 1.0/3.0;
			if (j > 0 && j < ntheta4 && j%2 == 1)
				temp *= 4.0/3.0;
			else if (j > 0 && j < ntheta4)
				temp *= 2.0/3.0;
			else
				temp *= 1.0/3.0;
			if (k > 0 && k < nrho && k%2 == 1)
				temp *= 4.0/3.0;
			else if (k > 0 && k < nrho)
				temp *= 2.0/3.0;
			else
				temp *= 1.0/3.0;
			value += coef*temp*drho*dtheta*dphi;
		 }
		}
	}

	theta1 = -3.0*M_PI/4.0;
	theta2 = -M_PI/4.0;
	C = -SS;
	for (i = 0; i <= nphi; i++)
	{
		phi = i*dphi;
		for (j = 0; j <= ntheta4; j++)
		{
		 theta = theta1+j*dtheta;
		 drho = sin(theta)*sin(phi)/C/nrho;
		 for (k = 0; k <= nrho; k++)
		 {
			rho = k*drho;
			v[0] = 0.0;
			v[1] = 0.0;
			v[2] = 0.0;
			for (s = 0; s < sol.dnum; s++)
			{
				w[0] = cos(theta)*sin(phi)-rho*(sol.d[s][0]-midpt[0]);
				w[1] = sin(theta)*sin(phi)-rho*(sol.d[s][1]-midpt[1]);
				w[2] = cos(phi)-rho*(sol.d[s][2]-midpt[2]);
				gradw3 = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2])*
						(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
				v[0] += sol.Q[s]*w[0]/gradw3;
				v[1] += sol.Q[s]*w[1]/gradw3;
				v[2] += sol.Q[s]*w[2]/gradw3;
			}
			temp = (v[0]*v[0]+v[1]*v[1]+v[2]*v[2])*sin(phi);
			if (i > 0 && i < nphi && i%2 == 1)
				temp *= 4.0/3.0;
			else if (i > 0 && i < nphi)
				temp *= 2.0/3.0;
			else
				temp *= 1.0/3.0;
			if (j > 0 && j < ntheta4 && j%2 == 1)
				temp *= 4.0/3.0;
			else if (j > 0 && j < ntheta4)
				temp *= 2.0/3.0;
			else
				temp *= 1.0/3.0;
			if (k > 0 && k < nrho && k%2 == 1)
				temp *= 4.0/3.0;
			else if (k > 0 && k < nrho)
				temp *= 2.0/3.0;
			else
				temp *= 1.0/3.0;
			value += coef*temp*drho*dtheta*dphi;
		 }
		}
	}

	C = SS;
	du = 1.0/C/nu;
	for (i = 0; i <= nx; i++)
	{
		x = -SS+i*dx;
		for (j = 0; j <= ny; j++)
		{
		 y = -SS+j*dy;
		 for (k = 0; k <= nu; k++)
		 {
			u = k*du;
			v[0] = 0.0;
			v[1] = 0.0;
			v[2] = 0.0;
			for (s = 0; s < sol.dnum; s++)
			{
				w[0] = u*(x-(sol.d[s][0]-midpt[0]));
				w[1] = u*(y-(sol.d[s][1]-midpt[1]));
				w[2] = 1.0-u*(sol.d[s][2]-midpt[2]);
				gradw3 = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2])*
						(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
				v[0] += u*sol.Q[s]*w[0]/gradw3;
				v[1] += u*sol.Q[s]*w[1]/gradw3;
				v[2] += u*sol.Q[s]*w[2]/gradw3;
			}
			temp = (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
			if (i > 0 && i < nx && i%2 == 1)
				temp *= 4.0/3.0;
			else if (i > 0 && i < nx)
				temp *= 2.0/3.0;
			else
				temp *= 1.0/3.0;
			if (j > 0 && j < ny && j%2 == 1)
				temp *= 4.0/3.0;
			else if (j > 0 && j < ny)
				temp *= 2.0/3.0;
			else
				temp *= 1.0/3.0;
			if (k > 0 && k < nu && k%2 == 1)
				temp *= 4.0/3.0;
			else if (k > 0 && k < nu)
				temp *= 2.0/3.0;
			else
				temp *= 1.0/3.0;
			value += coef*temp*dx*dy*du;
		 }
		}
	}

	C = -SS;
	du = -1.0/C/nu;
	for (i = 0; i <= nx; i++)
	{
		x = -SS+i*dx;
		for (j = 0; j <= ny; j++)
		{
		 y = -SS+j*dy;
		 for (k = 0; k <= nu; k++)
		 {
			u = 1.0/C+k*du;
			v[0] = 0.0;
			v[1] = 0.0;
			v[2] = 0.0;
			for (s = 0; s < sol.dnum; s++)
			{
				w[0] = u*(x-(sol.d[s][0]-midpt[0]));
				w[1] = u*(y-(sol.d[s][1]-midpt[1]));
				w[2] = 1.0-u*(sol.d[s][2]-midpt[2]);
				gradw3 = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2])*
						(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
				v[0] += u*sol.Q[s]*w[0]/gradw3;
				v[1] += u*sol.Q[s]*w[1]/gradw3;
				v[2] += u*sol.Q[s]*w[2]/gradw3;
			}
			temp = (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
			if (i > 0 && i < nx && i%2 == 1)
				temp *= 4.0/3.0;
			else if (i > 0 && i < nx)
				temp *= 2.0/3.0;
			else
				temp *= 1.0/3.0;
			if (j > 0 && j < ny && j%2 == 1)
				temp *= 4.0/3.0;
			else if (j > 0 && j < ny)
				temp *= 2.0/3.0;
			else
				temp *= 1.0/3.0;
			if (k > 0 && k < nu && k%2 == 1)
				temp *= 4.0/3.0;
			else if (k > 0 && k < nu)
				temp *= 2.0/3.0;
			else
				temp *= 1.0/3.0;
			value += coef*temp*dx*dy*du;
		 }
		}
	}

//   cout << "value = " << value << endl;
	return value;
}


void insidebox(double R, double SS, int N, SoluteStruct &sol);
