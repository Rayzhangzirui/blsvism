//when gpu does not support double, but this will give large error
//#define double float
// copy from cfang.cl
void sub2coord(double *x, int *index, double *Dx, double *A, int dim)
{
   int r;

   for (r = 0; r < dim; r++)
      x[r] = A[r]+index[r]*Dx[r];
}
// copy from cfang.cl
void getcurrentindex(int *index, int dim)
{
   int i;

   for (i = 0; i < dim; i++)
      index[i] = get_global_id(i);
}


//given 3d index and axis (x=0,y=1,z=2), convert to 1D index loc, evalute 4d array A, used for GB
double evalarray4d(__global double *A, int *index, int *Nx, int axis)
{
   int dim = 3;
   int loc = index[0];
   int i, factor = 1;

   for (i = 1; i < dim; i++)
   {
      factor *= Nx[i-1]+1;
      loc += factor*index[i];    //loc = index[0] + (Nx[0]+1)*index[1] + (Nx[0]+1)*(Nx[1]+1)*index[2]
   }
   loc += axis*(Nx[0]+1)*(Nx[1]+1)*(Nx[1]+1);
   return A[loc];
}

void setvalarray4d(__global double *A, int *index, int *Nx, int axis, double value)
{
   int dim = 3;
   int loc = index[0];
   int i, factor = 1;

   for (i = 1; i < dim; i++)
   {
      factor *= Nx[i-1]+1;
      loc += factor*index[i];
   }
   loc += axis*(Nx[0]+1)*(Nx[1]+1)*(Nx[1]+1);
   A[loc] = value;
}

// v is 4d array, vector field of Qi (x-xi)/|x-xi|^3,v[0/1/2]=component in x/y/z
//noinit = index of particle, also current step number of initialization
//solx/soly/solz = coordinate of particle, sole = epsilon, sols = sigma,
//nx/ny/nz = number of grid point in x/y/z direction,
//dx/dy/dz = spacing of grid in x/y/z direction,
//ax/ay/az = left boundary in x/y/z direction
__kernel void getGBcomponent(__global double *v, int noinit, double solx, double soly, 
                  double solz, double Q, int nx, int ny, int nz, 
                  double dx, double dy, double dz, double ax, double ay, double az)
{
   int dim = 3;
   double gradw3,temp;
   int index[3], Nx[3];
   double d[3], Dx[3], A[3], x[3],w[3],v_component[3];
   double tol = 1.0e-10;
   int i;

   Nx[0] = nx;
   Nx[1] = ny;
   Nx[2] = nz;
   Dx[0] = dx;
   Dx[1] = dy;
   Dx[2] = dz;
   A[0] = ax;
   A[1] = ay;
   A[2] = az;
   d[0] = solx;
   d[1] = soly;
   d[2] = solz;
   getcurrentindex(index, 3);
   for (i = 0; i < dim && index[i] <= Nx[i]; i++);
   // if index[i]<Nx[i],i=0,1,2, the index is valid, and i = 3.
   if (i >= dim){
      sub2coord(x,index,Dx,A,dim);
      //if not the first step, get GB
      if (noinit){
         v_component[0] = evalarray4d(v,index,Nx,0);
         v_component[1] = evalarray4d(v,index,Nx,1);
         v_component[2] = evalarray4d(v,index,Nx,2);
      }
      else{
         v_component[0] = 0.0;
         v_component[1] = 0.0;
         v_component[2] = 0.0;
      }
      w[0] = x[0]-solx;
      w[1] = x[1]-soly;
      w[2] = x[2]-solz;
      gradw3 = w[0]*w[0]+w[1]*w[1]+w[2]*w[2];
      gradw3 *= sqrt(gradw3);
      if (gradw3 < 1.0e-10)
         gradw3 = 1.0e-10;
      for (int axis = 0; axis < dim; axis++){
         temp = Q*w[axis]/gradw3;
         setvalarray4d(v,index,Nx,axis,v_component[axis]+ temp);
      }
      //printf("index %i,%i,%i:gradw3 = %f:v0 = %f:v1 = %f:v2 = %f \n",
      //   index[0],index[1],index[2],gradw3,temp_debug[0],temp_debug[1],temp_debug[2]);//debug
   }
}


__kernel void getCFAcomp(__global double *v, __global double *soldata,
                  int dnum, int nx, int ny, int nz, 
                  double dx, double dy, double dz, double ax, double ay, double az)
{
   const int dim = 3;
   double gradw3,temp;
   int index[dim], Nx[dim];
   double d[dim], Dx[dim], A[dim], x[dim],w[dim], cfa_comp[dim];
   double tol = 1.0e-10;
   int i;

   Nx[0] = nx;
   Nx[1] = ny;
   Nx[2] = nz;
   Dx[0] = dx;
   Dx[1] = dy;
   Dx[2] = dz;
   A[0] = ax;
   A[1] = ay;
   A[2] = az;
   
   getcurrentindex(index,dim);
   if (index[0]<=Nx[0]&&index[1]<=Nx[1]&&index[2]<=Nx[2]){
      sub2coord(x,index,Dx,A,dim);//x is coordinate of current grid point
      cfa_comp[0] = 0.0;
      cfa_comp[1] = 0.0;
      cfa_comp[2] = 0.0;
      
      for (int i = 0;i<dnum;i++){// for each atom
         for (int d =0; d < dim; d++){
            w[d] = x[d]-soldata[d*dnum+i]; // w[d] = x[d]-xi[d]
         }

         double gradw3 =  pow(w[0]*w[0] + w[1]*w[1] + w[2]*w[2], 1.5);
         if (gradw3 < tol){
            gradw3 = tol;
         }

         for (int d =0; d < dim; d++){
            cfa_comp[d] += soldata[5*dnum + i] * w[d] / gradw3; //soldata[5*dnum + i] = sol(i).Q
         }
      }
      // set cfa_comp to 4d CFAcomp array v
      for (int d =0; d < dim; d++){
         setvalarray4d(v,index,Nx,d,cfa_comp[d]);
      }
      //printf("index %i,%i,%i:gradw3 = %f:v0 = %f:v1 = %f:v2 = %f \n",index[0],index[1],index[2],gradw3,cfa_comp[0],cfa_comp[1],cfa_comp[2]);//debug
   }
   
   
   
}



