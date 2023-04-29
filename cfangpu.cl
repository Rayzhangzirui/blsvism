//when gpu does not support double, but this will give large error
//#define double float

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//////


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

void getcurrentindex(int *index, int dim);
double evalarray(__global double *A, int *index, int *Nx, int dim);
char cevalarray(__global char *A, int *index, int *Nx, int dim);
int ievalarray(__global int *A, int *index, int *Nx, int dim);
void setvalarray(__global double *A, int *index, int *Nx, int dim, double value);
void csetvalarray(__global char *A, int *index, int *Nx, int dim, char value);
void isetvalarray(__global int *A, int *index, int *Nx, int dim, int value);
// HERE
void sub2coord(double *x, int *index, double *Dx, double *A, int dim);
__kernel void getLJ(__global double *LJ, int noinit, double solx, double soly, 
                    double solz, double sole, double sols, int nx, int ny, int nz, 
                    double dx, double dy, double dz, double ax, double ay, double az);
__kernel void getphi(__global double *phi, int noinit, double solx, double soly,
                     double solz, double sols, int nx, int ny, int nz, double dx,
                     double dy, double dz, double ax, double ay, double az);
__kernel void getLJsol(__global double *LJsol, __global double *psol, int noinit, 
                       double solx, double soly, double solz, double sole, 
                       double sols, int N);
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

// sets initial grid position for this work-item:

void getcurrentindex(int *index, int dim)
{
   int i;

   for (i = 0; i < dim; i++)
      index[i] = get_global_id(i);
}

////////////////////////////////////////////////////////////////////

// retrieves data from "thearray" (that represents a grid of dimensions
// "dims") at position "position":

double evalarray(__global double *A, int *index, int *Nx, int dim)
{
   int loc = index[0];
   int i, factor = 1;

   for (i = 1; i < dim; i++)
   {
      factor *= Nx[i-1]+1;
      loc += factor*index[i];
   }

   return A[loc];
}

char cevalarray(__global char *A, int *index, int *Nx, int dim)
{
   int loc = index[0];
   int i, factor = 1;

   for (i = 1; i < dim; i++)
   {
      factor *= Nx[i-1]+1;
      loc += factor*index[i];
   }

   return A[loc];
}

int ievalarray(__global int *A, int *index, int *Nx, int dim)
{
   int loc = index[0];
   int i, factor = 1;

   for (i = 1; i < dim; i++)
   {
      factor *= Nx[i-1]+1;
      loc += factor*index[i];
   }

   return A[loc];
}

////////////////////////////////////////////////////////////////////

// puts "thedata" into "thearray" (that represents a grid of dimensions "dims)
// at position "position":
// A is 1d array that flatten the 3d grid. index = {i,j,k}, 0<=i<= Nx[0],0<=j<= Nx[1],0<=k<= Nx[2]
// put value at A[i + (Nx[0]+1)j + (Nx[0]+1)(Nx[1]+1)k ]
// note Nx is number of cell, Nx+1 is number of grid point, from 0 to Nx
void setvalarray(__global double *A, int *index, int *Nx, int dim, double value)
{
   int loc = index[0];
   int i, factor = 1;

   for (i = 1; i < dim; i++)
   {
      factor *= Nx[i-1]+1;
      loc += factor*index[i];
   }

   A[loc] = value;
}

void csetvalarray(__global char *A, int *index, int *Nx, int dim, char value)
{
   int loc = index[0];
   int i, factor = 1;

   for (i = 1; i < dim; i++)
   {
      factor *= Nx[i-1]+1;
      loc += factor*index[i];
   }

   A[loc] = value;
}

void isetvalarray(__global int *A, int *index, int *Nx, int dim, int value)
{
   int loc = index[0];
   int i, factor = 1;

   for (i = 1; i < dim; i++)
   {
      factor *= Nx[i-1]+1;
      loc += factor*index[i];
   }

   A[loc] = value;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

void sub2coord(double *x, int *index, double *Dx, double *A, int dim)
{
   int r;

   for (r = 0; r < dim; r++)
      x[r] = A[r]+index[r]*Dx[r];
}

__kernel void getLJ(__global double *LJ, int noinit, double solx, double soly, 
                    double solz, double sole, double sols, int nx, int ny, int nz, 
                    double dx, double dy, double dz, double ax, double ay, double az)
{
   const int dim = 3;
   double x[3], temp, value;
   int index[3], Nx[3];
   double d[3], Dx[3], A[3];
   double tol = 1.0e-10;
   int i;
   //average sigma and epsilon with water performed outside of kernel
   
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
   getcurrentindex(index,dim);
   for (i = 0; i < dim && index[i] <= Nx[i]; i++);//if index is valid, index[i]<Nx[i] for i=1,2,3, then i=3
   if (i >= dim)
   {
      sub2coord(x,index,Dx,A,dim);
      if (noinit)
         value = evalarray(LJ,index,Nx,dim);
      else
        value = 0.0;
      temp = 0.0;
      for (i = 0; i < dim; i++)
         temp += (x[i]-d[i])*(x[i]-d[i]);
      if (temp < tol)
         temp = tol;
      //temp = 4.0*sole*(exp((double)6.0*log(sols*sols/temp))-exp((double)3.0*log(sols*sols/temp)));
      temp = 4.0*sole*(pow(sols*sols/temp,6.0)-pow(sols*sols/temp,3.0));
      setvalarray(LJ,index,Nx,dim,value+temp);
      //printf("%d %d %d %f %f %f\n",index[0],index[1],index[2],x[0],x[1],x[2]);
   }
}
// noinit should be i for i=0:sol.num-1. For the 0 th element, !noinit=true, and phi is set as 1.
// then for i=0:sol.num-1, we change phi based its distance to i-th particle

__kernel void getphi(__global double *phi, int noinit, double solx, double soly,
                     double solz, double sols, int nx, int ny, int nz, double dx,
                     double dy, double dz, double ax, double ay, double az)
{
   const int dim = 3;
   double x[3], temp;
   int index[3], Nx[3];
   double d[3], Dx[3], A[3];
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

   getcurrentindex(index,dim);
   for (i = 0; i < dim && index[i] <= Nx[i]; i++);
   if (i >= dim)
   {
      sub2coord(x,index,Dx,A,dim);

      if (!noinit)
         setvalarray(phi,index,Nx,dim,1.0);

      if (evalarray(phi,index,Nx,dim) > 0.0)
      {
         temp = 0.0;
         for (i = 0; i < dim; i++)
            temp += (x[i]-d[i])*(x[i]-d[i]);
         if (temp <= sols*sols)
            setvalarray(phi,index,Nx,dim,-1.0);
      }
   }
}


//for each atom in large protein, calculate it's LJ with all atom of small protein
//LJsol is size of larger protein, psol is smaller solute
// solx soly solz sole sols are data from psol
__kernel void getLJsol(__global double *LJsol, __global double *psol, int noinit, 
                       double solx, double soly, double solz, double sole, 
                       double sols, int N)
{
   const int dim = 3, dim2 = 2, dim1 = 1;
   double x[3], temp;
   int index[1], Nx[2], sindex[2];
   double d[3];
   double epsilon, sigma, tol = 1.0e-10;
   int i;

   getcurrentindex(index,dim1);
   if (index[0] < N)
   {
      // the first particle is not counted
      if (!noinit)
         LJsol[index[0]] = 0.0;
      //else
      //{
      d[0] = solx;
      d[1] = soly;
      d[2] = solz;

      Nx[0] = N-1;
      Nx[1] = dim-1;

      sindex[0] = index[0];
      for (i = 0; i < dim; i++)
      {
         sindex[1] = i;
         x[i] = evalarray(psol,sindex,Nx,dim2);
      }
      sindex[1] = dim;
      epsilon = sqrt(evalarray(psol,sindex,Nx,dim2)*sole);
      sindex[1] = dim+1;
      sigma = 0.5*(evalarray(psol,sindex,Nx,dim2)+sols);
      
      temp = 0.0;
      for (i = 0; i < dim; i++)
         temp += (x[i]-d[i])*(x[i]-d[i]);
      if (temp < tol)
         temp = tol;

      //temp = 4.0*epsilon*(exp((double)6.0*log(sigma*sigma/temp))-
      //                    exp((double)3.0*log(sigma*sigma/temp)));// power is 6 and 3 if temp is not sqrt
      temp = 4.0*epsilon*(pow(sigma*sigma/temp,6.0)-
                          pow(sigma*sigma/temp,3.0));
      LJsol[index[0]] += temp; 
         //printf("epsilon = %f, sigma = %f, r = %f",epsilon,sigma,sqrt(temp));
      //}
   }
}

//another way to compute LJ field
//soldata is global array for solute data, each kernel take one grid point, and calculate LJ from all atom
// nx ny nz for grid index
// dx dy dz for grid spacing
// ax ay az for grid negative end point
__kernel void getLJm2(__global double *LJ, __global double *soldata,
                    int dnum, int nx, int ny, int nz, 
                    double dx, double dy, double dz, double ax, double ay, double az)
{
   const int dim = 3;
   double x[3];
   int index[3], Nx[3];
   double  Dx[3], A[3];
   double tol = 1.0e-10;

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
      double ljtotal = 0.0; //sum of LJ at grid
      sub2coord(x,index,Dx,A,dim);//x is coordinate of current grid point
      for (int i = 0;i<dnum;i++){
         double sigma = soldata[4*dnum+i];
         double eps = soldata[3*dnum+i];
         //distance squ red between grid and particle   
         double dist_square =  (x[0]-soldata[0*dnum+i])*(x[0]-soldata[0*dnum+i])+
                        (x[1]-soldata[1*dnum+i])*(x[1]-soldata[1*dnum+i])+
                        (x[2]-soldata[2*dnum+i])*(x[2]-soldata[2*dnum+i]);
         if (dist_square < tol){
            dist_square = tol;
         }
         ljtotal = ljtotal + 4.0*eps*(pow(sigma*sigma/dist_square,6.0)-pow(sigma*sigma/dist_square,3.0));
      }
   //printf("%d %d %d %f %f %f %f %f %f %f\n",index[0],index[1],index[2],x[0],x[1],x[2], sigma,eps,ljtotal,dist_square);
   setvalarray(LJ,index,Nx,dim,ljtotal);
   }
}

__kernel void getphi2(__global double *phi, __global double *soldata, int tightfit,
                    int dnum, int nx, int ny, int nz, 
                    double dx, double dy, double dz, double ax, double ay, double az)
{
   const int dim = 3;
   double x[3];
   int index[3], Nx[3];
   double  Dx[3], A[3];
   double tol = 1.0e-10;
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
      if (tightfit==0){
         setvalarray(phi,index,Nx,dim,-1.0);//set inside
      }else{
         setvalarray(phi,index,Nx,dim,1.0);// set outside
         sub2coord(x,index,Dx,A,dim);
         for (int i = 0;i<dnum;i++){
            double sigma = soldata[4*dnum+i];
            //distance squ red between grid and particle   
            double dist_square =  (x[0]-soldata[0*dnum+i])*(x[0]-soldata[0*dnum+i])+
                           (x[1]-soldata[1*dnum+i])*(x[1]-soldata[1*dnum+i])+
                           (x[2]-soldata[2*dnum+i])*(x[2]-soldata[2*dnum+i]);
            if (dist_square < sigma*sigma){
               setvalarray(phi,index,Nx,dim,-1.0);//set inside
            }
         }   
      }
   }
}

// spherical coordinate subscript [index] to cartisian coordinate [x], r is inversed rho
//void sphsub2coord(double *sphx, int *index, double *Dx, double *A)
//{
//      sphx[0] = A[0] + index[0] * Dx[0] + Dx[0]/2;// avoid 0 and 2pi
//      sphx[1] = A[1] + index[1] * Dx[1] + Dx[1]/2;//
//      sphx[2] = A[2] + index[2] * Dx[2] + Dx[2]/2; // avoid r = 0
//
//}
//
//// spherical coordinate subscript [index] to cartisian coordinate [x], r is inversed rho
//void sph2cart(double *x, double *sphx){
//
//      x[0] =  sin(x[1]) * cos(x[0])/r;
//      x[1] =  sin(x[1]) * sin(x[0])/r;
//      x[2] =  cos(x[1])/r;
//}

// return rho^-4 in Cartesian coordinate, the integral of this function outside unit sphere should be 4pi
inline double testfcn(double x, double y, double z){
   return pow(x*x+y*y+z*z,-2.0);
}

// field is a 1 d array, dimension ntheta * nphi * nr, midpoint rule
__kernel void ousidebox_integrate(__global double *field, int nx, int ny, int nz, 
               double dtheta, double dphi, double dr,
               double ax, double ay, double az)
{
   const int dim = 3;
   int Nx[3];
   double Dx[3];
   double A[3];
   // use midpoint rule, points are at cell center
   Nx[0] = nx-1;
   Nx[1] = ny-1;
   Nx[2] = nz-1;
   Dx[0] = dtheta;
   Dx[1] = dphi;
   Dx[2] = dr;
   A[0] = ax;
   A[1] = ay;
   A[2] = az;
 
   int ntheta = get_global_id(0);
   int nphi = get_global_id(1);
   int nr = get_global_id(2);
   
   int index[3] = {ntheta, nphi, nr}; // grid index

   double theta = A[0] + index[0] * Dx[0] + Dx[0]/2;// avoid 0 and 2pi
   double phi = A[1] + index[1] * Dx[1] + Dx[1]/2;//
   double r = A[2] + index[2] * Dx[2] + Dx[2]/2; // avoid r = 0
   
   double x =  sin(phi) * cos(theta)/r;
   double y =  sin(phi) * sin(theta)/r;
   double z =  cos(phi)/r;
   // if x in side box
//   if (x[0]>ax && x[0]<-ax && x[1]>ay && x[1]<-ay && x[2]>ax && x[2]<-ax){
//      v = 0.0;
//   } else{
//      v = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
//   }
   double v ;
   v = testfcn(x,y,z)*sin(phi)/pow(r,4.0) * dtheta* dphi * dr; //r^2
   setvalarray(field, index, Nx, dim, v);
//   printf("ntheta = %d, nphi = %d, nr = %d, theta = %f, phi = %f, r = %f, x = %f, y = %f,z = %f, v =  %f\n",
//      ntheta,nphi,nr,theta,phi,r, x, y,z, v );
}

//helper function to get data from 1d soldata
inline double atomx(__global double *soldata, int dnum, int i){
   return soldata[0*dnum+i];
}

inline double atomy(__global double *soldata, int dnum, int i){
   return soldata[1*dnum+i];
}

inline double atomz(__global double *soldata, int dnum, int i){
   return soldata[2*dnum+i];
}

inline double atomq(__global double *soldata, int dnum, int i){
   return soldata[5*dnum+i];
}

// field is a 1 d array, dimension ntheta * nphi * nr, midpoint rule
__kernel void getousideboxGB(__global double *field, __global double *soldata,
               int dnum, int nx, int ny, int nz, 
               double dtheta, double dphi, double dr,
               double atheta, double aphi, double ar,
               double ax, double ay, double az)
{
   const int dim = 3;
   int Nx[3];
   double Dx[3];
   double A[3];
   Nx[0] = nx-1;
   Nx[1] = ny-1;
   Nx[2] = nz-1;
   Dx[0] = dtheta;
   Dx[1] = dphi;
   Dx[2] = dr;
   A[0] = ax;
   A[1] = ay;
   A[2] = az;
 
   int ntheta = get_global_id(0);
   int nphi = get_global_id(1);
   int nr = get_global_id(2);
   int index[3] = {ntheta, nphi, nr}; // spherical grid index

   double theta = atheta + ntheta * dtheta + dtheta/2;// avoid 0 and 2pi
   double phi = aphi + nphi * dphi + dphi/2;//
   double r = ar + nr * dr + dr/2; // avoid r = 0
   
   double x =  sin(phi) * cos(theta)/r;
   double y =  sin(phi) * sin(theta)/r;
   double z =  cos(phi)/r;

   double integralfactor = sin(phi)/pow(r,4.0) * dtheta* dphi * dr;

   double v1 = 0;
   double v2 = 0;
   double v3 = 0;
   double totalv;
   // if inside box
   if (x>ax && x<-ax && y>ay && y<-ay && z>az && z<-az){
      totalv = 0.0;
   } else{
      for (int i = 0; i<dnum;i++){
         double q = atomq(soldata, dnum, i);
         double w1 = x - atomx(soldata, dnum, i);
         double w2 = y - atomy(soldata, dnum, i);
         double w3 = z - atomz(soldata, dnum, i);
         double gradw3 = pow(w1*w1 + w2*w2 + w3*w3, 1.5);
         v1 += q*w1/gradw3;
         v2 += q*w2/gradw3;
         v3 += q*w3/gradw3;
      }
      totalv = (v1*v1 + v2*v2 + v3*v3);
   }
   totalv *= integralfactor;

   setvalarray(field, index, Nx, dim, totalv);
//   printf("ntheta = %d, nphi = %d, nr = %d, theta = %f, phi = %f, r = %f, x = %f, y = %f,z = %f, totalv =  %f\n",
//      ntheta,nphi,nr,theta,phi,r, x, y,z, totalv);
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
   int index[3], Nx[3];
   double Dx[3], A[3], x[3],w[3], cfa_comp[3];
   double tol = 1.0e-10;

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
