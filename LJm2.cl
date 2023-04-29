

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
__kernel void getLJ2(__global double *LJ, int noinit, double solx, double soly, 
                    double solz, double sole, double sols, int nx, int ny, int nz, 
                    double dx, double dy, double dz, double ax, double ay, double az);

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

__kernel void getLJm2(__global double *LJ,__global double *soldata,
                    int dnum, int nx, int ny, int nz, 
                    double dx, double dy, double dz, double ax, double ay, double az)
{
   int dim = 3;
   double x[dim];
   int index[dim], Nx[dim];
   double  Dx[dim], A[dim];
   double tol = 1.0e-10;
   //average sigma and epsilon with water
   double sigma, eps;

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
      double dist_square = 0.0;//distance squred between grid and particle
      double temp = 0.0; //sum of LJ at grid
      sub2coord(x,index,Dx,A,dim);//x is coordinate of current grid point
      for (int i = 0;i<dnum;i++){
         sigma = soldata[4*dnum+i];
         eps = soldata[3*dnum+i];
         
         dist_square =  (x[0]-soldata[0*dnum+i])*(x[0]-soldata[0*dnum+i])+
                        (x[1]-soldata[1*dnum+i])*(x[1]-soldata[1*dnum+i])+
                        (x[2]-soldata[2*dnum+i])*(x[2]-soldata[2*dnum+i]);
         if (dist_square < tol){
            dist_square = tol;
      }
      
      temp = temp + 4.0*eps*(pow(sigma*sigma/dist_square,6)-pow(sigma*sigma/dist_square,3));
      }
      //printf("%d %d %d %f %f %f %f %f %f %f\n",index[0],index[1],index[2],x[0],x[1],x[2], sigma,eps,temp,dist_square);
      setvalarray(LJ,index,Nx,dim,temp);
   }

}
