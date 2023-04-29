#ifndef UTIL_H
#define UTIL_H
#define CL_SILENCE_DEPRECATION //opencl deprecated in macOS 10.14 Mojave

#include <iostream>
#include <iomanip>
#include <fstream>
#include <array>
#include <cmath>
#include <chrono>
#include <vector>

using namespace std;
using Coord = std::array<double,3>;
using Sub = std::array<int,3>;

struct GridData
{
   static const int dim = 3;
   int nx[dim] = {0,0,0};//number of cell in each dimension
   //Note: number of grid point is nx+1, indexed from 0 to nx, in ***matrix
   //CAREFUl: in currently assuming each dimensions are the same
   double dx[dim] = {0,0,0}; //cell size
   double a[dim] = {0,0,0};//min x y z
   double b[dim] = {0,0,0};// max x y z
   int N;
   double dV;
   GridData(){};
   GridData(int gridnum,double gridbound){
      nx[0] = gridnum;
      nx[1] = nx[0];
      nx[2] = nx[0];
      a[0] = -gridbound;
      a[1] = -gridbound;
      a[2] = -gridbound;
      b[0] = gridbound;
      b[1] = gridbound;
      b[2] = gridbound;
      dx[0] = (b[0]-a[0])/nx[0];
      dx[1] = (b[1]-a[1])/nx[1];
      dx[2] = (b[2]-a[2])/nx[2];
      N = (nx[0]+1)*(nx[1]+1)*(nx[2]+1);
      dV = dx[0] * dx[1] * dx[2];
   }
   // uneven grid for spherical
   GridData(int nx0, int nx1, int nx2, double a0, double a1, double a2, double b0, double b1, double b2){
      nx[0] = nx0;
      nx[1] = nx1;
      nx[2] = nx2;
      a[0] = -a0;
      a[1] = -a1;
      a[2] = -a2;
      b[0] = b0;
      b[1] = b1;
      b[2] = b2;
      dx[0] = (b[0]-a[0])/nx[0];
      dx[1] = (b[1]-a[1])/nx[1];
      dx[2] = (b[2]-a[2])/nx[2];
      N = (nx[0]+1)*(nx[1]+1)*(nx[2]+1);
      dV = dx[0] * dx[1] * dx[2];
   }

   // subscript [i,j,k] to linear index
   int sub2ind(Sub sub){
      return sub[2] + sub[1] * (nx[1]+1) + sub[0] * (nx[1]+1) * (nx[0]+1);
   }


   Sub ind2sub(int ind){
      int i = ind/ ((nx[1]+1) * (nx[0]+1));
      int j = (ind % ((nx[1]+1) * (nx[0]+1)))/(nx[1]+1);
      int k = ind % (nx[1]+1);
   return {{i,j,k}};
   }
};


struct Energy{
   double ljss = 0; // lj solute-solute
   double bond = 0; // bond
   double bend = 0; // bending
   double surf = 0; // surface
   double ljswbox = 0; // lj solute water in box
   double ljswout = 0; // lj solute water outside
   double eswbox = 0; // electrostatic solute water in box
   double eswout = 0; // electrostatic solute water outside


   double total(){
      return ljss + bond + bend + surf + ljswbox + ljswout + eswbox + eswout;
   }

   double ljsw(){
      return  ljswbox + ljswout;
   }

   double esw(){
      return eswbox + eswout;
   }

   double solsol(){
      return ljss + bond + bend;
   }

   void Print(std::ostream& out){
      out
      <<"ljss = " << ljss 
      <<", bond = " << bond 
      <<", bend = " << bend 
      <<", surf = " << surf 
      <<", ljswbox = " << ljswbox
      <<", ljswout = " << ljswout
      <<", total = " << total()
      << std::endl;
   }

   void ConcisePrint(ostream& out){
      out<<ljss <<" "<< bond <<" "<< bend <<" "
      << surf <<" "
      << ljswbox <<" "<< ljswout <<" "
      << eswbox <<" "<< eswout <<" "
      << total()<< std::endl;
   }

};



inline double signum(double x)
{
   if (x >= 0.0)
      return 1.0;
   else
      return -1.0;
}

//subscritp to coordinate
inline void sub2coord(double *x, int *index, GridData &grid)
{
   for (int r = 0; r < grid.dim; r++)
      x[r] = grid.a[r]+index[r]*grid.dx[r];
}

template <typename T>
inline T evalarray(T ***A, int i, int j, int k)
{
   return A?A[i][j][k]:0;
}

template <typename T>
inline T evalarray(T ***A, Sub sub)
{
   return A?A[sub[0]][sub[1]][sub[2]]:0;
}

template <typename T>
inline T evalarray(T ***A, int *index)
{
   return evalarray(A,index[0],index[1],index[2]);
}

template <typename T>
inline void setvalarray(T ***A, int *index, T value)
{
   A[index[0]][index[1]][index[2]] = value;
}


// 2d double array size [row+1 col+1]
template<typename T>
inline T **matrix(int row, int col){
   T **m;
   m = new T*[row+1];
   for (int i = 0; i <= row; i++){
      m[i] = new T[col+1] ();
   }
   return m;
}

// free 2d matrix for sol.d
template <typename T> 
inline void free_matrix(T **m, int row, int col)
{
   if (!m){
      return;
   }

   for (int i = 0; i <= row; i++)
   {
      delete [] m[i];
   }
   delete [] m;
}

// 3d double array size [row+1 col+1 fr+1]
// indices from 0 to row
template <typename T> 
inline T ***matrix(int row, int col, int fr)
{
   T ***m;
   m = new T**[row+1];
   
   for (int i = 0; i <= row; i++)
   {
      m[i] = new T*[col+1];
      for (int j = 0; j <= col; j++)
      {
       m[i][j] = new T[fr+1]();
      }
   }

   return m;
}


// free matrix 3d
template<typename T>
inline void free_matrix(T ***m, int row, int col, int fr)
{
   if (!m){
      return;
   }
   int i, j;

   for (i = 0; i <= row; i++){
      for (j = 0; j <= col; j++){
         delete [] (m[i][j]);
      }
      delete [] (m[i]);
   }
   delete [] m;
}


template <typename T> 
inline T ***matrix(const GridData& grid)
{
   return matrix<T>(grid.nx[0],grid.nx[1],grid.nx[2]);
}

template <typename T> 
inline void free_matrix(T ***m, const GridData& grid)
{
   return free_matrix<T>(m, grid.nx[0], grid.nx[1], grid.nx[2]);
}



template <typename T> 
inline T ****matrix4d(int row, int col, int fr)
{
   const int dim = 3;
   T**** m;
   m = new T*** [row+1];
   for (int i = 0; i <= row; i++){
      m[i] = new T** [col+1];
      for (int j = 0; j <= col; j++){
         m[i][j] = new T* [fr+1];
         for (int k = 0;k <= fr ; k++){
            m[i][j][k] = new T[dim] ();
         }
      }
   }
   return m;
}

// free 4d matrix
template<typename T>
inline void free_matrix4d(T ****m, int row, int col, int fr){
   for (int i = 0; i <= row; i++){
      for (int j = 0; j <= col; j++){
         for (int k = 0; k <= fr ;k++ ){
            delete [] (m[i][j][k]);
         }
         delete [] (m[i][j]);
      }
      delete [] (m[i]);
   }
   delete [] m;
}


// subscript (ijk) to coordinate
inline Coord sub2coord(Sub sub, const GridData& grid){
   return Coord {grid.a[0]+sub[0]*grid.dx[0],
                     grid.a[1]+sub[1]*grid.dx[1],
                     grid.a[2]+sub[2]*grid.dx[2] };
}

// index ind to subscript [i,j,k]
inline Sub ind2sub(const GridData& grid, int ind){
   int i = ind/ ((grid.nx[1]+1) * (grid.nx[0]+1));
   int j = (ind % ((grid.nx[1]+1) * (grid.nx[0]+1)))/(grid.nx[1]+1);
   int k = ind % (grid.nx[1]+1);
   return {{i,j,k}};
}




// subscript [i,j,k] to linear index
inline int sub2ind(const GridData& grid, Sub sub){
   return sub[2] + sub[1] * (grid.nx[1]+1) + sub[0] * (grid.nx[1]+1) * (grid.nx[0]+1);
}

inline vector<int> sub2ind(const GridData& grid, vector<Sub> subs){
   vector<int> v(0, subs.size());
   for(int i = 0; i < subs.size(); i ++){
    v[i] = sub2ind(grid, subs[i]);
   }
   return v;
}

//set field to a constant, used to make sol.GB=0 or sol.LJ = 1
template <typename T>
inline void set_field(const GridData &grid, T ***matrix, T phi){
   for (int i=0;i<=grid.nx[0];i++){
      for (int j=0;j<=grid.nx[1];j++){
         for (int k=0;k<=grid.nx[2];k++){
               matrix[i][j][k]=phi;
         }
      }
   }
}

// write field to file
template <typename T>
inline void write_field(string fname,const GridData &grid,T ***matrix){
   ofstream outputfile(fname);
   for (int i=0;i<=grid.nx[0];i++){
     for (int j=0;j<=grid.nx[1];j++){
       for (int k=0;k<=grid.nx[2];k++){
            outputfile << matrix[i][j][k] << " ";
       }
     }
   }
}


// print field for debug
template <typename T>
inline void print_field( T ***f,int row, int col, int fr){
   for (int i=0; i <= row; i++){
      for (int j=0; j <= col; j++){
         for (int k=0; k <= fr; k++){
            cout<<"("<<i<<","<<j<<","<<k<<") "<<f[i][j][k]<<"\n";
         }
      }
   }
}


template <typename T>
inline void print_field( T **f, int row, int col){
   for (int i=0;i<=row;i++){
      for (int j=0;j<=col;j++){
         cout<<"("<<i<<","<<j<<") "<<f[i][j]<<"\n";
      }
   }
}


//copy from f1 to f2
template <typename T>
inline void copy_field(const GridData &grid, T ***source, T ***dest){
   for (int i=0; i <= grid.nx[0]; i++){
      for (int j=0; j <= grid.nx[1]; j++){
         for (int k=0; k <= grid.nx[2]; k++){
            dest[i][j][k] = source[i][j][k];
         }
      }
   }
}


inline Sub Add(Sub a, Sub b){
  return Sub {a[0]+b[0],a[1]+b[1],a[2]+b[2]};
}

inline bool Equal(Sub a, Sub b){
  return (a[0]==b[0])&&(a[1]==b[1])&&(a[2]==b[2]);
}


//2-norm
inline double dist2(Coord x, Coord y){
   return sqrt(pow(x[0]-y[0],2) + pow(x[1]-y[1],2) + pow(x[2]-y[2],2));
}

//1-norm
inline double dist1(Coord x, Coord y){
   return abs(x[0]-y[0]) + abs(x[1]-y[1]) + abs(x[2]-y[2]);
}

//Infinity norm
inline double distInf(Coord x, Coord y){
   // return max( {abs(x[0]-y[0]) , abs(x[1]-y[1]) , abs(x[2]-y[2])} );
   return max(max( abs(x[0]-y[0]) , abs(x[1]-y[1])) , abs(x[2]-y[2]) );
}



//check if index out of bound
inline bool outofbound(int ind[], const GridData &grid){
   for (int d = 0; d < grid.dim; d++){
      if(ind[d]<0||ind[d]>grid.nx[d]){
         return true;
      }
   }
   return false;
}

inline bool outofbound(Sub ind, const GridData &grid){
   return outofbound(ind.data(),grid);
}

inline bool FileExist (const string& name) {
    ifstream f(name.c_str());
    return f.good();
}

// a timer that prints the during after getting out of scope.
// https://stackoverflow.com/questions/31391914/timing-in-an-elegant-way-in-c
class TimerScoped {

private:
    const chrono::steady_clock::time_point start;
    const string name;
    const clock_t c_start;
   

public:
    TimerScoped( const string & name ) : name( name ), start( chrono::steady_clock::now() ), c_start( clock() ) {
    }
    ~TimerScoped() {
        const auto end(chrono::steady_clock::now());
        const auto duration = chrono::duration_cast<chrono::nanoseconds>( end - start ).count() * 1e-9;

        // your_algorithm
      clock_t c_end = clock();

      double cputime =  ((double)c_end - c_start) / CLOCKS_PER_SEC;
        cout <<setprecision(6)<< name << " duration: " << duration << "/" << cputime<< "s" << std::endl;
    }
};



#endif


