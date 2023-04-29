#ifndef SURF_H
#define SURF_H

#include "util.h"
#include "globals.h"


struct SurfaceStruct
{
   double ***phi; //phi = -1 inside, phi = 1 outside
   double p = 0;
   double delta;
   GridData grid;
   SurfaceStruct();//default constructor
   SurfaceStruct(const GridData &grid, double gamma0 = ::gamma0);//constructor
   SurfaceStruct(const SurfaceStruct &s);//copy constructor
   
   friend void swap(SurfaceStruct &first, SurfaceStruct &second);//swap
   
   SurfaceStruct& operator = (SurfaceStruct rhs);// assign 
   
   // void merge(const SurfaceStruct &first);//merge with another surface
   void read(std::string file);//read surface from file

   ~SurfaceStruct();//destructor

   void flip(int i, int j, int k); // flip at (i,j,k)
   void flip(Sub sub);
   void flip(int index);
   
   bool is_inside(Sub sub);
   bool is_outside(Sub sub);

   bool near_interface(Sub sub);
   vector<Sub> get_nbrs(Sub sub);

};

#endif