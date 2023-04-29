#ifndef SOLUTE_H
#define SOLUTE_H

#include "Vector.h"
#include "Rotation.h"
#include <vector>
#include <string>

extern double c25sigma;
extern double c25eps;

struct SoluteStruct
{
   int dnum;//number of atom
   int *index;//index from file
   double **d;//matrix size dnum x 3
   double *epsilon;//LJ epsilon
   double *sigma;//LJ sigma
   double *Q;//charge
   double rho0 = 0.0333;//solvent number density A^-3
   int estat;//electrostatic state: 0 without cfa, 1 with cfa
   double epsilonex = 80.0;// exterior (solvent) relative permittivity 80
   double epsilonin = 1.0;// interior (vacumm) relative permittivity 1
   double epsilon0 = 0.00014319;// 0.00014319
   
   // constant from C25 shenggao
   double sigLJss=3.73;//LJss sigma
   double epsLJss=0.5856;//LJss eps: kj/mol
   double AA=3347.200;//kJ/mol/AA^2
   double r0=1.53;//AA
   double thet0=111.*M_PI/180.;//radian, 111 degree
   double BB=462.000;//kJ/mol/rad^2 
   // the unit is  462 kJ/mol/deg^2 in paper
   // but from this source, https://aip.scitation.org/doi/10.1063/5.0015184. 462 kJ/mol/rad^2 seems to be more reasonable

   // for compatibility with LT's code
   double ***LJ;
   double ***GB;

   // empty constructor
   SoluteStruct();
   //constructor
   SoluteStruct(int n);

   //copy constructor
   SoluteStruct(const SoluteStruct& s);

   //destructor
   ~SoluteStruct();

   // read data, starting and ending index
   void readdata(std::string filename, int start, int end);

   // read data, only read position
   void readpos(std::string filename);
   
   // get centroid, arithmetic average of x y z coordinates
   va::Vector centroid() const;
   
   // get radius of gyration
   double RadiusGyration() const;

   // Print first n atom
   void Print(int n) const;

   // write
   void Write(std::string fname) const;

   // roate solute
   void rotate(va::Rotation R);
   void rotate(va::Vector u, double th);
   void rotate(double ax, double ay, double az, double th);

   // translate solute
   void translate(va::Vector v) ;

   // center the solute
   void Centering();

   //swap
   friend void swap(SoluteStruct &first, SoluteStruct &second);
   
   // assign 
   SoluteStruct& operator= (SoluteStruct rhs);
   
   // pairwise distance
   double min_pairwise_dist(const SoluteStruct& s2) const;
   
   // center to center distance
   double c2c_dist(const SoluteStruct& s2) const;
   
   //get center of two solute (average of min max)
   static va::Vector get_center(std::initializer_list<SoluteStruct> solutes);

   double MaxCoord() const;
   
};


SoluteStruct CreateMethane(int n = 1, double charge = 0.0, double eps = c25eps, double sigma = c25sigma);

SoluteStruct CreateC25(int n);

SoluteStruct SolFromFile(std::string filename);

#endif