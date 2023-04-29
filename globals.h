#ifndef GLOBAL_H
#define GLOBAL_H

#include <string>
#include "cxxopts.hpp"

// tolerance for small number
const double TOL = 1.0e-15;

// Geometry

extern int GRIDNUM;
extern double GRIDBOUND;

// phyiscs constants

extern double gamma0; // surface tension
extern double T; // temperature


// Zhou et al PB VISM, test and applicatoin
const double rho0 = 0.0333;//solvent number density A^-3
const double epsilonex = 80.0;// exterior (solvent) relative permittivity 80
const double epsilonin = 1.0;// interior (vacumm) relative permittivity 1
const double epsilon0 = 0.00014319;// 0.00014319

// water lj spc/e
// https://www.nist.gov/mml/csd/chemical-informatics-research-group/spce-water-reference-calculations-10a-cutoff
// Energy Units Converter
// https://www.weizmann.ac.il/oc/martin/tools/hartree.html
// averaging is performed in energy calculation with water. LJ between solute should not be averaged

// paper: Effects of lengthscales and attractions on the collapse of hydrophobic polymers in water, use epsilon_ww = 0.6502kJ/mol
// simga is not mentioned
extern double WATERSIG; // AA
extern double WATEREPS; // kJ/mol

/* 
in shenggao's code, no LJ interaction with 1st and 2nd neibhbor,
the interaction with 3rd neighbor is 1/2 the usual LJ
*/

// if the interaction with the LJSKIP+1 neighor is halved.
extern int HALFLJ; 

// how many atom to skip. e.g. LJSKIP = 1, do not consider direct neighbor
extern int SKIPLJ; 

// VISM paramter
extern bool TIGHTFIT;
extern int KERNELRAD; // kernel radius, multiples of h
extern bool CFA;

// IO
extern std::string SOLFILE; // file for solute
extern std::string OUTPUTFILE; // output file

// MC parameter
extern int SEED; // random engin seed
extern int MAXSTEP; // max step of 
extern double MAXTRANS;
extern int INFOSTEP; // interval to print infomation
extern int METHOD; // method for MC
// method 0, crude mc, perturb atom, find vism energy, accept or reject
// method 1, 2phase fluctuatoin of surface
// method 2, 2phase flow the surface

extern bool WATEROFF; // do not include water
extern int SOLMCSTEP; // # of swipe of all atoms
extern int VISMMCSTEP; // # of random perturbation of surface
extern int WRITESTEP; // frequency to write position and surface

extern cxxopts::Options options;

void CmdLine(int argc, char *argv[]);
void PrintOptions(std::ostream& out);

#endif