#include "util.h"
#include "globals.h"
#include <sstream>


// Geometry
int GRIDNUM;
double GRIDBOUND;

// physics constants
double gamma0;
double T;

int HALFLJ;
int LJSKIP;

double WATERSIG = 3.166; // AA
double WATEREPS = 0.6502; // kJ/mol

// WATERSIG = 3.1536; //AA used in MC-VISm
// WATEREPS = 0.2601; //kbt

// VISM parameters
bool TIGHTFIT;
int KERNELRAD;
bool CFA;

// IO
string SOLFILE;
string OUTPUTFILE; // output

// MC parameters
int SEED;
int MAXSTEP;
double MAXTRANS;
int METHOD;
bool WATEROFF; // 

int INFOSTEP;
int SOLMCSTEP; // # of swipe of all atoms
int VISMMCSTEP; // # of random perturbation of surface
int WRITESTEP;

// surface tension
// double gamma0 = 0.43;// kj/mol/A^2, from MS thesis fig 7
// double gamma0 = 0.33;// kj/mol/A^2, from MS thesis fig 9
// double gamma0 = 2.479 * 0.1315; // = 0.3259885 PB paper, Zhou, gamma = 0.1315 kbT/A^2, 1 kbT = 2.479 kJ/mol

cxxopts::Options options("vism", "vism calculation");

void CmdLine(int argc, char *argv[]){
    options.allow_unrecognised_options();    

    options.add_options()
    	// Geometry
        ("gridbound", "boundary", cxxopts::value<double>()->default_value("30.0"))
        ("gridnum", "grid number", cxxopts::value<int>()->default_value("60"))

        // constants
        ("gamma0", "gamma",cxxopts::value<double>()->default_value("0.3259885"))
        ("temperature", "temperature",cxxopts::value<double>()->default_value("300"))
        ("halflj", "halflj",cxxopts::value<int>()->default_value("1"))
        ("ljskip", "ljskip",cxxopts::value<int>()->default_value("2"))

        //VISM
        ("tightfit", "tightfit",cxxopts::value<bool>()->default_value("true"))
        ("kernelrad", "kernalrad",cxxopts::value<int>()->default_value("3"))
        ("cfa", "cfa",cxxopts::value<bool>()->default_value("true"))

        // IO
        ("outputfile", "output file",cxxopts::value<string>()->default_value(""))
        ("solfile", "solute file",cxxopts::value<string>()->default_value(""))

        // MC
        ("seed", "seed", cxxopts::value<int>()->default_value("0"))
        ("maxstep", "mc step", cxxopts::value<int>()->default_value("10"))
        ("maxtrans", "max translation", cxxopts::value<double>()->default_value("0.05"))
        ("method", "mc method", cxxopts::value<int>()->default_value("1"))
        ("infostep", "info step", cxxopts::value<int>()->default_value("10"))
        ("wateroff", "water off", cxxopts::value<bool>()->default_value("false"))

        ("solmcstep", "sol mc step", cxxopts::value<int>()->default_value("1"))
        ("vismmcstep", "vism mc step", cxxopts::value<int>()->default_value("10"))
        ("writestep", "write step", cxxopts::value<int>()->default_value("-1"))
    ;

    auto result = options.parse(argc, argv);

    // Geometry
    GRIDBOUND = result["gridbound"].as<double>();
    GRIDNUM = result["gridnum"].as<int>();

    // constants
    gamma0 = result["gamma0"].as<double>();
    T = result["temperature"].as<double>();
    
    HALFLJ = result["halflj"].as<int>();
    LJSKIP = result["ljskip"].as<int>();

    //VISM
	KERNELRAD = result["kernelrad"].as<int>();
	TIGHTFIT = result["tightfit"].as<bool>();
    CFA = result["cfa"].as<bool>();

    //IO
    OUTPUTFILE = result["outputfile"].as<string>();
    SOLFILE = result["solfile"].as<string>();

    // MC
    SEED = result["seed"].as<int>();
    MAXSTEP = result["maxstep"].as<int>();
    MAXTRANS = result["maxtrans"].as<double>();
    METHOD = result["method"].as<int>();
    INFOSTEP = result["infostep"].as<int>();
    WATEROFF = result["wateroff"].as<bool>();

    SOLMCSTEP = result["solmcstep"].as<int>();
    VISMMCSTEP = result["vismmcstep"].as<int>();
    WRITESTEP = result["writestep"].as<int>();

}


void PrintOptions(ostream& out){
    // geometry
	out<<"# gridbound " 		<<GRIDBOUND		<<endl;
	out<<"# gridnum "	 		<<GRIDNUM		<<endl;
    // physics
	out<<"# gamma0 "	 		<<gamma0		<<endl;
	out<<"# temperature "	 	<<T		<<endl;
    out<<"# halflj "            <<HALFLJ     <<endl;
    out<<"# ljskip "            <<LJSKIP     <<endl;
	// vism
	out<<"# tightfit "			<<TIGHTFIT		<<endl;
	out<<"# kernelrad " 		<<KERNELRAD		<<endl;
    out<<"# cfa "         <<CFA     <<endl;
    
    //io
	out<<"# outputfile " 		<<OUTPUTFILE	<<endl;
    out<<"# solfile "           <<SOLFILE       <<endl;
    //mc
	out<<"# seed "				<<SEED			<<endl;
	out<<"# maxstep "			<<MAXSTEP		<<endl;
	out<<"# maxtrans "			<<MAXTRANS		<<endl;
	out<<"# method "			<<METHOD		<<endl;
	out<<"# info_step "			<<INFOSTEP		<<endl;
    out<<"# wateroff "          <<WATEROFF      <<endl;

    out<<"# solmcstep "            <<SOLMCSTEP        <<endl;
    out<<"# vismmcstep "           <<VISMMCSTEP        <<endl;
    out<<"# writestep "            <<WRITESTEP        <<endl;
	
}