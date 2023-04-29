/*
Used in binary paper for 1 atom test.
*/
#include <iomanip>
#include <cmath>
#include "vism_addon.h"
#include "vism.h"
#include "fftconv.h"

#include "globals.h"
#include "utiltest.h"

#include <spdlog/spdlog.h>
// #define FMT_HEADER_ONLY
// #include "fmt/core.h"
// #include "fmt/format.h"
// #include "fmt/ostream.h"



using namespace std;

extern int KERNELRAD;


// 1 atom
int main(int argc, char **argv) {
    
    spdlog::set_level(spdlog::level::info); // Set global log level to debug
    spdlog::set_pattern("%v");
    
    CmdLine(argc, argv);
    PrintOptions(cout);

    cxxopts::Options options("vismfft", "vismfft");
    options.allow_unrecognised_options();
    options.add_options()
        ("r0", "init rad", cxxopts::value<double>()->default_value("-1"))
        ("q,charge", "charge",cxxopts::value<double>()->default_value("1.0"))
    ;
    auto result = options.parse(argc, argv);
    double r0 = result["r0"].as<double>();
    double q = result["q"].as<double>();
    
    double sigma = 3.5, epsilon = 0.3; 
    double tau = 0;
    ::WATERSIG = sigma;
    ::WATEREPS = epsilon;
    ::gamma0 = 0.174;

    // SoluteStruct sol = SolFromFile(SOLFILE);
    SoluteStruct sol = CreateMethane(1,q,epsilon,sigma);

    // double d = 6;
    // SoluteStruct sol = CreateMethane(2, q, epsilon, sigma);
    // sol.d[0][0] = -d/2;
    // sol.d[1][0] = d/2; 
        

    sol.Centering();
    
    GridData grid(GRIDNUM,GRIDBOUND);
    sol.LJ = matrix<double>(grid);
    sol.GB = matrix<double>(grid);
        
    
    KernelStruct kernel(grid, KERNELRAD);

    // surface
    Coord center = {0,0,0};
    SurfaceStruct surf(grid);
    if (r0 > 0){
        create_sphere(surf, r0, center);
        spdlog::info("create sphere or radius {}\n",r0);
    } else {
        getinitsurf(surf, sol, grid, TIGHTFIT);
        spdlog::info("initsurf tight = {}\n", TIGHTFIT);
    }
        
    
    cal_LJ(sol, sol.LJ, grid);
    seqCFA(sol,sol.GB,grid);	
    
    // fft 
    FFTWStruct fftw;
    getinitfftw(fftw, kernel, grid);

    Energy e = VismEnergyOnly(sol, sol.LJ, sol.GB, surf, kernel, grid);
    spdlog::info( "initial energy surf={}, ljswbox={}, eswbox = {} total={}\n", e.surf, e.ljswbox, e.eswbox,  e.total());
    
    // Let's flow with fft
    binaryflowfftw(surf, sol, kernel, fftw, grid);
    e = VismEnergyOnly(sol, sol.LJ, sol.GB, surf, kernel, grid);
    spdlog::info( "fft energy surf={}, ljswbox={}, eswbox = {} total={}\n", e.surf, e.ljswbox, e.eswbox,  e.total());

    // reset surface
    if (r0 > 0){
        create_sphere(surf, r0, center);
    } else {
        getinitsurf(surf, sol, grid, TIGHTFIT);
    }
    FlipInfo info = binaryflowheaplistinterfaceonly(surf, sol, sol.LJ, sol.GB, kernel, grid);
    
    e = VismEnergyOnly(sol, sol.LJ, sol.GB, surf, kernel, grid);
    spdlog::info( "flip energy surf={}, ljswbox={}, eswbox = {} total={}\n", e.surf, e.ljswbox, e.eswbox,  e.total());

    free_matrix<double>(sol.LJ, grid);
    free_matrix<double>(sol.GB, grid);

}
