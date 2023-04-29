/*
Used in binary paper for 2 atom test. Output surface
*/

#include <iomanip>
#include <cmath>
#include "vism_addon.h"
#include "vism.h"

#include "globals.h"
#include "utiltest.h"

#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_sinks.h>
#include <fmt/core.h>

using namespace std;

extern int KERNELRAD;

#include <filesystem>
namespace fs = std::filesystem;



int main(int argc, char **argv) {

    CmdLine(argc, argv);
    PrintOptions(cout);
    
    // extra options
    cxxopts::Options options("2atom", "2atom");
    options.allow_unrecognised_options();
    options.add_options()
        ("c,code", "which test", cxxopts::value<int>()->default_value("0"))
    ;
    // code = 0, two atoms are at (d/2,0,0) and (-d/2,0,0)
    // code = 1, two atoms are at (+/-) d (2/sqrt(3)) (1,1,1), at diagonal
    // used to investigate grid effect

    auto result = options.parse(argc, argv);

    // code = 0: two atom at x-axis; code = 1: along (1,1,1)
    int code = result["code"].as<int>();
    
    const double sigma = 3.5, epsilon = 0.3;

    ::gamma0 = 0.174;
    ::WATERSIG = sigma;
    ::WATEREPS = epsilon;
    double q = 1;

    string FITSTR = (TIGHTFIT)?"tight":"loose";

    double kscale = 3;
    vector<double> ds = {4,6,8}; // spacing between two atom
    double gridbound = 10; //default 10
    int gridnum = 100; // default 100

    // log to console and file
    string logname = fmt::format("bvismtests/2atom/2atom_q{:1}_v{:1}_{}.txt",q,code,FITSTR);
    auto console_sink = std::make_shared<spdlog::sinks::stdout_sink_mt>();
    console_sink->set_pattern("%v");
    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(logname, true);
    file_sink->set_pattern("%v");
    spdlog::logger logger("multi_sink", {console_sink, file_sink});

    // header
    logger.info("{:>10}, {:>10}, {:>10}, {:>10}, {:>10}, {:>10}, {:>10}, {:>10}, {:>10}, {:>10}, {:>10}",
                "d", "surf", "ljswbox", "ljswout", "ljsw", "eswbox", "eswout", "esw", "total","flip","time");


    SoluteStruct sol = CreateMethane(2, q, epsilon, sigma);
 
    GridData grid(gridnum,gridbound);
    double ***lj = matrix<double>(grid);
    double ***cfa = matrix<double>(grid);
    

    SurfaceStruct surf(grid);

    double krad = kscale * sqrt(grid.dx[0]);
    KernelStruct kernel(grid, krad);


    for (double d: ds ){
        if (code == 0){
            sol.d[0][0] = -d/2;
            sol.d[1][0] = d/2;    
        } else {
            for(int i=0; i<3; i++){
                sol.d[0][i] = -d/2/sqrt(3);
                sol.d[1][i] = d/2/sqrt(3);                    
            }
        }
        
        auto [e,info]= VismEnergyInitFlow(sol, lj, cfa, surf, kernel ,grid, TIGHTFIT);
        logger.info("{:10}, {:10.6f}, {:10.6f}, {:10.6f}, {:10.6f}, {:10.6f}, {:10.6f}, {:10.6f}, {:10.6f}, {:10}, {:10.6f}",
                 d, e.surf, e.ljswbox, e.ljswout, e.ljsw(),e.eswbox, e.eswout, e.esw(), e.total(),info.numflip, info.time);
        
        string surfname = fmt::format("bvismtests/2atom/2atom_surf_q{:1}_d{:1}_v{:1}.txt",q,d,code);
        write_field(surfname, grid, surf.phi);
    }
    
    free_matrix<double>(lj, grid);
    free_matrix<double>(cfa, grid);

}



