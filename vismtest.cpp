/*
Use in VISM paper for biomolecule test.
*/
#include "kernel.h"
#include "surf.h"
#include "vism.h"

#include "vism_addon.h"
#include "globals.h"
#include "utiltest.h"

#include <fstream>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_sinks.h>

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>

using namespace std;

#include <filesystem>
namespace fs = std::filesystem;

extern cxxopts::Options options;

bool WSURF;
bool CENTERING;

string FITSTR;

// file = a file containing a list of positions
vector<string> parse_files(string file){
    ifstream infile(file);
    vector<string> vf;
    
    string line;
    while(getline(infile,line)){
        vf.push_back(line);
    }
    return vf;
}


void ComputeVismEnergy(string solfile, string surfile){

    SoluteStruct sol = SolFromFile(solfile);

    if(CENTERING){
        sol.Centering();
    }

    GridData grid(GRIDNUM,GRIDBOUND);

    double ***lj = matrix<double>(grid);
    double ***cfa = CFA? matrix<double>(grid): nullptr;
    
    SurfaceStruct surf(grid);

    KernelStruct kernel(grid, KERNELRAD);

    Energy e;
    FlipInfo info;

    spdlog::info(
    "{}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n",
     "surf","ljswbox","ljswout","ljsw","eswbox","eswout","esw","total","flip","time");

    if(surfile.empty()){
        auto [e1,info1] = VismEnergyInitFlow(sol, lj, cfa, surf, kernel, grid, TIGHTFIT);
        e = e1;
        info = info1;
    } else {
        surf.read(surfile);
        cal_LJ(sol, lj, grid);
        if(cfa){
            seqCFA(sol,cfa,grid);   
        }
        e = VismEnergyOnly(sol, lj, cfa, surf, kernel, grid);
        spdlog::info(
            "{:>.9f}, {:>.9f}, {:>.9f}, {:>.9f}, {:>.9f}, {:>.9f}, {:>.9f}, {:>.9f}, {:}, {:>.6f}",
            e.surf, e.ljswbox, e.ljswout, e.ljsw(), e.eswbox, e.eswout, e.esw(), e.total(), 0, 0.0);

        info = binaryflowheaplistinterfaceonly(surf, sol, lj, cfa, kernel, grid);
        e = VismEnergyOnly(sol, lj, cfa, surf, kernel, grid);
    }

    

    if(WSURF){
        string surfname = fmt::format("{}_{}_gn{}_k{}_surf",SOLFILE,FITSTR,GRIDNUM, KERNELRAD);
        write_field(surfname, grid, surf.phi);
    }
    
    // cout<<sol.RadiusGyration()<<" ";


    spdlog::info(
        "{:>.9f}, {:>.9f}, {:>.9f}, {:>.9f}, {:>.9f}, {:>.9f}, {:>.9f}, {:>.9f}, {:}, {:>.6f}",
         e.surf, e.ljswbox, e.ljswout, e.ljsw(), e.eswbox, e.eswout, e.esw(), e.total(), info.numflip, info.time);

    free_matrix<double>(lj, grid);
    free_matrix<double>(cfa, grid);
}



int main(int argc, char **argv) {

    options.add_options()
        ("center", "shift to center", cxxopts::value<bool>()->default_value("true"))
        ("insurf", "surface files", cxxopts::value<string>()->default_value(""))
        ("weng", "energy files", cxxopts::value<bool>()->default_value("false"))
        ("wsurf", "write surface", cxxopts::value<bool>()->default_value("false"));


    CmdLine(argc, argv);

    auto result = options.parse(argc, argv);
    string insurf = result["insurf"].as<string>();
    WSURF = result["wsurf"].as<bool>();
    bool weng = result["weng"].as<bool>();
    CENTERING = result["center"].as<bool>();

    cout<<"use water epsilon unit kbt"<<endl;
    WATERSIG = 3.1536;
    WATEREPS = 0.2601;
    gamma0 = 0.174;

    FITSTR = (TIGHTFIT)?"tight":"loose";

    PrintOptions(cout);

    // setup logging
    auto console_sink = std::make_shared<spdlog::sinks::stdout_sink_mt>();
    console_sink->set_pattern("%v");

    std::vector<std::shared_ptr<spdlog::sinks::sink>> sink_vec = {console_sink};

    if(weng){
        string filename;
        if(insurf.empty()){

            filename = fmt::format("{}_{}_gn{}_k{}_eng",SOLFILE,FITSTR,GRIDNUM, KERNELRAD);
        } else{
            filename = fmt::format("{}_k{}_eng",insurf, KERNELRAD);
        }
        auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(filename, true);
        file_sink->set_pattern("%v");
        sink_vec.push_back(file_sink);
    }
    auto combined_logger = std::make_shared<spdlog::logger>("multi_sink", sink_vec.begin(), sink_vec.end());
    spdlog::set_default_logger(combined_logger);
    // end of set up logger
    
    ComputeVismEnergy(SOLFILE,insurf);
    return 0;
}



