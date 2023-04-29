// test convergence of area using kernel2.cpp
#include <string>
#include <iomanip>
#include "cxxopts.hpp"
#include "surf.h"
#include "kernel2.h"
#include "Vector.h"
#include "utiltest.h"

#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_sinks.h>
#include <fmt/core.h>

#include <memory>

using namespace std;

#include <random>
std::random_device dev;
std::mt19937 rng(dev());

using Vector3d = va::Vector;

// inside false, outside true
bool sphere(Vector3d x, Vector3d o, double r){
    return (norm(x-o)>r);
}

// sign distance function of sphere
double SignDistSphere(Vector3d x, Vector3d o, double r){
    return (norm(x-o)-r);
}

// quadratic distance function of sphere
double QuadSphere(Vector3d x, Vector3d o, double r){
    return pow(norm(x-o)-r,2);
}

// ellipsoid with center o, side a
double LsEllipsoid(Vector3d p, Vector3d o, Vector3d a){
    return (pow( (p.x()-o.x()) / a.x(),2) 
          + pow( (p.y()-o.y()) / a.y(),2) 
          + pow( (p.z()-o.z()) / a.z(),2)) - 1;
}

double LsTorus(Vector3d p, Vector3d o, double R, double r){
    return pow(sqrt(pow(p.x()-o.x(),2) + pow(p.y()-o.y(),2))-R,2) + pow(p.z()-o.z(),2) - pow(r,2);
}



int main(int argc, char *argv[]){

    // Parse arguments

    cxxopts::Options options("vism", "vism calculation");

    options.add_options()
        // Geometry
        ("b,bound", "gridbound", cxxopts::value<double>()->default_value("1"))
        ("R,radius", "radius of sphere/torus", cxxopts::value<double>()->default_value("0.5"))
        ("a,start", "start grid size", cxxopts::value<int>()->default_value("20"))
        ("h,inc", "increment", cxxopts::value<int>()->default_value("20"))
        ("o,ratio", "ratio", cxxopts::value<int>()->default_value("1"))
        ("n,num", "number of increament", cxxopts::value<int>()->default_value("2"))
        ("r,nrand", "number of rand center", cxxopts::value<int>()->default_value("0"))
        
        ("surfcode", "surf code ", cxxopts::value<int>()->default_value("0"))
        ("fcncode", "fcn code ", cxxopts::value<int>()->default_value("0"))

        // kernel radius  = scale * h
        ("kscale", "kernalrad",cxxopts::value<double>()->default_value("1.0"))
        ("kcode", "code ", cxxopts::value<int>()->default_value("0"))

        // IO
        ("w,write", "write output",cxxopts::value<bool>()->default_value("true"))

        ("p,help", "Print usage")

    ;

    auto result = options.parse(argc, argv);

    // Geometry
    double gridbound = result["bound"].as<double>();    
    double R = result["radius"].as<double>();    
    int a = result["start"].as<int>();
    int h = result["inc"].as<int>();
    int o = result["ratio"].as<int>(); // ratio of geometri sequence
    int n = result["num"].as<int>();
    int nrand = result["nrand"].as<int>();
    
    int surfcode = result["surfcode"].as<int>();
    int fcncode = result["fcncode"].as<int>();

    // kernel
    double kscale = result["kscale"].as<double>();    
    int kcode = result["kcode"].as<int>();
    
    //IO
    bool wo = result["write"].as<bool>();

    if (result.count("help"))
    {
      cout << options.help() << endl;
      exit(0);
    }

    // Setup

    uniform_real_distribution<double> udist(0.0, 0.1);
    
    vector<int> grids;
    if(o==1){
        // Arithmetic sequence
        grids = LinearSpacedArray(a,h,n);    
    } else {
        // geometric sequence
        grids =  GeometricArray(a,o,n);    
    }
    
    double eA,eI,aA,aI; // exact/approx area/integral
    function<double(Vector3d)> ls; // level set
    function<bool(Vector3d)> bls; // bianry level set
    function<double(Vector3d)> fcn; // smooth function
    
    vector<Vector3d> cts;
    cts.push_back(Vector3d(0,0,0));
    for (int i = 0; i < nrand ; i ++){
        cts.push_back(Vector3d(udist(rng),udist(rng),udist(rng)));
    }


    // Setup logger

    string logname = fmt::format("bvismtests/area/area_surf{:1}_kscale{:1}_kcode{:1}_fcn{:1}_a{}_h{}_n{}_nr{}_b{}_R{}.txt",surfcode, kscale, kcode,fcncode, a, h, n, nrand, gridbound, R);
    auto console_sink = std::make_shared<spdlog::sinks::stdout_sink_mt>();
    console_sink->set_pattern("%v");

    std::vector<std::shared_ptr<spdlog::sinks::sink>> sink_vec = {console_sink};

    if(wo){
        cout<<"write output"<<endl;
        auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(logname, true);
        file_sink->set_pattern("%v");
        sink_vec.push_back(file_sink);
    }
    
    spdlog::logger logger("multi_sink", sink_vec.begin(), sink_vec.end());
    // header
    logger.info("{:>10}, {:>10}, {:>10}","trial", "gn", "err","relerr");
    
    for (int gn : grids){

        GridData grid(gn, gridbound);
        double rkrad = kscale * sqrt(grid.dx[0]); // real-valued kernel radius
        Kernel kernel;          
        InitKernelCode(kernel, grid, rkrad, kcode);

        for (int i = 0; i < cts.size(); i ++){
            Vector3d center = cts[i];
            // set up exact solutions

            if(surfcode == 0){
                double sr = R; // radius of sphere
                eA = 4.0 * M_PI * sr * sr; // exact area
                ls = [center,sr](Vector3d x) {return SignDistSphere(x, center, sr);};

            } else if (surfcode == 1){
                eA = 3.125273723338068; // surface area of ellipsoid with axis 0.4 0.5 0.6
                Vector3d axis(0.4,0.5,0.6);
                ls = [center,axis](Vector3d x) {return LsEllipsoid(x, center, axis);};
            } else if ( surfcode == 2){
                double r = 0.25;
                eA = 4*M_PI*M_PI*R*r;
                ls = [center,R,r](Vector3d x) {return LsTorus(x, center, R, r);};
            }
            else{
                cerr<<"underfined surface"<<endl;
                exit(1);
            }

            bls = [&](Vector3d x){return ls(x)>0;};// binary level set function

            // Set up function, smoothly extended off interface
            if(fcncode == 0){
                // integrate the constant 1 function, should be area
                eI = eA;
                fcn = [&](Vector3d x) {return 1;};
            } else if (fcncode == 1) {
                // integrate the zero level set function, should be 0
                eI = 0;
                fcn = ls;
            } else{
                cerr<<"undefined function"<<endl;
                exit(1);
            }
        
            
            // aI = ComputeSurfaceIntegral(kernel, grid, bls, fcn);// approximate surface integral
            aI = ComputeAreaFcn(kernel, grid, bls);// approximate area        
            // cout<<setprecision(16)<<gn<<","<<aA - eA<<","<<aI - eI<<endl;
            double err = aI - eI;
            double relerr = err/eI;
            logger.info("{:10}, {:10}, {:10}, {:10}", i, gn, err, relerr);
        }
        logger.flush();


    }
    

    return 0;
}