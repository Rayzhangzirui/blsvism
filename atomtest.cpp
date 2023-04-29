/*
Used in binary paper for 1 atom test.
*/
#include <iomanip>
#include <cmath>
#include "vism_addon.h"
#include "vism.h"

#include "globals.h"
#include "utiltest.h"

#define FMT_HEADER_ONLY
#include "fmt/core.h"
#include "fmt/format.h"
#include "fmt/ostream.h"

using namespace std;

extern int KERNELRAD;

#include <filesystem>
namespace fs = std::filesystem;

#include <random>
std::random_device dev;
std::mt19937 rng(dev());


// compare two sphere
pair<double,double> CompareSphere(SurfaceStruct& surf, Coord center){
    double min_rad_out = surf.grid.b[0]; // min radius in outside region
    double max_rad_in = 0; // max radius in insdie region

    #pragma omp parallel for collapse(3)
    for(int i = 0; i <= surf.grid.nx[0]; i++){
        for(int j = 0; j <= surf.grid.nx[1]; j++){
            for(int k = 0; k <= surf.grid.nx[2]; k++){
                Sub sub = {i,j,k};
                Coord x = sub2coord( sub, surf.grid); // coord of grid
                double rad = dist2(x,center);
                if(surf.is_inside(sub) && rad > max_rad_in){
                    max_rad_in = rad;
                }
                if(surf.is_outside(sub) && rad < min_rad_out){
                    min_rad_out = rad;
                }
            }
        }
    }
    return {max_rad_in,min_rad_out};
}

// output difference of two spherical surface with, as = approx surface, es = exact surface
// if approx surface has different sign, output distance to center
void WriteSphereDiff(SurfaceStruct& as, SurfaceStruct& es, Coord center, double r, string filename){
    ofstream fo(filename);
    
    fo<<setprecision(9);
    fo<<center[0]<<","<<center[1]<<","<<center[2]<<","<<r<<"\n";
    for(int i = 0; i <= as.grid.nx[0]; i++){
        for(int j = 0; j <= as.grid.nx[1]; j++){
            for(int k = 0; k <= as.grid.nx[2]; k++){
                if(as.phi[i][j][k]*es.phi[i][j][k]<0){
                    Coord x = sub2coord(array<int,3>{i,j,k}, as.grid);
                    fo<<i<<","<<j<<","<<k<<","<<dist2(x,center)-r<<"\n";
                }
            }   
        }   
    }
    
}


// 1 atom
int main(int argc, char **argv) {

    CmdLine(argc, argv);
    PrintOptions(cout);

    cxxopts::Options options("convergce", "one atom convergence test");
    options.allow_unrecognised_options();
    
    options.add_options()
        ("a,start", "start grid size", cxxopts::value<int>()->default_value("20"))
        ("h,inc", "increment", cxxopts::value<int>()->default_value("20"))
        // ratio of geometri sequence. if 1, then arithmetic sequence
        ("o,ratio", "ratio", cxxopts::value<int>()->default_value("1"))
        ("n,num", "number of increament", cxxopts::value<int>()->default_value("2"))
        ("r,nrand", "number of rand center", cxxopts::value<int>()->default_value("0"))
        ("k,kscale", "kscale",cxxopts::value<double>()->default_value("3"))
        ("m,dname", "dir name",cxxopts::value<string>()->default_value("test"))
        ("w,wsurf", "write surface",cxxopts::value<bool>()->default_value("false"))
        ("q,charge", "charge",cxxopts::value<double>()->default_value("1.0"))
    ;

    auto result = options.parse(argc, argv);

    int a = result["start"].as<int>();
    int h = result["inc"].as<int>(); 
    int o = result["ratio"].as<int>(); 
    int n = result["num"].as<int>();
    int nrand = result["nrand"].as<int>();
    bool ws = result["wsurf"].as<bool>();
    string dname = result["dname"].as<string>();
    double kscale = result["kscale"].as<double>();
    double q = result["charge"].as<double>();


    // provide sub directory name for the test
    fs::path p("./bvismtests");
    p /= dname;
    if (fs::exists(p)){
        cout<< "save to "<<p<<endl;
    }else{
        cout<<"dir does not exist"<<endl;
        exit(0);
    }

    uniform_real_distribution<double> udist(0.0, 0.1);

    // paramters
    double sigma = 3.5, epsilon = 0.3; 
    double tau = 0;
    ::WATERSIG = sigma;
    ::WATEREPS = epsilon;
    int type;

    // Sequence of grid
    vector<int> grids;
    if(o==1){
        // Arithmetic sequence
        grids = LinearSpacedArray(a,h,n);    
    } else {
        // geometric sequence
        grids =  GeometricArray(a,o,n);    
    }
    
    

    // gridnums and random centers
    vector<Coord> cts;
    cts.push_back({0,0,0});
    for (int i = 0; i < nrand ; i ++){
        cts.push_back({udist(rng),udist(rng),udist(rng)});
    }

    // log to console and file
    string filename = fmt::format("./bvismtests/{}/atom_conv_q{:1}_a{}_h{}_o{}_n{}_nr{}_g{}_k{}_gb{}.txt",dname,q, a, h, o, n, nrand, gamma0, kscale,GRIDBOUND);
    ofstream output(filename);

    // header
    string header_tpl =  "{:>5}, {:>5}, {:>5}, {:>12}, {:>12}, {:>12}, {:>12}, {:>12}, {:>12}, {:>12}, {:>12}, {:>12}, {:>12} {:>12} {:>12}\n";
    fmt::print(output, header_tpl,"type", "gn", "n", "surf", "ljswbox", "ljswout", "ljsw", "eswbox", "eswout", "esw", "total", "rmaxin", "rminout","flip","time");
    
    // row template
    string row_tpl =  "{:>5}, {:>5}, {:>5}, {:>.9f}, {:>.9f}, {:>.9f}, {:>.9f}, {:>.9f}, {:>.9f}, {:>.9f}, {:>.9f}, {:>.9f}, {:>.12f}, {:}, {:>.6f}\n";

    // theoretical results
    double r = VismOptimalRadius(sigma, epsilon, q, gamma0, tau, rho0, epsilonin, epsilonex);
    auto [Gsurf,Gvdw,Gelec,Gtot] = VismSphereEnergy(r, sigma, epsilon, q, gamma0, tau, rho0, epsilonin, epsilonex);
    fmt::print(output, row_tpl, -1, 0, 0, Gsurf, NAN, NAN,  Gvdw, NAN, NAN, Gelec, Gtot, r, r, 0, 0.0);

    SoluteStruct sol = CreateMethane(1,q,epsilon,sigma);

    for (int gridnum : grids){
        cout<<"gridnum "<<gridnum<<endl;

        GridData grid(gridnum,GRIDBOUND);
        double ***lj = matrix<double>(grid);
        double ***cfa = matrix<double>(grid);
            

        
        double krad = kscale * sqrt(grid.dx[0]);
        KernelStruct kernel(grid, krad);


        for (int i = 0; i < cts.size(); i++){
            cout<<"i "<<i<<endl;
            sol.d[0][0] = cts[i][0];
            sol.d[0][1] = cts[i][1];
            sol.d[0][2] = cts[i][2];

            // flow to find surface
            SurfaceStruct surf(grid);
            
            type = 0;
            auto [e,info] = VismEnergyInitFlow(sol, lj, cfa, surf, kernel, grid);
            auto [maxin,minout] = CompareSphere(surf, cts[i]);
            fmt::print(output,row_tpl,
             type, gridnum, i,  e.surf, e.ljswbox, e.ljswout, e.ljsw(),
             e.eswbox, e.eswout, e.esw(), e.total(), maxin, minout, info.numflip, info.time);    


            // use optimal radius as exact surface (esurf)
            SurfaceStruct esurf(grid);
            type = 1;
            create_sphere(esurf, r, cts[i]);
            Energy e2 = VismEnergyOnly(sol, lj, cfa, esurf, kernel, grid); // energy using exact surface
            auto [maxin2,minout2] = CompareSphere(esurf, cts[i]); // compare flow surf and esurf
            
            fmt::print(output,row_tpl,
             type, gridnum, i,  e2.surf, e2.ljswbox, e2.ljswout, e2.ljsw(),
             e2.eswbox, e2.eswout, e2.esw(), e2.total(), maxin2, minout2, 0, 0.0);

            // output detailed difference of two surface
            string surfname = fmt::format("./bvismtests/{}/gn{}_i{}",dname,gridnum,i);
            WriteSphereDiff(surf, esurf, cts[i], r,  surfname);

            if(ws){
                // write full surface
                string fullsurf = fmt::format("./bvismtests/{}/fullsuf_gn{}_i{}",dname,gridnum,i);
                write_field(fullsurf, surf.grid, surf.phi);    
            }
            
        }
        
        
        free_matrix<double>(lj, grid);
        free_matrix<double>(cfa, grid);

    }

}



