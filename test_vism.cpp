#include <iomanip>
#include <gtest/gtest.h>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_sinks.h>
#include <fmt/core.h>


#include "vism_addon.h"
#include "globals.h"

#include "utiltest.h"
#include "kernel.h"

using namespace std;

extern int KERNELRAD;


#include <random>
std::random_device dev;
std::mt19937 rng(dev());


const double sigma = 3.5, epsilon = 0.3, q = 13, tau = 0;


TEST(surf, conv) {

    uniform_real_distribution<double> udist(0.0, 0.1);

    const double kscale = 3;
    const double gridbound = 1;
    const double sphere_rad = 0.5;
    const int numrand = 5;
    double exact_area = 4 * M_PI * sphere_rad * sphere_rad;
    vector<double> gridnums = {20,40,80};

    spdlog::info("{:<10} {:<10} {:<10} {:<10}", "gn", "krad", "mae", "maxe");

    for (int gridnum : gridnums ){
        GridData grid(gridnum,gridbound);
        SurfaceStruct surf(grid);
        double krad = kscale * sqrt(grid.dx[0]);
        KernelStruct kernel(grid, krad);

        vector<double> relerr;

        for( int i = 0; i < numrand; i ++){
            // random center
            Coord random_center = {udist(rng),udist(rng),udist(rng)};
            create_sphere( surf, sphere_rad, random_center);
            double apprx_area = getlengthlist(surf.phi, kernel, grid);
            relerr.push_back(abs(apprx_area - exact_area)/exact_area);
        }

        double sum = std::accumulate(relerr.begin(), relerr.end(), 0.0);
        double mae = sum / relerr.size();
        double maxe = *max_element(relerr.begin(), relerr.end());

        spdlog::info("{:<10} {:<10.2f} {:<10.6f} {:<10.6f} ",gridnum, krad, mae, maxe);
    }

}

// 1 atom
TEST(vism, methane) {

    uniform_real_distribution<double> udist(0.0, 0.1);
    // paramters
    ::gamma0 = 0.174;
    ::WATERSIG = sigma;
    ::WATEREPS = epsilon;
    double gridbound = 5;
    double kscale = 3;
    int type;
    int nrand = 5; // number of random centers
    vector<int> grids = {20,40,80};


    // theoretical results
    double r = VismOptimalRadius(sigma, epsilon, q, gamma0, tau, rho0, epsilonin, epsilonex);
    auto [Gsurf,Gvdw,Gelec,Gtot] = VismSphereEnergy(r, sigma, epsilon, q, gamma0, tau, rho0, epsilonin, epsilonex);
    
    spdlog::info("{:<7} {:<7.6f} {:<7.6f} {:<7.6f} {:<7.6f}", 0, Gsurf, Gvdw, Gelec, Gtot);

    // header
    SoluteStruct sol = CreateMethane(1,q,epsilon,sigma);

    for (int gridnum : grids){
        GridData grid(gridnum,gridbound);
        double ***lj = matrix<double>(grid);
        double ***cfa = matrix<double>(grid);

        double krad = kscale * sqrt(grid.dx[0]);
        KernelStruct kernel(grid, krad);

        // flow to find surface
        SurfaceStruct surf(grid);
        Energy e = VismEnergyInitFlow(sol, lj, cfa, surf, kernel, grid);

        spdlog::info("{:<7} {:<7.6f} {:<7.6f} {:<7.6f} {:<7.6f} ", gridnum, e.surf, e.ljsw(), e.esw(), e.total());
        
        
        
        free_matrix<double>(lj, grid);
        free_matrix<double>(cfa, grid);

    }
    
}


// two atom, look at surface, compare with continuous vism
TEST(vism, 2atom) {

    // paramters
    
    ::gamma0 = 0.174;
    ::WATERSIG = sigma;
    ::WATEREPS = epsilon;

    double kscale = 3;
    vector<double> ds = {4,6,8}; // spacing between two atom
    double gridbound = 10;
    int gridnum = 100;

    // log to console and file
    string logname = fmt::format("bvismtests/atom_conv_q{:.1f}.txt",q);
    auto console_sink = std::make_shared<spdlog::sinks::stdout_sink_mt>();
    console_sink->set_pattern("%v");
    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(logname, true);
    console_sink->set_pattern("%v");
    spdlog::logger logger("multi_sink", {console_sink, file_sink});

    // header
    logger.info("{:>10}, {:>10}, {:>10}, {:>10}, {:>10}, {:>10}, {:>10}, {:>10} {:>10}",
        "d", "surf", "ljswbox", "ljswout", "ljsw", "eswbox", "eswout", "esw", "total");


    SoluteStruct sol = CreateMethane(2, q,epsilon,sigma);
 
    GridData grid(gridnum,gridbound);
    double ***lj = matrix<double>(grid);
    double ***cfa = matrix<double>(grid);
    

    SurfaceStruct surf(grid);
    double krad = kscale * sqrt(grid.dx[0]);
    KernelStruct kernel(grid, krad);


    for (double d: ds ){
        sol.d[0][0] = -d/2;
        sol.d[1][0] = d/2;
        Energy e = VismEnergyInitFlow(sol, lj, cfa, surf, kernel ,grid);
        logger.info("{:10}, {:10.6f}, {:10.6f}, {:10.6f}, {:10.6f}, {:10.6f}, {:10.6f}, {:10.6f} {:10.6f}",
                 d, e.surf, e.ljswbox, e.ljswout, e.ljsw(),
                 e.eswbox, e.eswout, e.esw(), e.total());
        
        string surfname = fmt::format("bvismtests/2atom_surf_d{:.1f}.txt",d);
        write_field(surfname, grid, surf.phi);
    }
    
    free_matrix<double>(lj, grid);
    free_matrix<double>(cfa, grid);
}



int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);

    spdlog::set_level(spdlog::level::info); // Set global log level to debug
    spdlog::set_pattern("%v");

    return RUN_ALL_TESTS();
}

