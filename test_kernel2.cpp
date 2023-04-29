#include "kernel2.h"
#include "surf.h"
#include "globals.h"
#include "utiltest.h"

#include <gtest/gtest.h>
#include <numeric>

#include <fmt/os.h>
#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>

#include <random>
std::random_device dev;
std::mt19937 rng(dev());

using namespace std;

TEST(kernel,a0){
	int n = 100;
	double indicator_apprx = ComputeA0Midpoint(n, KERNEL_IND);
	cout<< "n = "<< n<< " indicator_apprx = " << indicator_apprx << " error = " << indicator_apprx - A0_IND<<endl;

	
	double sin2_apprx = ComputeA0Midpoint(n, KERNEL_SIN2 );
	cout<< "n = "<< n<< " sin2_apprx = " << sin2_apprx << " error = " << sin2_apprx - A0_SIN2<<endl;

	double cos_apprx = ComputeA0Midpoint(n, KERNEL_COS );
	cout<< "n = "<< n<< " cos_apprx = " << cos_apprx << " error = " << cos_apprx - A0_COS<<endl;
}


// fix kernel radius, grid size, look at area on sphere with different radius
TEST(mc, kernel) {

	
    uniform_real_distribution<double> udist(0.0,1.0);

    GridData grid(GRIDNUM,GRIDBOUND);
    
    SurfaceStruct surf(grid);

    double rmin = GRIDBOUND/3.0;
    double rmax = 2 * rmin;
    int nr = 2;

    vector<double> rs = LinearSpacedArray(rmin, (rmax-rmin)/nr,nr);// different radius of sphere
    fmt::print("{}\n",rs);
        
    Kernel kernel(grid, 3*sqrt(grid.dx[0]));

    vector<double> err;


    for( double r : rs){
        double eA = 4 * M_PI * r * r;
        for( int i = 0; i < 10; i ++){
            // random center
            Coord center = {udist(rng),udist(rng),udist(rng)};
            create_sphere( surf, r, center);
            double aA = ComputeAreaGrid(surf.phi, kernel, grid);
            err.push_back(abs(aA-eA)/eA);
        }


        double sum = std::accumulate(err.begin(), err.end(), 0.0);
        double mae = sum / err.size();

        double sq_sum = std::inner_product(err.begin(), err.end(), err.begin(), 0.0);
        double std = std::sqrt(sq_sum / err.size() - mae * mae);

        fmt::print("{:>7.3f}, {:>7.3f} {:>7.3f}, {:<7.3f}\n",kernel.radius, r,mae,std);
        
    }

}


// fix radius of sphere, randomly shift center, look at convergence of area
TEST(mc, converge) {

    uniform_real_distribution<double> udist(0.0, 0.1);

    double gridbound = 1;
    double sphere_rad = 0.5;
    double exact_area = 4 * M_PI * sphere_rad * sphere_rad;
    vector<double> gridnums = {20,40,80};

    for (int gridnum : gridnums ){
        GridData grid(gridnum,gridbound);
        SurfaceStruct surf(grid);
        Kernel kernel;

        InitKernelCode(kernel, grid, 1.5*sqrt(grid.dx[0]),1);

        vector<double> relerr;
        fmt::print("{:<7} {:<7.2f}", gridnum, kernel.radius);

        for( int i = 0; i < 10; i ++){
            // random center
            Coord random_center = {udist(rng),udist(rng),udist(rng)};
            create_sphere( surf, sphere_rad, random_center);
            double apprx_area = ComputeAreaGrid(surf.phi, kernel, grid);
            relerr.push_back(abs(apprx_area - exact_area)/exact_area);
        }

        double sum = std::accumulate(relerr.begin(), relerr.end(), 0.0);
        double mae = sum / relerr.size();

        fmt::print("{:<7.4f}",mae);
    }

    fmt::print("\n");

}


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    CmdLine(argc, argv);
    return RUN_ALL_TESTS();
}

