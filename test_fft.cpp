#include "globals.h"
#include "kernel2.h"
#include "surf.h"
#include "utiltest.h"
#include "fftconv.h"

#include <gtest/gtest.h>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/ostr.h> // must be included
#include <spdlog/stopwatch.h>
#include "spdlog/cfg/argv.h" // for loading levels from argv

// only consider lj field, the correct surface should include all the negative cells
TEST(fft, area){

    double kscale = 3;
    GridData grid(GRIDNUM,GRIDBOUND);

    SurfaceStruct surf(grid);
    KernelStruct kernel(grid, kscale * sqrt(grid.dx[0]));
    
    FFTWStruct fftw;

    getinitfftw(fftw, kernel, grid);


    Coord center = {0.0,0.0,0.0};
    double r = 1.0;
    create_sphere( surf, r, center);
    

    double Afft = getlengthfftw(surf.phi, fftw, grid);

    double Aconv = getlengthlist(surf.phi, kernel, grid);

    fmt::print("fft area {:<7.4f}, conv area {:<7.4f}", Afft, Aconv);


}


int main(int argc, char **argv) {
    
    spdlog::set_level(spdlog::level::info); // Set global log level to debug
    spdlog::set_pattern("%v");
    spdlog::cfg::load_argv_levels(argc, argv);
    
    testing::InitGoogleTest(&argc, argv);

    CmdLine(argc, argv);
    return RUN_ALL_TESTS();
}
