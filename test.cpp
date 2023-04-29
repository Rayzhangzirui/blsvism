#include <gtest/gtest.h>
#include "mc.h"
#include "globals.h"
#include "vism_addon.h"
#include "vismgpu.h"
using namespace std;

extern int KERNELRAD;

TEST(vism, cmdline){
    PrintOptions(cout);
    SoluteStruct sol = SolFromFile(SOLFILE);
    sol.Print(sol.dnum);
}


TEST(vism, methane) {
	GridData grid(GRIDNUM,GRIDBOUND);
    SoluteStruct sol = CreateMethane();

    double ***lj = matrix<double>(grid);
    
    SurfaceStruct surf(grid);

    KernelStruct kernel(grid, KERNELRAD);

    getinitsurf(surf, sol, grid);
    cal_LJ(sol, lj, grid);
    
    cout<<"initial surf = " << GSurf(surf,kernel,grid)<<endl;
    cout<<"initial ljswbox = " <<  GLjSolWatBox(surf, lj, grid)<<endl;
    
    binaryflowheaplistinterfaceonly(surf, sol, lj, kernel, grid);

    cout<<"final surf = " << GSurf(surf,kernel,grid)<<endl;
    cout<<"final ljswbox = " <<  GLjSolWatBox(surf, lj, grid)<<endl;
    cout<<"ljoutsidebox = " << LJoutsidebox3D(sol, grid);

    free_matrix<double>(lj, grid);
}


// test random number generator
TEST(mc, rand){
    ofstream of1("matlabc25/test_uniform_real.txt");
    ofstream of2("matlabc25/test_uniform_normal.txt");
    ofstream of3("matlabc25/test_uniform_int.txt");
    for(int i = 0; i < 10000; i ++){
        of1<<UniformReal()<<endl;
        of2<<UniformNormal()<<endl;
        of3<<UniformInt(0,10)<<endl;
    }
}


TEST(vism,lj){
    // check relateive error of LJ computation, cpu and ocl
    GridData grid(GRIDNUM,GRIDBOUND);
    SoluteStruct sol = SolFromFile(SOLFILE);

    double ***lj1 = matrix<double>(grid);
    double ***lj2 = matrix<double>(grid);

    parLJm2(sol, lj1, grid);
    seqLJ(sol, lj2, grid);
    double RELERR =  1e-9;
    for (int i = 0;i <= grid.nx[0]; i++){
        for (int j = 0;j <= grid.nx[1]; j++){
            for (int k = 0; k <= grid.nx[2]; k++){
                EXPECT_TRUE(abs((lj2[i][j][k]-lj1[i][j][k])/lj1[i][j][k])<RELERR);
            }
        }
    }

    free_matrix<double>(lj1, grid);
    free_matrix<double>(lj2, grid);
    
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    CmdLine(argc, argv);
    return RUN_ALL_TESTS();
}

