#ifndef WRITER_H
#define WRITER_H

#include "mcenergy.h"


struct Writer{
    GridData grid;
    int info_intv = 1; // interval to output energy
    int write_surf_intv = 100; // interval to output position
    int time_est_n = 10; // estimate time every 1/10 of max step
    bool record_prob = false;
    double ***prob; // probability counter
    std::chrono::steady_clock::time_point begin; // start time
    string modelname;
    std::streambuf *buf;
    std::ofstream of;
    ostream out;

    Writer(string modelname, int info_intv, int write_surf_intv, int time_est_intv, bool record_prob, GridData grid);

    void WriteEnergy(SoluteStruct &sol, Energy &e);

    void RecordProbSurf(int k, SurfaceStruct &surf);

    void WriteOutput(int current_iter, SoluteStruct &sol, SurfaceStruct &surf, Energy &e);

    void WritePosAndSurface(int k, SoluteStruct &sol, SurfaceStruct &surf);

    ~Writer();
};

#endif