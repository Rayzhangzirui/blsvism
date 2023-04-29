#include "writer.h"
#include <chrono>

void IncProbCount(double*** prob, double*** phi, GridData& grid){
    // increment counter
    #pragma omp parallel for collapse(3)
    for(int i = 0; i <= grid.nx[0]; i++){
        for(int j = 0; j <= grid.nx[1]; j++){
            for(int k = 0; k <= grid.nx[2]; k++){
                if(phi[i][j][k] > 0){
                    prob[i][j][k] = prob[i][j][k] + 1;
                }
            }
        }
    }

}

// multiply the field by scalar
void MultiplyScalar(double ***f, double s, GridData& grid){
    #pragma omp parallel for collapse(3)
    for(int i = 0; i <= grid.nx[0]; i++){
        for(int j = 0; j <= grid.nx[1]; j++){
            for(int k = 0; k <= grid.nx[2]; k++){
                f[i][j][k] = f[i][j][k] * s;
            }
        }
    }
}




//https://stackoverflow.com/questions/22590821/convert-stdduration-to-human-readable-time
std::string beautify_duration(std::chrono::duration<double> input_seconds)
{
    using namespace std::chrono;
    typedef duration<int, std::ratio<86400>> days;
    auto d = duration_cast<days>(input_seconds);
    input_seconds -= d;
    auto h = duration_cast<hours>(input_seconds);
    input_seconds -= h;
    auto m = duration_cast<minutes>(input_seconds);
    input_seconds -= m;
    auto s = duration_cast<seconds>(input_seconds);

    auto dc = d.count();
    auto hc = h.count();
    auto mc = m.count();
    auto sc = s.count();

    std::stringstream ss;
    ss.fill('0');
    if (dc) {
        ss << d.count() << "d";
    }
    if (dc || hc) {
        if (dc) { ss << std::setw(2); } //pad if second set of numbers
        ss << h.count() << "h";
    }
    if (dc || hc || mc) {
        if (dc || hc) { ss << std::setw(2); }
        ss << m.count() << "m";
    }
    if (dc || hc || mc || sc) {
        if (dc || hc || mc) { ss << std::setw(2); }
        ss << s.count() << 's';
    }
    return ss.str();
}



/*
Writer class to deal with writing data
info_intv: interval to output info
write_surf_intv: interval to write surface
tiem_est_n: how frequentyl output estimated finish time
record_prob: record probability surface
*/

Writer::Writer(string modelname, int info_intv, int write_surf_intv, int time_est_n, bool record_prob, GridData grid)
:modelname(modelname),
info_intv(info_intv),
write_surf_intv(write_surf_intv),
time_est_n(time_est_n),
record_prob(record_prob),
grid(grid),
out(NULL)
{

    begin = std::chrono::steady_clock::now();

    // if OUTPUTFILE is empty or cout, ouput to cout
    //https://stackoverflow.com/questions/366955/obtain-a-stdostream-either-from-stdcout-or-stdofstreamfile
    if (OUTPUTFILE.empty() || OUTPUTFILE=="cout") {
        buf = std::cout.rdbuf();
    
    } else {
        of.open(modelname+"_eng.txt");
        buf = of.rdbuf();
    }
    
    out.rdbuf(buf);
    
    prob = matrix<double>(grid);

    PrintOptions(out);
    // write header
    out <<"# ljss bond bend surf ljswbox ljswout total rg"<<endl;
}




void Writer::WriteEnergy(SoluteStruct &sol, Energy &e){    
    out << std::fixed<<setprecision(3)
        << e.ljss<<" "<<e.bond<<" "<<e.bend<<" "<<e.surf<<" "<<e.ljswbox<<" "<<e.ljswout<<" "<<e.total()<<" "
        << sol.RadiusGyration()<<endl;
}

// record and write probability of bianry surface
void Writer::RecordProbSurf(int k, SurfaceStruct &surf){

    IncProbCount(prob, surf.phi, grid);

    if ( write_surf_intv>=1 && k%write_surf_intv==0 ){
        string probname = modelname + "_prob_"+ to_string(k); // name of prob
        MultiplyScalar(prob, (double) 1/write_surf_intv, grid);
        write_field(probname, grid, prob);
        set_field<double>(grid, prob, 0); // reset to zero
    }
}


void Writer::WritePosAndSurface(int k, SoluteStruct &sol, SurfaceStruct &surf){
    // write surface and position
    if ( write_surf_intv>=1 && k%write_surf_intv==0 ){
        string surfname = modelname + "_surf_"+ to_string(k); // name of surface file
        write_field(surfname, surf.grid, surf.phi);

        string posname = modelname + "_pos_"+ to_string(k); // name of position
        sol.Write(posname);
    }
}

void Writer::WriteOutput(int k, SoluteStruct &sol, SurfaceStruct &surf, Energy &e){
    // output energy
    if(info_intv>=1 && k%info_intv==0){
        WriteEnergy(sol, e);
    }

    // estimate remaining time
    bool every_est = (k%(MAXSTEP/time_est_n)==0);
    if( (k==MAXSTEP/1000) || every_est) {
        // current time
        auto now = std::chrono::system_clock::now();
        std::time_t now_time = std::chrono::system_clock::to_time_t(now);
        cout<<put_time(localtime(&now_time),"%Y-%m-%d %X, ");
        // output time estimation 100 times
        auto end = std::chrono::steady_clock::now(); // current time
        std::chrono::duration<double> duration = end-begin;
        cout << "step "<< k<< ", elapsed " << beautify_duration(duration) <<", est " << beautify_duration(duration*(MAXSTEP-k)/k) <<endl;
    }

    
    WritePosAndSurface(k, sol, surf);

    if(record_prob){
        RecordProbSurf(k, surf);
    }

}

// destructor, output total time
Writer::~Writer(){
    auto end = std::chrono::steady_clock::now(); // current time
    std::chrono::duration<double> duration = end-begin;
    out << "total time " << beautify_duration(duration) <<endl;

    free_matrix<double>(prob, grid);
}


