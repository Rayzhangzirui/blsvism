/* 
utitility for testing
*/
#ifndef TESTUTIL_H
#define TESTUTIL_H



double VismOptimalRadius( double sigma, double eps, double Q, double gamma0, double tau, double rhow, double epsm, double epsw){

    const double A = 22.11; // 4*pi*A = 1/(8*pi*0.00014319) from Wang et al 2012, one charged Particle equation
    double r = 1;
    double pi = M_PI;
    for (int k = 0; k < 30; k++){
        double Gp = 4*pi*(2*gamma0*r - 2*gamma0*tau + 4*rhow*eps*(pow(sigma,6)/pow(r,4)-pow(sigma,12)/pow(r,10))
                    +(A*Q*Q/pow(r,2))*(1/epsm-1/epsw));
        double Gpp = 4*pi*(2*gamma0+4*rhow*eps*(10*pow(sigma,12)/pow(r,11)-4*pow(sigma,6)/pow(r,5))
                    -(2*A*Q*Q/pow(r,3))*(1/epsm-1/epsw));
        r = r - Gp/Gpp;
    }

    return r;
}



// return multiple value 
// https://stackoverflow.com/questions/321068/returning-multiple-values-from-a-c-function
auto VismSphereEnergy(double r, double sigma, double eps, double Q, double gamma0, double tau, double rhow, double epsm, double epsw){
    const double A = 22.11; // 4*pi*A = 1/(8*pi*0.00014319) from Wang et al 2012, one charged Particle equation
    Energy e;
    struct solution{
        double Gsurf;
        double Gvdw;
        double Gelec;
        double Gtot;
        
    };

    double Gsurf = 4*M_PI*gamma0*(pow(r,2) - 2*tau*r);
    double Gvdw = 4*M_PI*4*rhow*eps*((pow(sigma,12)/9)/(pow(r,9))-(pow(sigma,6)/3)/(pow(r,3)));
    double Gelec = -(4*M_PI*A*Q*Q/r)*(1/epsm-1/epsw);
    double Gtot = 4*M_PI*(gamma0*(pow(r,2) - 2*tau*r)+4*rhow*eps*((pow(sigma,12)/9)/(pow(r,9))-(pow(sigma,6)/3)/(pow(r,3)))-(A*Q*Q/r)*(1/epsm-1/epsw));

    return solution{Gsurf,Gvdw,Gelec,Gtot};
}




// create sphere radius and center
void create_sphere(SurfaceStruct& surf, double radius, Coord center){
    #pragma omp parallel for collapse(3)
    for(int i = 0; i <= surf.grid.nx[0]; i++){
        for(int j = 0; j <= surf.grid.nx[1]; j++){
            for(int k = 0; k <= surf.grid.nx[2]; k++){
                array<double,3> x = sub2coord(array<int,3>{i,j,k},surf.grid);
                surf.phi[i][j][k]= (dist2(x,center)<radius)?-1:1;
            }
        }
    }
}




//a[0] = a, a[1] = a + h, ... a[N+1] = a + Nh
template <typename T>
inline std::vector<T> LinearSpacedArray(T a, T h, std::size_t N)
{
    std::vector<T> xs(N+1);
    for (int i = 0; i <= N; i ++) {
        xs[i] = a + h * i;
    }
    return xs;
}

//a[0] = a, a[1] = a*r, ... a[N+1] = a*r^N
template <typename T>
inline std::vector<T> GeometricArray(T a, T r, std::size_t N)
{
    std::vector<T> xs(N+1);
    xs[0] = a;
    for (int i = 1; i <= N; i ++) {
        xs[i] = xs[i-1]*r;
    }
    return xs;
}

#endif