#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include "time.h"
#include "ph_system.h"
#include "ph_evolve.h"
#include "ph_running.h"

using namespace std;

void plotvaryingJ(int N, int maxit, double tolerance, double h1, double h2, double Jmin, double Jmax, string filename);
void plotrandomuniformh(int averaging, int n, int maxit, double tolerance, double J, double hmin, double hmax, string filename);
void plotrandomuniformhJdivhpretty(int averaging, int n, int maxit, double tolerance, double J, double hmin, double hmax, string filename);



int main()
{
    //Constants for initializing PH_Running
    int N = 301;
    int maxit = 1e7;
    double tolerance = 1e-6;



    //Constants for initializing PH_evolve
    // Comment out as needed
    /*
    //plot_rangeJ

    double h1 = 1;
    double h2 = 1;
    double Jmin = 0;
    double Jmax = 100;
    string filename_varyingJ = "hom_h1_J0to100_res_301";
    plotvaryingJ(N, maxit, tolerance, h1, h2, Jmin, Jmax, filename_varyingJ);
    */


    //plot_..._rangeh, plot_randomuniform
    double J = 1;
    double hmin = 0.25;
    double hmax = 1000.0;

    //plot_randomuniform
    int n = 1001;                   // This should be adequate
    int averaging = 500;           // Should this?
    string filename_randomh = "randompot_h0p25to1000_J1_av500_res1001_filtered_cranked_improved";
    plotrandomuniformhJdivhpretty(averaging, n, maxit, tolerance, J, hmin, hmax, filename_randomh);

}

void plotvaryingJ(int N, int maxit, double tolerance, double h1, double h2, double Jmin, double Jmax, string filename)
{
    PH_Running run(N, maxit, tolerance, filename);
    run.plot_rangeJ(h1, h2, Jmin, Jmax);
}

void plotrandomuniformh(int averaging, int n, int maxit, double tolerance, double J, double hmin, double hmax, string filename)
{
    PH_Running runrandompotential(n, maxit, tolerance, filename);
    runrandompotential.plot_randomuniform(averaging, J, hmin, hmax);
}

void plotrandomuniformhJdivhpretty(int averaging, int n, int maxit, double tolerance, double J, double hmin, double hmax, string filename)
{
    PH_Running runrandompotential(n, maxit, tolerance, filename);
    runrandompotential.plot_randomuniform_Jdivh_pretty(averaging, J, hmin, hmax);
}
