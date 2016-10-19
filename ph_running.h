#ifndef PH_RUNNING_H
#define PH_RUNNING_H
#include <ph_system.h>
#include <ph_evolve.h>
#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include "time.h"
#include "ph_system.h"
#include "ph_evolve.h"
#include "ph_running.h"
using std::ofstream; using std::string;

class PH_Running
{
public:
    int N, maxit;
    double tolerance;
    ofstream coordinateFile, deltaEppFile, Eppexponentfile, betaFile, allthedataFile;
    ofstream diffFile;
    string filenamePrefix;
    PH_Evolve Iteration;

    // Initializing
    PH_Running();
    PH_Running(int N, int maxit, double tolerance, string filenamePrefix);


    void runit(double h1, double h2, double J);                 // To be called at each step


    void testit(double h1, double h2, double J);                // Test: Only runs one time.


    //Plotting functions
    void plot_rangeJ(double h1, double h2, double Jmin, double Jmax);
    //void plot_rangeJ_largeJ(double h1, double h2, double Jmin, double Jmax); // Only include this if things become REALLY SLOW.
    //void plot_rangeJ_normal(double h1, double h2, double Jmin, double Jmax); //This is basically just plot_rangeJ

    // Fixed fields
    void plot_hom_rangeh(double J, double hmin, double hmax);
    void plot_asym_rangeh(double J,double hmin, double hmax);

    // Random potential
    void plot_randomuniform(int averaging, double J, double hmin, double hmax);
    void plot_randomuniform_Jdivh_pretty(int averaging, double J, double hmin, double hmax);

    // Infinite temperature
    void plot_randomuniform_Jdivh_pretty_infinitetemp(int averaging, double J, double hmin, double hmax);
    void plot_hom_rangeh_infinitetemp(double J, double hmin, double hmax);
    void plot_asym_rangeh_infinitetemp(double J, double hmin, double hmax);
    void plot_lopsided_rangeh_infinitetemp(double J, double hmin, double hmax);

};

#endif // PH_RUNNING_H
