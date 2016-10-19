#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include "time.h"
#include "ph_system.h"
#include "ph_evolve.h"
#include "ph_running.h"

using namespace std;

void testweights(double h1, double h2, double J);
void plotvaryingJ(int N, int maxit, double tolerance, double h1, double h2, double Jmin, double Jmax, string filename);
void plotrandomuniformh(int averaging, int n, int maxit, double tolerance, double J, double hmin, double hmax, string filename);
void plotrandomuniformhJdivhpretty(int averaging, int n, int maxit, double tolerance, double J, double hmin, double hmax, string filename);
void plotrandomuniformhJdivhpretty_inftemp(int averaging, int n, int maxit, double tolerance, double J, double hmin, double hmax, string filename);
void plothomJdivhpretty_inftemp(int n, int maxit, double tolerance, double J, double hmin, double hmax, string filename);
void plotasymJdivhpretty_inftemp(int n, int maxit, double tolerance, double J, double hmin, double hmax, string filename);
void plotlopsJdivhpretty_inftemp(int n, int maxit, double tolerance, double J, double hmin, double hmax, string filename);
void plotrandomuniformhJdivhpretty_inftemp_sorted(int averaging, int n, int maxit, double tolerance, double J, double hmin, double hmax, string filename);


int main()
{   
    //Constants for initializing PH_Running
    //int N = 301;
    int maxit = 1e7;
    double tolerance = 1e-6;


    //PH_Running testrun(N, maxit, tolerance, "test");
    //testrun.testit(1,1,1);

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

    /* */

    //plot_..._rangeh, plot_randomuniform
    double J = 1;
    double hmin = 0.25;             //
    double hmax = 1000.0;

    //plot_randomuniform
    int n = 101;                   // Adequate? I don't need that high of a resolution...
    int averaging = 500;           // This is a bit too small, probably, but I am testing something.
    //string filename_randomh = "randompot_h1p6to2p5_J1_av100_res101_161016_test2";
    //string filename_randomh = "randompot_h0p0001to1000_J1_av1000_res1001_171016";
    string filename_randomh = "test_random_sorted_191016";
    //plotrandomuniformhJdivhpretty(averaging, n, maxit, tolerance, J, hmin, hmax, filename_randomh);
    plotrandomuniformhJdivhpretty(averaging, n, maxit, tolerance, J, hmin, hmax, filename_randomh);
    //plotrandomuniformhJdivhpretty_inftemp_sorted(averaging, n, maxit, tolerance, J, hmin, hmax, filename_randomh);
    //plotrandomuniformh(averaging, n, maxit, tolerance, J, hmin, hmax, filename_randomh);
    //plotlopsJdivhpretty_inftemp(n, maxit, tolerance, J, hmin, hmax, filename_randomh);




    /*
    int N = 1;
    double J  = 1.0;
    double h = 10000;
    //double diff1 = 0;
    //double diff2 = 0;
    //double h2 = -3.0;
    double av1, av2; // For averages
    double h1, h2;
    h2 = 0;

    //std::default_random_engine generator;                       // I asked the internet, and it replied
    //std::uniform_real_distribution<double> distribution(-h,h);
    for(int i=0; i<N; i++)
    {
        //double h1 = distribution(generator);
        //double h2 = distribution(generator);

        //cout << h1 << " " << h2 << endl;

        h1 = h;
        h2 = -h;

        bool finitetemp = false;

        PH_Evolve evolving(maxit, h1, h2, J, tolerance, finitetemp);



        //av1 += testsystem.walpha1;
        //av2 += testsystem.walpha2;

    }

    //av1 = av1/N;
    //av2 = av2/N;

    //diff1 = diff1/N;
    //diff2 = diff2/N;

    //cout << "Average weight: " << endl;
    //cout << "State 1: " << av1 << endl;
    //cout << "State 2: " << av2 << endl;

    //cout << "Average difference: " << endl;
    //cout << "State 1: " << diff1 << endl;
    //cout << "State 2: " << diff2 << endl;
    */





}

void testweights(double h1, double h2, double J)
{
    PH_System testsystem(h1, h2, J);

    cout << endl;
    cout << "h1 = " <<  h1 << "; h2 = " << h2 << "; J = " << J << endl;
    cout << "E++ = " <<  testsystem.Epp << "; lambda1 = " << testsystem.lambda1 << "; lambda2 = " << testsystem.lambda2 << "; Emm = " << testsystem.Emm << endl;
    cout << "State 1:" << endl;
    //cout << "From my first attempt: " << endl;
    cout << "Weight in matrices: " << testsystem.walpha1 << endl;
    //cout << "Weigth in thermal matrix: " << testsystem.wp1 << endl;
    //cout << "From my second attempt: " << endl;
    //cout << "Weight in thermal matrix: " << testsystem.walpha1 << endl;
    //cout << "Difference: " << (testsystem.walpha1-testsystem.wp1) << endl;


    cout << endl;
    cout << "State 2:" << endl;
    cout << "Weigth in matrices: " << testsystem.walpha2 << endl;
    //cout << "Weigth in thermal matrix: " << testsystem.wp2 << endl;
    //cout << "From my second attempt: " << endl;
    //cout << "Weight in thermal matrix: " << testsystem.walpha2 << endl;
    //cout << "Difference: " << (testsystem.walpha2-testsystem.wp2) << endl;



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

void plotrandomuniformhJdivhpretty_inftemp(int averaging, int n, int maxit, double tolerance, double J, double hmin, double hmax, string filename)
{
    PH_Running runrandompotential(n, maxit, tolerance, filename);
    runrandompotential.plot_randomuniform_Jdivh_pretty_infinitetemp(averaging, J, hmin, hmax);
}

void plothomJdivhpretty_inftemp(int n, int maxit, double tolerance, double J, double hmin, double hmax, string filename)
{
    PH_Running runhompotential(n, maxit, tolerance, filename);
    runhompotential.plot_hom_rangeh_infinitetemp(J, hmin, hmax);
}


void plotasymJdivhpretty_inftemp(int n, int maxit, double tolerance, double J, double hmin, double hmax, string filename)
{
    PH_Running runhompotential(n, maxit, tolerance, filename);
    runhompotential.plot_asym_rangeh_infinitetemp(J, hmin, hmax);
}

void plotlopsJdivhpretty_inftemp(int n, int maxit, double tolerance, double J, double hmin, double hmax, string filename)
{
    PH_Running runhompotential(n, maxit, tolerance, filename);
    runhompotential.plot_lopsided_rangeh_infinitetemp(J, hmin, hmax);
}


void plotrandomuniformhJdivhpretty_inftemp_sorted(int averaging, int n, int maxit, double tolerance, double J, double hmin, double hmax, string filename)
{
    PH_Running runrandompotential(n, maxit, tolerance, filename);
    runrandompotential.plot_randomuniform_Jdivh_pretty_infinitetemp_sorted(averaging, J, hmin, hmax);
}
