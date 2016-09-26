#include "ph_running.h"
#include <iostream>
#include <iomanip>
#include <armadillo>
#include <random>

using namespace std;
using namespace arma;

// So I guess I should do something like what I did in Python and plot J/|h|... Make an own function for that?
PH_Running::PH_Running(int N, int maxit,double tolerance, string filenamePrefix)
{
    this->tolerance=tolerance;
    this->N        =N;
    this->maxit    =maxit;
    this->filenamePrefix = filenamePrefix;
}


//Must print to file here somewhere. Not in running, maybe?
void PH_Running::runit(double h1, double h2, double J)
{
    Iteration = PH_Evolve(maxit, h1, h2, J, tolerance);
}





void PH_Running::plot_hom_rangeh(double J, double hmin, double hmax)
{
    double h1, h2, fraction;
    double a = clock();
    // Open coordinate file
    this->filenamePrefix = filenamePrefix;
    char *filename = new char[1000];                                    // File name can have max 1000 characters
    sprintf(filename, "%s_coordinateFile.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    coordinateFile.open(filename);
    delete filename;

    vec hs = zeros(N);

    hs = linspace(hmin, hmax, N);
    fraction = (N-1)/4.;
    for(int i=0; i<N; i++){
        if(i==fraction)                   cout <<"25%, time:" << clock()-a << endl;
        if(i==fraction*2)                 cout <<"50%, time:" << clock()-a << endl;
        if(i==fraction*3)                 cout <<"75%, time:" << clock()-a << endl;
        h1 =  hs[i];
        h2 =  hs[i];
        runit(h1, h2, J);
        coordinateFile << J << " " << abs(Iteration.delta_rho_Epp) << " " << abs(Iteration.delta_rho_Emm) << " " << abs(Iteration.delta_rho_2nd) << " " << abs(Iteration.delta_rho_3rd) << " " << endl;
      }
    coordinateFile.close();
}

void PH_Running::plot_rangeJ(double h1, double h2, double Jmin, double Jmax)
{
    double J, fraction;
    double a = clock();
    // Open coordinate file
    char *filename = new char[1000];                                    // File name can have max 1000 characters
    sprintf(filename, "%s_coordinateFile.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    coordinateFile.open(filename);
    delete filename;
    vec Js = zeros(N);
    Js = linspace(Jmin, Jmax, N);  // arma::linspace ?
    fraction = floor((N-1)/4.0);                   // Unit to check progress. Flooring just in case.
    for(int i=0; i<N; i++){
        if (i==fraction)                cout << "25%, time:" << clock()-a << endl;
        if (i==fraction*2)              cout << "50 %, time:"<< clock()-a << endl;
        if (i==fraction*3)              cout << "75%, time:" << clock()-a << endl;
        //Done with printing
        J = Js[i];
        runit(h1, h2, J);
        coordinateFile << J << " " << abs(Iteration.delta_rho_Epp) << " " << abs(Iteration.delta_rho_Emm) << " " << abs(Iteration.delta_rho_2nd) << " " << abs(Iteration.delta_rho_3rd) << " " << endl;}
    coordinateFile.close();
    }

void PH_Running::plot_asym_rangeh(double J, double hmin, double hmax)
{
    double h1, h2, fraction;
    double a = clock();

    // Open coordinate file
    this->filenamePrefix = filenamePrefix;
    char *filename = new char[1000];                                    // File name can have max 1000 characters
    sprintf(filename, "%s_coordinateFile.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    coordinateFile.open(filename);
    delete filename;

    vec hs = zeros(N);

    hs = linspace(hmin, hmax, N);
    fraction = (N-1)/4.;
    for(int i=0; i<N; i++){
        if(i==fraction)                   cout <<"25%, time:" << clock()-a << endl;
        if(i==fraction*2)                 cout <<"50%, time:" << clock()-a << endl;
        if(i==fraction*3)                 cout <<"75%, time:" << clock()-a << endl;
        h1 =   hs[i];
        h2 = - hs[i];
        runit(h1, h2, J);
        coordinateFile << J << " " << abs(Iteration.delta_rho_Epp) << " " << abs(Iteration.delta_rho_Emm) << " " << abs(Iteration.delta_rho_2nd) << " " << abs(Iteration.delta_rho_3rd) << " " << endl;
      }
    coordinateFile.close();
}


void PH_Running::plot_randomuniform(int averaging, double J, double hmin, double hmax)
{

    // Open file for printing delta_rhos to file
    this->filenamePrefix = filenamePrefix;
    char *filename1 = new char[1000];                                    // File name can have max 1000 characters
    sprintf(filename1, "%s_coordinateFile.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    coordinateFile.open(filename1);
    delete filename1;

    // Open file for errors?
    //char *filename2 = new char[1000];                               // File name can have max 1000 characters
    //sprintf(filename2, "%s_diffFile.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    //diffFile.open(filename2);
    //delete filename2;

    vec hs;
    double h, counter;


    hs = linspace(hmin, hmax, N);
    double rmsdelta_rho_Epp, rmsdelta_rho_Emm, rmsdelta_rho_2nd, rmsdelta_rho_3rd;
    double variance_delta_rho_Epp, variance_delta_rho_Emm, variance_delta_rho_2nd, variance_delta_rho_3rd;
    double stddelta_rho_Epp, stddelta_rho_Emm, stddelta_rho_2nd, stddelta_rho_3rd;
    double absmean_delta_rho_Epp, absmean_delta_rho_Emm, absmean_delta_rho_2nd, absmean_delta_rho_3rd;
    for(int i=0; i<N; i++)
    {
        cout << "i= " << i << endl;
        counter = 0;       // Should we fix the number of contributions?
        rmsdelta_rho_Epp = 0;
        rmsdelta_rho_Emm = 0;
        rmsdelta_rho_2nd = 0;
        rmsdelta_rho_3rd = 0;
        variance_delta_rho_Epp = 0;
        variance_delta_rho_Emm = 0;
        variance_delta_rho_2nd = 0;
        variance_delta_rho_3rd = 0;
        stddelta_rho_Epp = 0;
        stddelta_rho_Emm = 0;
        stddelta_rho_2nd = 0;
        stddelta_rho_3rd = 0;
        absmean_delta_rho_Epp = 0;
        absmean_delta_rho_Emm = 0;
        absmean_delta_rho_2nd = 0;
        absmean_delta_rho_3rd = 0;
        h = hs[i];
        std::default_random_engine generator;                       // I asked the internet, and it replied
        std::uniform_real_distribution<double> distribution(-h,h);
        while(counter<averaging)
        {
            double h1 = distribution(generator);
            double h2 = distribution(generator);
            runit(h1, h2, J);
            if(Iteration.nm_diff < 100*tolerance)
            {
            rmsdelta_rho_Epp += Iteration.delta_rho_Epp*Iteration.delta_rho_Epp;
            rmsdelta_rho_Emm += Iteration.delta_rho_Emm*Iteration.delta_rho_Emm;
            rmsdelta_rho_2nd += Iteration.delta_rho_2nd*Iteration.delta_rho_2nd;
            rmsdelta_rho_3rd += Iteration.delta_rho_3rd*Iteration.delta_rho_3rd;
            variance_delta_rho_Epp += Iteration.delta_rho_Epp;
            variance_delta_rho_Emm += Iteration.delta_rho_Emm;
            variance_delta_rho_2nd += Iteration.delta_rho_2nd;
            variance_delta_rho_3rd += Iteration.delta_rho_3rd;
            absmean_delta_rho_Epp += abs(Iteration.delta_rho_Epp);
            absmean_delta_rho_Emm += abs(Iteration.delta_rho_Emm);
            absmean_delta_rho_2nd += abs(Iteration.delta_rho_2nd);
            absmean_delta_rho_3rd += abs(Iteration.delta_rho_3rd);
            counter++;
            }
        }  // End while loop (samples of systems with h1,h2 from random uniform dist from [-h,h])

        // Variance (rmsdelta is still just the sum over the squares)
        // It's weird that I need to take the absolute value here. Something wrong?
        // But this is not really the variance? ... So what is it? is it useful at all?
        variance_delta_rho_Epp = abs(rmsdelta_rho_Epp - variance_delta_rho_Epp*variance_delta_rho_Epp)/averaging;
        variance_delta_rho_Emm = abs(rmsdelta_rho_Emm - variance_delta_rho_Emm*variance_delta_rho_Emm)/averaging;
        variance_delta_rho_2nd = abs(rmsdelta_rho_2nd - variance_delta_rho_2nd*variance_delta_rho_2nd)/averaging;
        variance_delta_rho_3rd = abs(rmsdelta_rho_3rd - variance_delta_rho_3rd*variance_delta_rho_3rd)/averaging;
        // Standard deviation
        // Not really the standard deviation?
        stddelta_rho_Epp = sqrt(variance_delta_rho_Epp);       // Could calculate this after reading from file
        stddelta_rho_Emm = sqrt(variance_delta_rho_Emm);
        stddelta_rho_2nd = sqrt(variance_delta_rho_2nd);
        stddelta_rho_3rd = sqrt(variance_delta_rho_3rd);
        // Root mean square
        rmsdelta_rho_Epp = sqrt(rmsdelta_rho_Epp/averaging);
        rmsdelta_rho_Emm = sqrt(rmsdelta_rho_Emm/averaging);
        rmsdelta_rho_2nd = sqrt(rmsdelta_rho_2nd/averaging);
        rmsdelta_rho_3rd = sqrt(rmsdelta_rho_3rd/averaging);
        // Mean
        absmean_delta_rho_Epp = absmean_delta_rho_Epp/averaging;
        absmean_delta_rho_Emm = absmean_delta_rho_Emm/averaging;
        absmean_delta_rho_2nd = absmean_delta_rho_2nd/averaging;
        absmean_delta_rho_3rd = absmean_delta_rho_3rd/averaging;



        //Print to file. This may be altered.
        // Print J/h instead? Or h/J?
        coordinateFile << J/h << " " << rmsdelta_rho_Epp << " " << rmsdelta_rho_Emm << " " << rmsdelta_rho_2nd << " " <<rmsdelta_rho_3rd << " " << variance_delta_rho_Epp << " " << variance_delta_rho_Emm << " " << variance_delta_rho_2nd << " " << variance_delta_rho_3rd << " " << stddelta_rho_Epp << " " << stddelta_rho_Emm << " " << stddelta_rho_2nd << " " << stddelta_rho_3rd << " " << absmean_delta_rho_Epp << " " << absmean_delta_rho_Emm << " " << absmean_delta_rho_2nd << " " << absmean_delta_rho_3rd << endl;
        // Must find a smarter way to print...

    } // End loop over i (values of h)
    coordinateFile.close();
    //diffFile.close();               // Will probably not use diffFile
}

void PH_Running::plot_randomuniform_Jdivh_pretty(int averaging, double J, double hmin, double hmax)
{

    // Open file for printing delta_rhos to file
    this->filenamePrefix = filenamePrefix;
    char *filename1 = new char[1000];                                    // File name can have max 1000 characters
    sprintf(filename1, "%s_coordinateFile.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    coordinateFile.open(filename1);
    delete filename1;

    // Open file for errors?
    //char *filename2 = new char[1000];                               // File name can have max 1000 characters
    //sprintf(filename2, "%s_diffFile.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    //diffFile.open(filename2);
    //delete filename2;

    vec Jdivhssmall, Jdivhslarge, hs;
    double h, counter;

    double minJ = J/hmax;
    double maxJ = J/hmin;
    // Must have minJ<0.5 && maxJ > 0.6.
    int littleN = floor(0.9*N);
    int deltaN = N - littleN;
    Jdivhssmall = linspace(minJ, 0.65, littleN);
    Jdivhslarge = linspace(0.65, maxJ, deltaN);
    hs = zeros(N);
    for(int j=0; j<littleN; j++)        hs[j] = 1./Jdivhssmall[j];
    int j = littleN-1;
    for(int k=0; k<deltaN; k++)
    {
        j++;
        hs[j] = 1./Jdivhslarge[k];
    }
    // Could have more values for small Jdivhs, but that would take some time to make and result in ugly code...

    double rmsdelta_rho_Epp, rmsdelta_rho_Emm, rmsdelta_rho_2nd, rmsdelta_rho_3rd;
    double variance_delta_rho_Epp, variance_delta_rho_Emm, variance_delta_rho_2nd, variance_delta_rho_3rd;
    double stddelta_rho_Epp, stddelta_rho_Emm, stddelta_rho_2nd, stddelta_rho_3rd;
    double absmean_delta_rho_Epp, absmean_delta_rho_Emm, absmean_delta_rho_2nd, absmean_delta_rho_3rd;
    for(int i=0; i<N; i++)
    {
        cout << "i= " << i << endl;
        counter = 0;       // Should we fix the number of contributions?
        rmsdelta_rho_Epp = 0;
        rmsdelta_rho_Emm = 0;
        rmsdelta_rho_2nd = 0;
        rmsdelta_rho_3rd = 0;
        variance_delta_rho_Epp = 0;
        variance_delta_rho_Emm = 0;
        variance_delta_rho_2nd = 0;
        variance_delta_rho_3rd = 0;
        stddelta_rho_Epp = 0;
        stddelta_rho_Emm = 0;
        stddelta_rho_2nd = 0;
        stddelta_rho_3rd = 0;
        absmean_delta_rho_Epp = 0;
        absmean_delta_rho_Emm = 0;
        absmean_delta_rho_2nd = 0;
        absmean_delta_rho_3rd = 0;
        h = hs[i];
        std::default_random_engine generator;                       // I asked the internet, and it replied
        std::uniform_real_distribution<double> distribution(-h,h);
        while(counter<averaging)
        {
            double h1 = distribution(generator);
            double h2 = distribution(generator);
            runit(h1, h2, J);
            if(Iteration.nm_diff < 100*tolerance)
            {
            rmsdelta_rho_Epp += Iteration.delta_rho_Epp*Iteration.delta_rho_Epp;
            rmsdelta_rho_Emm += Iteration.delta_rho_Emm*Iteration.delta_rho_Emm;
            rmsdelta_rho_2nd += Iteration.delta_rho_2nd*Iteration.delta_rho_2nd;
            rmsdelta_rho_3rd += Iteration.delta_rho_3rd*Iteration.delta_rho_3rd;
            variance_delta_rho_Epp += Iteration.delta_rho_Epp;
            variance_delta_rho_Emm += Iteration.delta_rho_Emm;
            variance_delta_rho_2nd += Iteration.delta_rho_2nd;
            variance_delta_rho_3rd += Iteration.delta_rho_3rd;
            absmean_delta_rho_Epp += abs(Iteration.delta_rho_Epp);
            absmean_delta_rho_Emm += abs(Iteration.delta_rho_Emm);
            absmean_delta_rho_2nd += abs(Iteration.delta_rho_2nd);
            absmean_delta_rho_3rd += abs(Iteration.delta_rho_3rd);
            counter++;
            }
        }  // End while loop (samples of systems with h1,h2 from random uniform dist from [-h,h])

        // Variance (rmsdelta is still just the sum over the squares)
        // It's weird that I need to take the absolute value here. Something wrong?
        // But this is not really the variance? ... So what is it? is it useful at all?
        variance_delta_rho_Epp = abs(rmsdelta_rho_Epp - variance_delta_rho_Epp*variance_delta_rho_Epp)/averaging;
        variance_delta_rho_Emm = abs(rmsdelta_rho_Emm - variance_delta_rho_Emm*variance_delta_rho_Emm)/averaging;
        variance_delta_rho_2nd = abs(rmsdelta_rho_2nd - variance_delta_rho_2nd*variance_delta_rho_2nd)/averaging;
        variance_delta_rho_3rd = abs(rmsdelta_rho_3rd - variance_delta_rho_3rd*variance_delta_rho_3rd)/averaging;
        // Standard deviation
        // Not really the standard deviation?
        stddelta_rho_Epp = sqrt(variance_delta_rho_Epp);       // Could calculate this after reading from file
        stddelta_rho_Emm = sqrt(variance_delta_rho_Emm);
        stddelta_rho_2nd = sqrt(variance_delta_rho_2nd);
        stddelta_rho_3rd = sqrt(variance_delta_rho_3rd);
        // Root mean square
        rmsdelta_rho_Epp = sqrt(rmsdelta_rho_Epp/averaging);
        rmsdelta_rho_Emm = sqrt(rmsdelta_rho_Emm/averaging);
        rmsdelta_rho_2nd = sqrt(rmsdelta_rho_2nd/averaging);
        rmsdelta_rho_3rd = sqrt(rmsdelta_rho_3rd/averaging);
        // Mean
        absmean_delta_rho_Epp = absmean_delta_rho_Epp/averaging;
        absmean_delta_rho_Emm = absmean_delta_rho_Emm/averaging;
        absmean_delta_rho_2nd = absmean_delta_rho_2nd/averaging;
        absmean_delta_rho_3rd = absmean_delta_rho_3rd/averaging;



        //Print to file. This may be altered.
        // Print J/h instead? Or h/J?
        coordinateFile << hs[i]/J << " " << rmsdelta_rho_Epp << " " << rmsdelta_rho_Emm << " " << rmsdelta_rho_2nd << " " <<rmsdelta_rho_3rd << " " << variance_delta_rho_Epp << " " << variance_delta_rho_Emm << " " << variance_delta_rho_2nd << " " << variance_delta_rho_3rd << " " << stddelta_rho_Epp << " " << stddelta_rho_Emm << " " << stddelta_rho_2nd << " " << stddelta_rho_3rd << " " << absmean_delta_rho_Epp << " " << absmean_delta_rho_Emm << " " << absmean_delta_rho_2nd << " " << absmean_delta_rho_3rd << endl;
        // Must find a smarter way to print...

    } // End loop over i (values of h)
    coordinateFile.close();
    //diffFile.close();               // Will probably not use diffFile
}

/*

void PH_Running::plot_rangeJ_normal(double h1, double h2, double Jmin, double Jmax)
{
    double J, fraction;
    double a = clock();
    // Open coordinate file
    char *filename = new char[1000];                                    // File name can have max 1000 characters
    sprintf(filename, "%s_coordinateFile.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    coordinateFile.open(filename);
    delete filename;
    vec Js = zeros(N);
    Js = linspace(Jmin, Jmax, N);  // arma::linspace ?
    fraction = floor((N-1)/4.0);                   // Unit to check progress. Flooring just in case.
    for(int i=0; i<N; i++){
        if (i==fraction)                cout << "25%, time:" << clock()-a << endl;
        if (i==fraction*2)              cout << "50 %, time:"<< clock()-a << endl;
        if (i==fraction*3)              cout << "75%, time:" << clock()-a << endl;
        //Done with printing
        J = Js[i];
        runit(h1, h2, J);
        coordinateFile << J << " " << Iteration.delta_rho_Epp << " " << Iteration.delta_rho_Emm << " " << Iteration.delta_rho_2nd << " " << Iteration.delta_rho_3rd << " " << endl;}
    coordinateFile.close();
    }

void PH_Running::plot_rangeJ_largeJ(double h1, double h2, double Jmin, double Jmax)
{
    double J, fraction;

    double a = clock();
    // Open coordinate file
    char *filename = new char[1000];                                     // File name can have max 1000 characters
    sprintf(filename, "%s_coordinateFile.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    coordinateFile.open(filename);
    delete filename;

    int smallN = int(floor(0.9*N));
    int deltaN = N - smallN;
    vec smallJs = zeros(smallN);
    vec largeJs = zeros(deltaN);
    smallJs = linspace(Jmin, 5, smallN);                                 // Granted, this works best if Jmin is much smaller than 5 (preferably zero)
    largeJs = linspace(5, Jmax, deltaN);   // 5+deltaN?
    fraction = floor((smallN-1)/4.0);                   // Unit to check progress. Flooring just in case.
    for(int i=0; i<deltaN; i++){
        if (i==fraction)                cout << "small 25% done, time:" << clock()-a << endl;
        if (i==fraction*2)              cout << "small 50% done, time:" << clock()-a << endl;
        if (i==fraction*3)              cout << "small 75% done, time:" << clock()-a << endl;
        //Done with printing
        J = smallJs[i];
        runit(h1, h2, J);
        coordinateFile << J << " " << Iteration.delta_rho_Epp << " " << Iteration.delta_rho_Emm << " " << Iteration.delta_rho_2nd << " " <<Iteration.delta_rho_3rd << " " << endl;}
    fraction = floor((deltaN-1)/4.0);                   // Unit to check progress. Flooring just in case.
    for(int i=0; i<deltaN; i++){
        if (i==fraction)                cout << "large 25% done, time:" << clock()-a << endl;
        if (i==fraction*2)              cout << "large 50% done, time:" << clock()-a << endl;
        if (i==fraction*3)              cout << "large 75% done, time:" << clock()-a << endl;
        //Done with printing
        J = largeJs[i];
        runit(h1, h2, J);
        coordinateFile << J << " " << Iteration.delta_rho_Epp << " " << Iteration.delta_rho_Emm << " " << Iteration.delta_rho_2nd << " " << Iteration.delta_rho_3rd << " " << endl;}
    coordinateFile.close();
    }


void PH_Running::plot_rangeJ(double h1, double h2, double Jmin, double Jmax)
{
    if((Jmax > 5) and (Jmin < 5))        plot_rangeJ_largeJ(h1, h2, Jmin, Jmax);
    else                                 plot_rangeJ_normal(h1, h2, Jmin, Jmax);
}

*/


