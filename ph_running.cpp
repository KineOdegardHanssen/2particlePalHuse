#include "ph_running.h"
#include <iostream>
#include <iomanip>
#include <armadillo>
#include <random>
#include <cmath>

using namespace std;
using namespace arma;

// So I guess I should do something like what I did in Python and plot J/|h|... Make an own function for that?
PH_Running::PH_Running()
{

}


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

void PH_Running::testit(double h1, double h2, double J)
{
    runit(h1, h2, J);
    cout << "|delta_rho| for:" << endl;
    cout << "The E++-state:" << abs(Iteration.delta_rho_Epp) << endl;
    cout << "The E---state:" << abs(Iteration.delta_rho_Emm) << endl;
    cout << "The lambda1-state:" << abs(Iteration.delta_rho_2nd) << endl;
    cout << "The lambda2-state:" << abs(Iteration.delta_rho_3rd) << endl;
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
    char *filename1 = new char[1000];                                     // File name can have max 1000 characters
    sprintf(filename1, "%s_coordinateFile.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    coordinateFile.open(filename1);
    delete filename1;

    // Open file for errors?
    //char *filename2 = new char[1000];                               // File name can have max 1000 characters
    //sprintf(filename2, "%s_diffFile.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    //diffFile.open(filename2);
    //delete filename2;

    // Open file for printing delta_rho_Epps to file, to investigate the divergence.
    char *filename2 = new char[1000];                                     // File name can have max 1000 characters
    sprintf(filename2, "%s_deltaEppFile.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    deltaEppFile.open(filename2);
    delete filename2;

    // Open file for printing e^(-beta*(Epp-evmin)) to file, to investigate the divergence.
    char *filename3 = new char[1000];                                     // File name can have max 1000 characters
    sprintf(filename3, "%s_EppexponentFile.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    Eppexponentfile.open(filename3);
    delete filename3;

    char *filename4 = new char[1000];                                     // File name can have max 1000 characters
    sprintf(filename4, "%s_betaFile.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    betaFile.open(filename4);
    delete filename4;

    char *filename5 = new char[1000];                                     // File name can have max 1000 characters
    sprintf(filename5, "%s_allthedataFile.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    allthedataFile.open(filename5);
    delete filename5;

    vec Jdivhssmall, Jdivhslarge, Jdivhs, hs;
    double h, counter;

    double minJ = J/hmax;
    double maxJ = J/hmin;
    cout << "minJ = " << minJ << "; maxJ = " << maxJ << endl;
    // Must have minJ<0.5 && maxJ > 0.6.
    int littleN = floor(0.9*N);
    int deltaN = N - littleN;
    Jdivhssmall = linspace(minJ, 0.65, littleN);
    Jdivhslarge = linspace(0.65, maxJ, deltaN);
    Jdivhs = zeros(N);
    hs = zeros(N);
    for(int j=0; j<littleN; j++)
    {
        hs[j] = J/Jdivhssmall[j];
        Jdivhs[j] = Jdivhssmall[j];
    }
    int j = littleN-1;
    for(int k=0; k<deltaN; k++)
    {
        j++;
        hs[j] = J/Jdivhslarge[k];
        Jdivhs[j] = Jdivhslarge[k];
        //cout << "j = " << j << "; hs[j] = " << hs[j] << endl;
    }
    // Could have more values for small Jdivhs, but that would take some time to make and result in ugly code...
    cout << max(Jdivhs) << endl;
    double rmsdelta_rho_Epp, rmsdelta_rho_Emm, rmsdelta_rho_2nd, rmsdelta_rho_3rd;
    double variance_delta_rho_Epp, variance_delta_rho_Emm, variance_delta_rho_2nd, variance_delta_rho_3rd;
    double stddelta_rho_Epp, stddelta_rho_Emm, stddelta_rho_2nd, stddelta_rho_3rd;
    double absmean_delta_rho_Epp, absmean_delta_rho_Emm, absmean_delta_rho_2nd, absmean_delta_rho_3rd;
    vec nos_of_Eppnans(N);
    for(int i=0; i<N; i++)    nos_of_Eppnans(i) = 0;
    for(int i=0; i<N; i++)
    {
        cout << "i= " << i << endl;
        deltaEppFile << " i = " << i << endl << endl;
        Eppexponentfile << "i = " << i << endl;
        betaFile << "i = " << i << endl;
        allthedataFile << "i = " << i << endl;
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
                /*
                if(isnan(Iteration.delta_rho_Epp)) // Only true if we have a nan.
                {
                    nos_of_Eppnans(i)++;
                    cout << "delta_rho_Epp is nan! i = " << i << endl;
                }
                */

                deltaEppFile << Iteration.delta_rho_Epp << " ";
                Eppexponentfile << exp(Iteration.betausedpp*(Iteration.minele-Iteration.Eppe)) << " ";
                betaFile << Iteration.betausedpp << " " << Iteration.betausedmm << " " << Iteration.betausedl1 << " " << Iteration.betausedl2 << endl;


                allthedataFile << "i = " << i << "; j = " << counter << ":" << endl;
                allthedataFile << "Betas, E++, E--, lambda+, lambda-: " << Iteration.betausedpp << " " << Iteration.betausedmm << " " << Iteration.betausedl1 << " " << Iteration.betausedl2 << endl;
                allthedataFile << "Eigenvalues, E++, E--, lambda+, lambda-: " << Iteration.Eppe << " " << Iteration.Emme << " " << Iteration.lambda1e << " " << Iteration.lambda2e << endl;
                allthedataFile << "Exponentials for Epp: (disregard if beta = 0): " << exp(Iteration.betausedpp*(Iteration.minele-Iteration.Eppe)) << " "<< exp(Iteration.betausedpp*(Iteration.minele-Iteration.Emme)) << " "<< exp(Iteration.betausedpp*(Iteration.minele-Iteration.lambda1e)) << " "<< exp(Iteration.betausedpp*(Iteration.minele-Iteration.lambda2e));
                allthedataFile << "Exponentials for Epp: (disregard if beta = 0): " << exp(Iteration.betausedmm*(Iteration.minele-Iteration.Eppe)) << " "<< exp(Iteration.betausedmm*(Iteration.minele-Iteration.Emme)) << " "<< exp(Iteration.betausedmm*(Iteration.minele-Iteration.lambda1e)) << " "<< exp(Iteration.betausedmm*(Iteration.minele-Iteration.lambda2e));
                allthedataFile << "Exponentials for Epp: (disregard if beta = 0): " << exp(Iteration.betausedl1*(Iteration.minele-Iteration.Eppe)) << " "<< exp(Iteration.betausedl2*(Iteration.minele-Iteration.Emme)) << " "<< exp(Iteration.betausedl1*(Iteration.minele-Iteration.lambda1e)) << " "<< exp(Iteration.betausedl1*(Iteration.minele-Iteration.lambda2e));
                allthedataFile << "Exponentials for Epp: (disregard if beta = 0): " << exp(Iteration.betausedl2*(Iteration.minele-Iteration.Eppe)) << " "<< exp(Iteration.betausedl1*(Iteration.minele-Iteration.Emme)) << " "<< exp(Iteration.betausedl2*(Iteration.minele-Iteration.lambda1e)) << " "<< exp(Iteration.betausedl2*(Iteration.minele-Iteration.lambda2e));
                allthedataFile << "Delta rhos, E++, E--, lambda+, lambda-: " <<  Iteration.delta_rho_Epp << " "<<  Iteration.delta_rho_Emm << " "<<  Iteration.delta_rho_2nd << " "<<  Iteration.delta_rho_3rd << " ";


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

        deltaEppFile << endl;
        deltaEppFile << "absmean: " << absmean_delta_rho_Epp << endl;
        deltaEppFile << "absmean, E++, E--, lambda+, lambda-: " << absmean_delta_rho_Epp << "" << absmean_delta_rho_Emm << " " << absmean_delta_rho_2nd << " "<< absmean_delta_rho_3rd << " "  << endl;

        //Print to file. This may be altered.
        // Print J/h instead? Or h/J?
        coordinateFile << Jdivhs[i] << " " << rmsdelta_rho_Epp << " " << rmsdelta_rho_Emm << " " << rmsdelta_rho_2nd << " " <<rmsdelta_rho_3rd << " " << variance_delta_rho_Epp << " " << variance_delta_rho_Emm << " " << variance_delta_rho_2nd << " " << variance_delta_rho_3rd << " " << stddelta_rho_Epp << " " << stddelta_rho_Emm << " " << stddelta_rho_2nd << " " << stddelta_rho_3rd << " " << absmean_delta_rho_Epp << " " << absmean_delta_rho_Emm << " " << absmean_delta_rho_2nd << " " << absmean_delta_rho_3rd << endl;
        // Must find a smarter way to print...

    } // End loop over i (values of h)
    coordinateFile.close();
    deltaEppFile.close();
    Eppexponentfile.close();
    betaFile.close();
    //diffFile.close();               // Will probably not use diffFile
    /*
    for(int i=0; i<N; i++)
    {
        cout << "i = " << i << "Number of NaNs in delta_rho_Epp: " <<nos_of_Eppnans(i) << endl;
    }
    */
}


void PH_Running::plot_randomuniform_Jdivh_pretty_infinitetemp(int averaging, double J, double hmin, double hmax)
{
    // Open file for printing delta_rhos to file
    char *filename1 = new char[1000];                                     // File name can have max 1000 characters
    sprintf(filename1, "%s_inftemp_coordinateFile.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    coordinateFile.open(filename1);
    delete filename1;

    // Open file for printing delta_rho_Epps to file, to investigate the divergence.
    /*
    char *filename2 = new char[1000];                                     // File name can have max 1000 characters
    sprintf(filename2, "%s_deltaEppFile.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    deltaEppFile.open(filename2);
    delete filename2;
    */

    bool finitetemp = false;


    vec Jdivhssmall, Jdivhslarge, Jdivhs, hs;
    double h, counter;

    double minJ = J/hmax;
    double maxJ = J/hmin;
    cout << "minJ = " << minJ << "; maxJ = " << maxJ << endl;
    // Must have minJ<0.5 && maxJ > 0.6.
    int littleN = floor(0.9*N);
    int deltaN = N - littleN;
    Jdivhssmall = linspace(minJ, 0.65, littleN);
    Jdivhslarge = linspace(0.65, maxJ, deltaN);
    Jdivhs = zeros(N);
    hs = zeros(N);
    for(int j=0; j<littleN; j++)
    {
        hs[j] = J/Jdivhssmall[j];
        Jdivhs[j] = Jdivhssmall[j];
    }
    int j = littleN-1;
    for(int k=0; k<deltaN; k++)
    {
        j++;
        hs[j] = J/Jdivhslarge[k];
        Jdivhs[j] = Jdivhslarge[k];
        //cout << "j = " << j << "; hs[j] = " << hs[j] << endl;
    }
    // Could have more values for small Jdivhs, but that would take some time to make and result in ugly code...
    cout << max(Jdivhs) << endl;
    double rmsdelta_rho_Epp, rmsdelta_rho_Emm, rmsdelta_rho_2nd, rmsdelta_rho_3rd;
    double variance_delta_rho_Epp, variance_delta_rho_Emm, variance_delta_rho_2nd, variance_delta_rho_3rd;
    double stddelta_rho_Epp, stddelta_rho_Emm, stddelta_rho_2nd, stddelta_rho_3rd;
    double absmean_delta_rho_Epp, absmean_delta_rho_Emm, absmean_delta_rho_2nd, absmean_delta_rho_3rd;
    vec nos_of_Eppnans(N);
    for(int i=0; i<N; i++)    nos_of_Eppnans(i) = 0;
    for(int i=0; i<N; i++)
    {
        cout << "i= " << i << endl;
        //deltaEppFile << " i = " << endl << endl;
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
            Iteration = PH_Evolve(maxit, h1, h2, J, tolerance, finitetemp);
            if(Iteration.nm_diff < 100*tolerance)
            {
                /*
                if(isnan(Iteration.delta_rho_Epp)) // Only true if we have a nan.
                {
                    nos_of_Eppnans(i)++;
                    cout << "delta_rho_Epp is nan! i = " << i << endl;
                }
                */

                //deltaEppFile << Iteration.delta_rho_Epp << " ";

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

        //deltaEppFile << endl;
        //deltaEppFile << "absmean: " << absmean_delta_rho_Epp << endl;

        //Print to file. This may be altered.
        // Print J/h instead? Or h/J?
        coordinateFile << Jdivhs[i] << " " << rmsdelta_rho_Epp << " " << rmsdelta_rho_Emm << " " << rmsdelta_rho_2nd << " " <<rmsdelta_rho_3rd << " " << variance_delta_rho_Epp << " " << variance_delta_rho_Emm << " " << variance_delta_rho_2nd << " " << variance_delta_rho_3rd << " " << stddelta_rho_Epp << " " << stddelta_rho_Emm << " " << stddelta_rho_2nd << " " << stddelta_rho_3rd << " " << absmean_delta_rho_Epp << " " << absmean_delta_rho_Emm << " " << absmean_delta_rho_2nd << " " << absmean_delta_rho_3rd << endl;
        // Must find a smarter way to print...

    } // End loop over i (values of h)
    coordinateFile.close();
}

void PH_Running::plot_hom_rangeh_infinitetemp(double J, double hmin, double hmax)
{
    double h, fraction;
    // Open coordinate file
    this->filenamePrefix = filenamePrefix;
    char *filename = new char[1000];                                    // File name can have max 1000 characters
    sprintf(filename, "%s_coordinateFile.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    coordinateFile.open(filename);
    delete filename;

    bool finitetemp = false;

    vec hs = zeros(N);

    hs = linspace(hmin, hmax, N);
    fraction = (N-1)/4.;
    for(int i=0; i<N; i++){
        h =  hs[i];
        Iteration = PH_Evolve(maxit, h, h, J, tolerance, finitetemp);
        coordinateFile << (J/h) << " " << abs(Iteration.delta_rho_Epp) << " " << abs(Iteration.delta_rho_Emm) << " " << abs(Iteration.delta_rho_2nd) << " " << abs(Iteration.delta_rho_3rd) << " " << endl;
      }
    coordinateFile.close();
}


void PH_Running::plot_asym_rangeh_infinitetemp(double J, double hmin, double hmax)
{
    double h, fraction;
    // Open coordinate file
    this->filenamePrefix = filenamePrefix;
    char *filename = new char[1000];                                    // File name can have max 1000 characters
    sprintf(filename, "%s_coordinateFile.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    coordinateFile.open(filename);
    delete filename;

    bool finitetemp = false;

    vec hs = zeros(N);

    hs = linspace(hmin, hmax, N);
    fraction = (N-1)/4.;
    for(int i=0; i<N; i++){
        h =  hs[i];
        Iteration = PH_Evolve(maxit, h, -h, J, tolerance, finitetemp);
        coordinateFile << (J/h) << " " << abs(Iteration.delta_rho_Epp) << " " << abs(Iteration.delta_rho_Emm) << " " << abs(Iteration.delta_rho_2nd) << " " << abs(Iteration.delta_rho_3rd) << " " << endl;
      }
    coordinateFile.close();
}

void PH_Running::plot_lopsided_rangeh_infinitetemp(double J, double hmin, double hmax)
{
    double h, fraction;
    // Open coordinate file
    this->filenamePrefix = filenamePrefix;
    char *filename = new char[1000];                                    // File name can have max 1000 characters
    sprintf(filename, "%s_coordinateFile.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    coordinateFile.open(filename);
    delete filename;

    bool finitetemp = false;

    vec hs = zeros(N);

    hs = linspace(hmin, hmax, N);
    fraction = (N-1)/4.;
    for(int i=0; i<N; i++){
        h =  hs[i];
        Iteration = PH_Evolve(maxit, h, -0.1*h, J, tolerance, finitetemp);
        coordinateFile << (J/h) << " " << abs(Iteration.delta_rho_Epp) << " " << abs(Iteration.delta_rho_Emm) << " " << abs(Iteration.delta_rho_2nd) << " " << abs(Iteration.delta_rho_3rd) << " " << endl;
      }
    coordinateFile.close();
}


void PH_Running::plot_randomuniform_Jdivh_pretty_infinitetemp_sorted(int averaging, double J, double hmin, double hmax)
{
    // Open file for printing delta_rhos to file
    char *filename1 = new char[1000];                                     // File name can have max 1000 characters
    sprintf(filename1, "%s_inftemp_coordinateFile.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    coordinateFile.open(filename1);
    delete filename1;

    // Open file for printing delta_rho_Epps to file, to investigate the divergence.
    /*
    char *filename2 = new char[1000];                                     // File name can have max 1000 characters
    sprintf(filename2, "%s_deltaEppFile.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    deltaEppFile.open(filename2);
    delete filename2;
    */

    bool finitetemp = false;


    vec Jdivhssmall, Jdivhslarge, Jdivhs, hs;
    double h, counter;

    double minJ = J/hmax;
    double maxJ = J/hmin;
    cout << "minJ = " << minJ << "; maxJ = " << maxJ << endl;
    // Must have minJ<0.5 && maxJ > 0.6.
    int littleN = floor(0.9*N);
    int deltaN = N - littleN;
    Jdivhssmall = linspace(minJ, 0.65, littleN);
    Jdivhslarge = linspace(0.65, maxJ, deltaN);
    Jdivhs = zeros(N);
    hs = zeros(N);
    for(int j=0; j<littleN; j++)
    {
        hs[j] = J/Jdivhssmall[j];
        Jdivhs[j] = Jdivhssmall[j];
    }
    int j = littleN-1;
    for(int k=0; k<deltaN; k++)
    {
        j++;
        hs[j] = J/Jdivhslarge[k];
        Jdivhs[j] = Jdivhslarge[k];
        //cout << "j = " << j << "; hs[j] = " << hs[j] << endl;
    }
    // Could have more values for small Jdivhs, but that would take some time to make and result in ugly code...
    cout << max(Jdivhs) << endl;
    double delta_rho_largest, delta_rho_2ndlargest, delta_rho_2ndsmallest, delta_rho_smallest;
    double rmsdelta_rho_largest, rmsdelta_rho_2ndlargest, rmsdelta_rho_2ndsmallest, rmsdelta_rho_smallest;
    double variance_delta_rho_largest, variance_delta_rho_2ndlargest, variance_delta_rho_2ndsmallest, variance_delta_rho_smallest;
    double stddelta_rho_largest, stddelta_rho_2ndlargest, stddelta_rho_2ndsmallest, stddelta_rho_smallest;
    double absmean_delta_rho_largest, absmean_delta_rho_2ndlargest, absmean_delta_rho_2ndsmallest, absmean_delta_rho_smallest;
    for(int i=0; i<N; i++)
    {
        cout << "i= " << i << endl;
        //deltaEppFile << " i = " << endl << endl;
        counter = 0;       // Should we fix the number of contributions?
        rmsdelta_rho_largest = 0;
        rmsdelta_rho_2ndlargest = 0;
        rmsdelta_rho_2ndsmallest = 0;
        rmsdelta_rho_smallest = 0;
        variance_delta_rho_largest = 0;
        variance_delta_rho_2ndlargest = 0;
        variance_delta_rho_2ndsmallest = 0;
        variance_delta_rho_smallest = 0;
        stddelta_rho_largest = 0;
        stddelta_rho_2ndlargest = 0;
        stddelta_rho_2ndsmallest = 0;
        stddelta_rho_smallest = 0;
        absmean_delta_rho_largest = 0;
        absmean_delta_rho_2ndlargest = 0;
        absmean_delta_rho_2ndsmallest = 0;
        absmean_delta_rho_smallest = 0;
        h = hs[i];
        std::default_random_engine generator;                       // I asked the internet, and it replied
        std::uniform_real_distribution<double> distribution(-h,h);
        while(counter<averaging)
        {
            double h1 = distribution(generator);
            double h2 = distribution(generator);
            Iteration = PH_Evolve(maxit, h1, h2, J, tolerance, finitetemp);
            if(Iteration.nm_diff < 100*tolerance)
            {
                // This is unelegant, but I do not bother making functions. That would be awkward.
                if(Iteration.Eppe==Iteration.maxele)               delta_rho_largest = Iteration.delta_rho_Epp;
                else if(Iteration.Emme==Iteration.maxele)          delta_rho_largest = Iteration.delta_rho_Emm;
                else if(Iteration.lambda1e==Iteration.maxele)      delta_rho_largest = Iteration.delta_rho_2nd;
                else if(Iteration.lambda2e==Iteration.maxele)      delta_rho_largest = Iteration.delta_rho_3rd;

                if(Iteration.Eppe==Iteration.minele)               delta_rho_smallest = Iteration.delta_rho_Epp;
                else if(Iteration.Emme==Iteration.minele)          delta_rho_smallest = Iteration.delta_rho_Emm;
                else if(Iteration.lambda1e==Iteration.minele)      delta_rho_smallest = Iteration.delta_rho_2nd;
                else if(Iteration.lambda2e==Iteration.minele)      delta_rho_smallest = Iteration.delta_rho_3rd;

                if(Iteration.Eppe==Iteration.midbigele)               delta_rho_2ndlargest = Iteration.delta_rho_Epp;
                else if(Iteration.Emme==Iteration.midbigele)          delta_rho_2ndlargest = Iteration.delta_rho_Emm;
                else if(Iteration.lambda1e==Iteration.midbigele)      delta_rho_2ndlargest = Iteration.delta_rho_2nd;
                else if(Iteration.lambda2e==Iteration.midbigele)      delta_rho_2ndlargest = Iteration.delta_rho_3rd;

                if(Iteration.Eppe==Iteration.midsmallele)               delta_rho_2ndsmallest = Iteration.delta_rho_Epp;
                else if(Iteration.Emme==Iteration.midsmallele)          delta_rho_2ndsmallest = Iteration.delta_rho_Emm;
                else if(Iteration.lambda1e==Iteration.midsmallele)      delta_rho_2ndsmallest = Iteration.delta_rho_2nd;
                else if(Iteration.lambda2e==Iteration.midsmallele)      delta_rho_2ndsmallest = Iteration.delta_rho_3rd;




                rmsdelta_rho_largest += delta_rho_largest*delta_rho_largest;
                rmsdelta_rho_2ndlargest += delta_rho_2ndlargest*delta_rho_2ndlargest;
                rmsdelta_rho_2ndsmallest += delta_rho_2ndsmallest*delta_rho_2ndsmallest;
                rmsdelta_rho_smallest += delta_rho_smallest*delta_rho_smallest;
                variance_delta_rho_largest += delta_rho_largest;
                variance_delta_rho_2ndlargest += delta_rho_2ndlargest;
                variance_delta_rho_2ndsmallest += delta_rho_2ndsmallest;
                variance_delta_rho_smallest += delta_rho_smallest;
                absmean_delta_rho_largest += abs(delta_rho_largest);
                absmean_delta_rho_2ndlargest += abs(delta_rho_2ndlargest);
                absmean_delta_rho_2ndsmallest += abs(delta_rho_2ndsmallest);
                absmean_delta_rho_smallest += abs(delta_rho_smallest);
                counter++;
            }
        }  // End while loop (samples of systems with h1,h2 from random uniform dist from [-h,h])

        // Variance (rmsdelta is still just the sum over the squares)
        // It's weird that I need to take the absolute value here. Something wrong?
        // But this is not really the variance? ... So what is it? is it useful at all?
        variance_delta_rho_largest = abs(rmsdelta_rho_largest - variance_delta_rho_largest*variance_delta_rho_largest)/averaging;
        variance_delta_rho_2ndlargest = abs(rmsdelta_rho_2ndlargest - variance_delta_rho_2ndlargest*variance_delta_rho_2ndlargest)/averaging;
        variance_delta_rho_2ndsmallest = abs(rmsdelta_rho_2ndsmallest - variance_delta_rho_2ndsmallest*variance_delta_rho_2ndsmallest)/averaging;
        variance_delta_rho_smallest = abs(rmsdelta_rho_smallest - variance_delta_rho_smallest*variance_delta_rho_smallest)/averaging;
        // Standard deviation
        // Not really the standard deviation?
        stddelta_rho_largest = sqrt(variance_delta_rho_largest);       // Could calculate this after reading from file
        stddelta_rho_2ndlargest = sqrt(variance_delta_rho_2ndlargest);
        stddelta_rho_2ndsmallest = sqrt(variance_delta_rho_2ndsmallest);
        stddelta_rho_smallest = sqrt(variance_delta_rho_smallest);
        // Root mean square
        rmsdelta_rho_largest = sqrt(rmsdelta_rho_largest/averaging);
        rmsdelta_rho_2ndlargest = sqrt(rmsdelta_rho_2ndlargest/averaging);
        rmsdelta_rho_2ndsmallest = sqrt(rmsdelta_rho_2ndsmallest/averaging);
        rmsdelta_rho_smallest = sqrt(rmsdelta_rho_smallest/averaging);
        // Mean
        absmean_delta_rho_largest = absmean_delta_rho_largest/averaging;
        absmean_delta_rho_2ndlargest = absmean_delta_rho_2ndlargest/averaging;
        absmean_delta_rho_2ndsmallest = absmean_delta_rho_2ndsmallest/averaging;
        absmean_delta_rho_smallest = absmean_delta_rho_smallest/averaging;

        //deltaEppFile << endl;
        //deltaEppFile << "absmean: " << absmean_delta_rho_Epp << endl;

        //Print to file. This may be altered.
        // Print J/h instead? Or h/J?
        coordinateFile << Jdivhs[i] << " " << rmsdelta_rho_largest << " " << rmsdelta_rho_2ndlargest << " " << rmsdelta_rho_2ndsmallest << " " <<rmsdelta_rho_smallest << " " << variance_delta_rho_largest << " " << variance_delta_rho_2ndlargest << " " << variance_delta_rho_2ndsmallest << " " << variance_delta_rho_smallest << " " << stddelta_rho_largest << " " << stddelta_rho_2ndlargest << " " << stddelta_rho_2ndsmallest << " " << stddelta_rho_smallest << " " << absmean_delta_rho_largest << " " << absmean_delta_rho_2ndlargest << " " << absmean_delta_rho_2ndsmallest << " " << absmean_delta_rho_smallest << endl;
        // Must find a smarter way to print...

    } // End loop over i (values of h)
    coordinateFile.close();
}

void PH_Running::plot_randomuniform_Jdivh_pretty_sorted(int averaging, double J, double hmin, double hmax)
{
    // Open file for printing delta_rhos to file
    char *filename1 = new char[1000];                                     // File name can have max 1000 characters
    sprintf(filename1, "%s_inftemp_coordinateFile.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    coordinateFile.open(filename1);
    delete filename1;

    // Open file for printing delta_rho_Epps to file, to investigate the divergence.
    /*
    char *filename2 = new char[1000];                                     // File name can have max 1000 characters
    sprintf(filename2, "%s_deltaEppFile.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
    deltaEppFile.open(filename2);
    delete filename2;
    */

    bool finitetemp = true;


    vec Jdivhssmall, Jdivhslarge, Jdivhs, hs;
    double h, counter;

    double minJ = J/hmax;
    double maxJ = J/hmin;
    cout << "minJ = " << minJ << "; maxJ = " << maxJ << endl;
    // Must have minJ<0.5 && maxJ > 0.6.
    int littleN = floor(0.9*N);
    int deltaN = N - littleN;
    Jdivhssmall = linspace(minJ, 0.65, littleN);
    Jdivhslarge = linspace(0.65, maxJ, deltaN);
    Jdivhs = zeros(N);
    hs = zeros(N);
    for(int j=0; j<littleN; j++)
    {
        hs[j] = J/Jdivhssmall[j];
        Jdivhs[j] = Jdivhssmall[j];
    }
    int j = littleN-1;
    for(int k=0; k<deltaN; k++)
    {
        j++;
        hs[j] = J/Jdivhslarge[k];
        Jdivhs[j] = Jdivhslarge[k];
        //cout << "j = " << j << "; hs[j] = " << hs[j] << endl;
    }
    // Could have more values for small Jdivhs, but that would take some time to make and result in ugly code...
    cout << max(Jdivhs) << endl;
    double delta_rho_largest, delta_rho_2ndlargest, delta_rho_2ndsmallest, delta_rho_smallest;
    double rmsdelta_rho_largest, rmsdelta_rho_2ndlargest, rmsdelta_rho_2ndsmallest, rmsdelta_rho_smallest;
    double variance_delta_rho_largest, variance_delta_rho_2ndlargest, variance_delta_rho_2ndsmallest, variance_delta_rho_smallest;
    double stddelta_rho_largest, stddelta_rho_2ndlargest, stddelta_rho_2ndsmallest, stddelta_rho_smallest;
    double absmean_delta_rho_largest, absmean_delta_rho_2ndlargest, absmean_delta_rho_2ndsmallest, absmean_delta_rho_smallest;
    for(int i=0; i<N; i++)
    {
        cout << "i= " << i << endl;
        //deltaEppFile << " i = " << endl << endl;
        counter = 0;       // Should we fix the number of contributions?
        rmsdelta_rho_largest = 0;
        rmsdelta_rho_2ndlargest = 0;
        rmsdelta_rho_2ndsmallest = 0;
        rmsdelta_rho_smallest = 0;
        variance_delta_rho_largest = 0;
        variance_delta_rho_2ndlargest = 0;
        variance_delta_rho_2ndsmallest = 0;
        variance_delta_rho_smallest = 0;
        stddelta_rho_largest = 0;
        stddelta_rho_2ndlargest = 0;
        stddelta_rho_2ndsmallest = 0;
        stddelta_rho_smallest = 0;
        absmean_delta_rho_largest = 0;
        absmean_delta_rho_2ndlargest = 0;
        absmean_delta_rho_2ndsmallest = 0;
        absmean_delta_rho_smallest = 0;
        h = hs[i];
        std::default_random_engine generator;                       // I asked the internet, and it replied
        std::uniform_real_distribution<double> distribution(-h,h);
        while(counter<averaging)
        {
            double h1 = distribution(generator);
            double h2 = distribution(generator);
            Iteration = PH_Evolve(maxit, h1, h2, J, tolerance, finitetemp);
            if(Iteration.nm_diff < 100*tolerance)
            {
                // This is unelegant, but I do not bother making functions. That would be awkward.
                if(Iteration.Eppe==Iteration.maxele)               delta_rho_largest = Iteration.delta_rho_Epp;
                else if(Iteration.Emme==Iteration.maxele)          delta_rho_largest = Iteration.delta_rho_Emm;
                else if(Iteration.lambda1e==Iteration.maxele)      delta_rho_largest = Iteration.delta_rho_2nd;
                else if(Iteration.lambda2e==Iteration.maxele)      delta_rho_largest = Iteration.delta_rho_3rd;

                if(Iteration.Eppe==Iteration.minele)               delta_rho_smallest = Iteration.delta_rho_Epp;
                else if(Iteration.Emme==Iteration.minele)          delta_rho_smallest = Iteration.delta_rho_Emm;
                else if(Iteration.lambda1e==Iteration.minele)      delta_rho_smallest = Iteration.delta_rho_2nd;
                else if(Iteration.lambda2e==Iteration.minele)      delta_rho_smallest = Iteration.delta_rho_3rd;

                if(Iteration.Eppe==Iteration.midbigele)               delta_rho_2ndlargest = Iteration.delta_rho_Epp;
                else if(Iteration.Emme==Iteration.midbigele)          delta_rho_2ndlargest = Iteration.delta_rho_Emm;
                else if(Iteration.lambda1e==Iteration.midbigele)      delta_rho_2ndlargest = Iteration.delta_rho_2nd;
                else if(Iteration.lambda2e==Iteration.midbigele)      delta_rho_2ndlargest = Iteration.delta_rho_3rd;

                if(Iteration.Eppe==Iteration.midsmallele)               delta_rho_2ndsmallest = Iteration.delta_rho_Epp;
                else if(Iteration.Emme==Iteration.midsmallele)          delta_rho_2ndsmallest = Iteration.delta_rho_Emm;
                else if(Iteration.lambda1e==Iteration.midsmallele)      delta_rho_2ndsmallest = Iteration.delta_rho_2nd;
                else if(Iteration.lambda2e==Iteration.midsmallele)      delta_rho_2ndsmallest = Iteration.delta_rho_3rd;


                rmsdelta_rho_largest += delta_rho_largest*delta_rho_largest;
                rmsdelta_rho_2ndlargest += delta_rho_2ndlargest*delta_rho_2ndlargest;
                rmsdelta_rho_2ndsmallest += delta_rho_2ndsmallest*delta_rho_2ndsmallest;
                rmsdelta_rho_smallest += delta_rho_smallest*delta_rho_smallest;
                variance_delta_rho_largest += delta_rho_largest;
                variance_delta_rho_2ndlargest += delta_rho_2ndlargest;
                variance_delta_rho_2ndsmallest += delta_rho_2ndsmallest;
                variance_delta_rho_smallest += delta_rho_smallest;
                absmean_delta_rho_largest += abs(delta_rho_largest);
                absmean_delta_rho_2ndlargest += abs(delta_rho_2ndlargest);
                absmean_delta_rho_2ndsmallest += abs(delta_rho_2ndsmallest);
                absmean_delta_rho_smallest += abs(delta_rho_smallest);
                counter++;
            }
        }  // End while loop (samples of systems with h1,h2 from random uniform dist from [-h,h])

        // Variance (rmsdelta is still just the sum over the squares)
        // It's weird that I need to take the absolute value here. Something wrong?
        // But this is not really the variance? ... So what is it? is it useful at all?
        variance_delta_rho_largest = abs(rmsdelta_rho_largest - variance_delta_rho_largest*variance_delta_rho_largest)/averaging;
        variance_delta_rho_2ndlargest = abs(rmsdelta_rho_2ndlargest - variance_delta_rho_2ndlargest*variance_delta_rho_2ndlargest)/averaging;
        variance_delta_rho_2ndsmallest = abs(rmsdelta_rho_2ndsmallest - variance_delta_rho_2ndsmallest*variance_delta_rho_2ndsmallest)/averaging;
        variance_delta_rho_smallest = abs(rmsdelta_rho_smallest - variance_delta_rho_smallest*variance_delta_rho_smallest)/averaging;
        // Standard deviation
        // Not really the standard deviation?
        stddelta_rho_largest = sqrt(variance_delta_rho_largest);       // Could calculate this after reading from file
        stddelta_rho_2ndlargest = sqrt(variance_delta_rho_2ndlargest);
        stddelta_rho_2ndsmallest = sqrt(variance_delta_rho_2ndsmallest);
        stddelta_rho_smallest = sqrt(variance_delta_rho_smallest);
        // Root mean square
        rmsdelta_rho_largest = sqrt(rmsdelta_rho_largest/averaging);
        rmsdelta_rho_2ndlargest = sqrt(rmsdelta_rho_2ndlargest/averaging);
        rmsdelta_rho_2ndsmallest = sqrt(rmsdelta_rho_2ndsmallest/averaging);
        rmsdelta_rho_smallest = sqrt(rmsdelta_rho_smallest/averaging);
        // Mean
        absmean_delta_rho_largest = absmean_delta_rho_largest/averaging;
        absmean_delta_rho_2ndlargest = absmean_delta_rho_2ndlargest/averaging;
        absmean_delta_rho_2ndsmallest = absmean_delta_rho_2ndsmallest/averaging;
        absmean_delta_rho_smallest = absmean_delta_rho_smallest/averaging;

        //deltaEppFile << endl;
        //deltaEppFile << "absmean: " << absmean_delta_rho_Epp << endl;

        //Print to file. This may be altered.
        // Print J/h instead? Or h/J?
        coordinateFile << Jdivhs[i] << " " << rmsdelta_rho_largest << " " << rmsdelta_rho_2ndlargest << " " << rmsdelta_rho_2ndsmallest << " " <<rmsdelta_rho_smallest << " " << variance_delta_rho_largest << " " << variance_delta_rho_2ndlargest << " " << variance_delta_rho_2ndsmallest << " " << variance_delta_rho_smallest << " " << stddelta_rho_largest << " " << stddelta_rho_2ndlargest << " " << stddelta_rho_2ndsmallest << " " << stddelta_rho_smallest << " " << absmean_delta_rho_largest << " " << absmean_delta_rho_2ndlargest << " " << absmean_delta_rho_2ndsmallest << " " << absmean_delta_rho_smallest << endl;
        // Must find a smarter way to print...

    } // End loop over i (values of h)
    coordinateFile.close();
}


