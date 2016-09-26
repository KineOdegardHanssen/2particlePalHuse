#include "ph_system.h"
#include <iostream>
#include <iomanip>
#include <armadillo>
#include <cmath>

using namespace std;

PH_System::PH_System()
{

}

PH_System::PH_System(double senditinh1, double senditinh2, double senditinJ)
{
    this->h1 = senditinh1;
    this->h2 = senditinh2;
    this->J  =  senditinJ;
    this->beta=2;  // default value

    //cout << "h1: " << senditinh1 << ", h2:  "<< senditinh2 << ", J: "<< senditinJ << endl; // Seems fine

    if(J!=0)
    {
        constructHamiltonian();
        weights();
        Jis0 = false;
    }
    else
    {
        constructHamiltonian_Jis0();
        wp1 = 1; wp2 = 0;
        Jis0 = true;
    }

    //cout << "Epp = " << Epp << "; Emm = " << Emm << "; lambda1 = " << lambda1 << "; lambda2 = " << lambda2 << endl;

}

    void PH_System::constructHamiltonian()
    {
        this->Epp = (h1+h2+J)/2.0;
        this->Emm = (-h1-h2+J)/2.0;
        a1 = (h1-h2-J)/2.0;
        b1 = (-h1+h2-J)/2.0;
        c1 = J;
        this->lambda1 =(a1+b1+sqrt((a1+b1)*(a1+b1) - 4.0*(a1*b1-c1*c1)))/2.0;
        this->lambda2 =(a1+b1-sqrt((a1+b1)*(a1+b1) - 4.0*(a1*b1-c1*c1)))/2.0;

        elements = zeros(4);
        elements(0) = Epp;
        elements(1) = Emm;
        elements(2) = lambda1;
        elements(3) = lambda2;
        this->maxel = max(elements);
        this->minel = min(elements);
    }


    void PH_System::constructHamiltonian_Jis0()
    {
        this->Epp = (h1+h2)/2.0;
        this->Emm = -Epp;
        this->lambda1 = (h1-h2)/2.0;
        this->lambda2 = -lambda1;

        elements = zeros(4);
        elements(0) = Epp;
        elements(1) = Emm;
        elements(2) = lambda1;
        elements(3) = lambda2;
        this->maxel = max(elements);
        this->minel = min(elements);
    }


    void PH_System::weights()
    {
        va1y = -((1.0*a1-lambda1)/c1);   // First: Finding the relevant elements of the basis vectors
        va2x = -((1.0*b1-lambda2)/c1);

        this->walpha1 = 1.0/(1.0+va1y*va1y); // Weight of |++><++| in en.eig.state with E=lambda1

        this->walpha2 = va2x*va2x/(1.0+va2x*va2x); // Weight of |++><++| in en.eig.state with E=lambda2

        double r2 = (c1-va1y*a1)/(1.0-va1y*va2x);
        double r1 = a1 - va2x*r2;
        this->wp1 = r1/lambda1;
        this->wp2 = r2*va2x/lambda2;
    }
