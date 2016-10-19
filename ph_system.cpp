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
        Epp = (h1+h2+J)/2.0;
        Emm = (-h1-h2+J)/2.0;
        a1 = (h1-h2-J)/2.0;
        b1 = (-h1+h2-J)/2.0;
        c1 = J;
        lambda1 =(a1+b1+sqrt((a1-b1)*(a1-b1) + 4.0*c1*c1))/2.0;
        lambda2 =(a1+b1-sqrt((a1-b1)*(a1-b1) + 4.0*c1*c1))/2.0;

        elements = zeros(4);
        elements(0) = Epp;
        elements(1) = Emm;
        elements(2) = lambda1;
        elements(3) = lambda2;
        maxel = max(elements);
        minel = min(elements);

        double j=0;
        vec middles(2);
        for(int i=0; i<4;i++)
        {
            if(elements(i)!=maxel && elements(i)!=minel)
            {
                middles(j) = elements(i);
                j++;
            }
        }
        midelbig   =   max(middles);
        midelsmall = min(middles);

    }


    void PH_System::constructHamiltonian_Jis0()
    {
        Epp = (h1+h2)/2.0;
        Emm = -Epp;
        lambda1 = (h1-h2)/2.0;
        lambda2 = -lambda1;

        elements = zeros(4);
        elements(0) = Epp;
        elements(1) = Emm;
        elements(2) = lambda1;
        elements(3) = lambda2;
        maxel = max(elements);
        minel = min(elements);
    }


    void PH_System::weights()
    {
        double vapx_un = -((1.0*b1-lambda1)/c1);    // Unnormalized elements
        double vamy_un = -((1.0*a1-lambda2)/c1);
        double normpsq = 1.0/(1.0+vapx_un*vapx_un);   // To prevent round-off errors
        double normmsq = 1.0/(1.0+vamy_un*vamy_un);
        double vapxsq = vapx_un*vapx_un*normpsq;   // First: Finding the relevant elements of the basis vectors
        double vamxsq = normmsq;

        walpha1 = vapxsq; // Weight of |++><++| in en.eig.state with E=lambda1

        walpha2 = vamxsq; // Weight of |++><++| in en.eig.state with E=lambda2

        wp1 = walpha1;  // These are equal. Found analytically.
        wp2 = walpha2;  // These are equal. Found analytically.

        //cout << "v_{+x}: " << vapx_un << endl;
        //cout << "v_{-y}: " << vamy_un << endl;

        //cout << "E++ = " << Epp << "; E-- = " << Emm << endl;
        //cout << "lambda+ = " << lambda1 << "; lambda- = " << lambda2 << endl;
    }
