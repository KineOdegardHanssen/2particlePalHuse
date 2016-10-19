#include "ph_evolve.h"
#include <typeinfo>
#include <iostream>
#include <iomanip>
#include <ph_system.h>


PH_Evolve::PH_Evolve()
{

}

PH_Evolve::PH_Evolve(int maxit, double sendinh1, double sendinh2, double sendinJ, double tolerance)
{
    this->tolerance=tolerance;
    this->maxit=maxit;
    finitetemp = false;

    //cout << "h1: " << sendinh1 << ", h2:  "<< sendinh2 << ", J: "<< sendinJ << endl; // Seems fine

    System = PH_System(sendinh1, sendinh2, sendinJ);

    Eppe = System.Epp;
    Emme = System.Emm;
    lambda1e = System.lambda1;
    lambda2e = System.lambda2;
    minele = System.minel;
    maxele = System.maxel;

    if(System.Jis0 == false)        calculate_delta_rhos();
    else                            calculate_delta_rhos_Jis0();
} // End initializer

PH_Evolve::PH_Evolve(int maxit, double sendinh1, double sendinh2, double sendinJ, double tolerance, bool finitetemp)
{
    this->tolerance=tolerance;
    this->maxit=maxit;
    this->finitetemp = finitetemp;

    //cout << "h1: " << sendinh1 << ", h2:  "<< sendinh2 << ", J: "<< sendinJ << endl; // Seems fine

    System = PH_System(sendinh1, sendinh2, sendinJ);

    Eppe = System.Epp;
    Emme = System.Emm;
    lambda1e = System.lambda1;
    lambda2e = System.lambda2;
    minele = System.minel;
    maxele = System.maxel;

    if(finitetemp)
    {
        if(System.Jis0 == false)
        {
            calculate_delta_rhos();
        }
        else
        {
            calculate_delta_rhos_Jis0();
        }
    }
    else    calculate_delta_rhos_infinitetemp();
}

void PH_Evolve::calculate_delta_rhos()  //should this be in running_plotting?
{
    rhoQpp = 1.0;                       //should these be here? It's practical, but...
    rhoQmm = 0.0;
    rhoQ2nd = System.walpha1;
    rhoQ3rd = System.walpha2;

    if((Eppe == maxele) || (Eppe == minele)){             // || is the or operator
        betausedpp = 0; // Default value when set like this
        this->delta_rho_Epp = 0.0;}//rhoQpp - 1 ;} // Extremal
    else{
        newtonsmethod(Eppe);
        betausedpp = System.beta; // Default value when set
        this->delta_rho_Epp = rhoQpp - rhoth();
        //cout << "E++ is not extremal. Delta_rho_Epp = " << delta_rho_Epp << endl;
    }
    if((Emme == maxele) || (Emme == minele)){
        betausedmm = 0; // Default value when set like this
        this->delta_rho_Emm = 0.0;}
    else{
        newtonsmethod(Emme);
        betausedpp = System.beta; // Default value when set like this
        this->delta_rho_Emm = rhoQmm - rhoth();}
    if((lambda1e == maxele) || (lambda1e == minele)){
        betausedl1 = 0; // Default value when set like this
        this->delta_rho_2nd = 0.0;}
    else{
        newtonsmethod(lambda1e);
        betausedl1 = System.beta; // Default value when set like this
        this->delta_rho_2nd = rhoQ2nd - rhoth();}
    if((lambda2e == maxele) || (lambda2e == minele)){
        betausedl2 = 0; // Default value when set like this
        this->delta_rho_3rd = 0.0;}
    else{
        newtonsmethod(lambda2e);
        betausedl2 = System.beta; // Default value when set like this
        this->delta_rho_3rd = rhoQ3rd - rhoth();}
}

void PH_Evolve::calculate_delta_rhos_Jis0()
{
    rhoQpp = 1; rhoQmm = 0; rhoQ2nd = 1; rhoQ3rd = 0;

    if((Eppe == maxele) || (Eppe == minele))       delta_rho_Epp = 0; //rhoQpp - 1; // Extremal
    else{
        newtonsmethod(Eppe);                                          //Call things to find beta
        this->delta_rho_Epp = rhoQpp - rhoth_Jis0();}
    if((Emme == maxele) || (Emme == minele))       delta_rho_Emm = rhoQmm;
    else{
        newtonsmethod(Emme);
        this->delta_rho_Emm = rhoQmm - rhoth_Jis0();}
    if((lambda1e == maxele) || (lambda1e == minele))   delta_rho_2nd = rhoQ2nd - 1;
    else{
        newtonsmethod(lambda1e);
        this->delta_rho_2nd = rhoQ2nd - rhoth_Jis0();}
    if((lambda2e == maxele) || (lambda2e == minele))   delta_rho_3rd = rhoQ3rd;
    else{
        newtonsmethod(lambda2e);
        this->delta_rho_3rd = rhoQ3rd - rhoth_Jis0();
      }
 }

void PH_Evolve::calculate_delta_rhos_infinitetemp()  //should this be in running_plotting?
{
    rhoQpp = 1.0;                       //should these be here? It's practical, but...
    rhoQmm = 0.0;
    rhoQ2nd = System.walpha1;
    rhoQ3rd = System.walpha2;
    System.beta = 0;          // Should beta really be a property of System?

    double rhothermal = 0.25*(1 + rhoQ2nd + rhoQ3rd);


    // The system is at infinite temperature, beta=0. No need to check which energy is the largest
    this->delta_rho_Epp = rhoQpp - rhothermal;
    this->delta_rho_Emm = rhoQmm - rhothermal;
    this->delta_rho_2nd = rhoQ2nd - rhothermal;
    this->delta_rho_3rd = rhoQ3rd - rhothermal;

    /*
    cout << "rhothermal: " << rhothermal << endl;
    cout << "rhoQ2nd = " << rhoQ2nd << "; rhoQ3rd = " << rhoQ3rd << endl;
    cout << "delta_rho_Epp = " << delta_rho_Epp << "; delta_rho_Emm = " << delta_rho_Emm << endl;
    cout << "delta_rho_2nd = " << delta_rho_2nd << "; delta_rho_3rd = " << delta_rho_3rd << endl;
    */

}


double PH_Evolve::condition(double eigenenergy, double beta)
{
    double conddifference = 0.0;
    if(Eppe == minele)            conddifference = condition_Eppmin(eigenenergy, beta);
    if(Emme == minele)            conddifference = condition_Emmmin(eigenenergy, beta);
    if(lambda1e == minele)        conddifference = condition_lambda1min(eigenenergy, beta);
    if(lambda2e == minele)        conddifference = condition_lambda2min(eigenenergy, beta);
    return conddifference;
 }

double PH_Evolve::condition_derivative(double eigenenergy, double beta)
{
    double condderdifference = 0.0;
    if(Eppe == minele)            condderdifference = condition_derivative_Eppmin(eigenenergy, beta);
    if(Emme == minele)            condderdifference = condition_derivative_Emmmin(eigenenergy, beta);
    if(lambda1e == minele)        condderdifference = condition_derivative_lambda1min(eigenenergy, beta);
    if(lambda2e == minele)        condderdifference = condition_derivative_lambda2min(eigenenergy, beta);
    return condderdifference;
 }


double PH_Evolve::condition_Eppmin(double eigenenergy, double beta)
{
    double energyaverage, Z;
    energyaverage = lambda1e*exp(-(lambda1e-Eppe)*beta)+lambda2e*exp(-(lambda2e-Eppe)*beta)+Eppe+Emme*exp(-(Emme-Eppe)*beta);
    Z =  1 + exp(-(lambda1e-Eppe)*beta) + exp(-(lambda2e-Eppe)*beta) + exp(-(Emme-Eppe)*beta);
    return energyaverage - Z*eigenenergy;
}

double PH_Evolve::condition_derivative_Eppmin(double eigenenergy, double beta)
{
    double energyaverage_der, Z_der;
    energyaverage_der = -lambda1e*lambda1e*exp(-(lambda1e-Eppe)*beta)-lambda2e*lambda2e*exp(-(lambda2e-Eppe)*beta)-Eppe*Eppe-Emme*Emme*exp(-(Emme-Eppe)*beta);
    Z_der =  -Eppe -lambda1e*exp(-(lambda1e-Eppe)*beta) -lambda2e*exp(-(lambda2e-Eppe)*beta) -Emme*exp(-(Emme-Eppe)*beta);
    return energyaverage_der - Z_der*eigenenergy;
}

double PH_Evolve::condition_Emmmin(double eigenenergy, double beta)
{
    double energyaverage, Z;
    energyaverage = lambda1e*exp(-(lambda1e-Emme)*beta)+lambda2e*exp(-(lambda2e-Emme)*beta)+Emme+Eppe*exp(-(Eppe-Emme)*beta);
    Z =  1 + exp(-(lambda1e-Emme)*beta) + exp(-(lambda2e-Emme)*beta) + exp(-(Eppe-Emme)*beta);
    return energyaverage - Z*eigenenergy;
}

double PH_Evolve::condition_derivative_Emmmin(double eigenenergy, double beta)
{
    double energyaverage_der, Z_der;
    energyaverage_der = -lambda1e*lambda1e*exp(-(lambda1e-Emme)*beta)-lambda2e*lambda2e*exp(-(lambda2e-Emme)*beta)-Emme*Emme-Eppe*Eppe*exp(-(Eppe-Emme)*beta);
    Z_der =  -Emme -lambda1e*exp(-(lambda1e-Emme)*beta) -lambda2e*exp(-(lambda2e-Emme)*beta) -Eppe*exp(-(Eppe-Emme)*beta);
    return energyaverage_der - Z_der*eigenenergy;
}

double PH_Evolve::condition_lambda1min(double eigenenergy, double beta)
{
    double energyaverage, Z;
    energyaverage = lambda1e+lambda2e*exp(-(lambda2e-lambda1e)*beta)+Emme*exp(-(Emme-lambda1e)*beta)+Eppe*exp(-(Eppe-lambda1e)*beta);
    Z =  1 + exp(-(lambda2e-lambda1e)*beta) + exp(-(Emme-lambda1e)*beta) + exp(-(Eppe-lambda1e)*beta);
    return energyaverage - Z*eigenenergy;
}

double PH_Evolve::condition_derivative_lambda1min(double eigenenergy, double beta)
{
    double energyaverage_der, Z_der;
    energyaverage_der = -Eppe*Eppe*exp(-(Eppe-lambda1e)*beta)-lambda2e*lambda2e*exp(-(lambda2e-lambda1e)*beta)-lambda1e*lambda1e-Emme*Emme*exp(-(Emme-lambda1e)*beta);
    Z_der =  -lambda1e -Emme*exp(-(Emme-lambda1e)*beta) -lambda2e*exp(-(lambda2e-lambda1e)*beta) -Eppe*exp(-(Eppe-lambda1e)*beta);
    return energyaverage_der - Z_der*eigenenergy;
}

double PH_Evolve::condition_lambda2min(double eigenenergy, double beta)
{
    double energyaverage, Z;
    energyaverage = lambda2e+lambda1e*exp(-(lambda1e-lambda2e)*beta)+Emme*exp(-(Emme-lambda2e)*beta)+Eppe*exp(-(Eppe-lambda2e)*beta);
    Z =  1 + exp(-(lambda1e-lambda2e)*beta) + exp(-(Emme-lambda2e)*beta) + exp(-(Eppe-lambda2e)*beta);
    return energyaverage - Z*eigenenergy;
}

double PH_Evolve::condition_derivative_lambda2min(double eigenenergy, double beta)
{
    double energyaverage_der, Z_der;
    energyaverage_der = -Eppe*Eppe*exp(-(Eppe-lambda2e)*beta)-lambda1e*lambda1e*exp(-(lambda1e-lambda2e)*beta)-lambda2e*lambda2e-Emme*Emme*exp(-(Emme-lambda2e)*beta);
    Z_der =  -lambda2e -Emme*exp(-(Emme-lambda2e)*beta) -lambda1e*exp(-(lambda1e-lambda2e)*beta) -Eppe*exp(-(Eppe-lambda2e)*beta);
    return energyaverage_der - Z_der*eigenenergy;
}






void PH_Evolve::newtonsmethod(double eigenenergy)
{
    double betaattempt = 2.0; // Arbitrary small value, hopefully close to our zero point.
    double diff = 1000.0;
    //double prevdiff = 1000.0;      // This implementation leads to an error for set values of h1, h2 :(
    //double prevprevdiff = 1000.0;  // I just have to live with a slow program
    int i = 0;
    double fbetan = 0.0;
    while((diff > tolerance) && (i < maxit)){
        fbetan = condition(eigenenergy, betaattempt);
        betaattempt -= fbetan/condition_derivative(eigenenergy, betaattempt);
        diff = abs(fbetan);
        i++;
        //if(i == maxit) cout<< "NB! Max no. of iterations exceeded in newtonsmethod. Diff="<< diff << endl;
        //prevprevdiff = prevdiff;
        //prevdiff = diff;
        //if(abs(prevprevdiff - diff) < 1e-15)  break; // To prevent the program from bouncing back and forth
    }
    System.beta = betaattempt;
    this->nm_diff = diff;      // Should try to use something like this to filter out bad betas...
}



double PH_Evolve::rhoth()
{
    double rhopp;
    if(Eppe == minele)            rhopp = rhoth_Eppmin();
    if(Emme == minele)            rhopp = rhoth_Emmmin();
    if(lambda1e == minele)        rhopp = rhoth_lambda1min();
    if(lambda2e == minele)        rhopp = rhoth_lambda2min();
    return rhopp;
}

double PH_Evolve::rhoth_Jis0()
{
    double rhoppJis0;
    if(Eppe == minele)            rhoppJis0 = rhoth_Jis0_Eppmin();
    if(Emme == minele)            rhoppJis0 = rhoth_Jis0_Emmmin();
    if(lambda1e == minele)        rhoppJis0 = rhoth_Jis0_lambda1min();
    if(lambda2e == minele)        rhoppJis0 = rhoth_Jis0_lambda2min();
    return rhoppJis0;
}

double PH_Evolve::rhoth_Eppmin()
{
    double beta = System.beta;
    double rhopp, Z;
    Z = 1 + exp(-(lambda1e-Eppe)*beta) + exp(-(lambda2e-Eppe)*beta) + exp(-(Emme-Eppe)*beta);
    rhopp = (1 + System.wp1*exp(-(lambda1e-Eppe)*beta) + System.wp2*exp(-(lambda2e-Eppe)*beta))/Z;
    return rhopp;
}

double PH_Evolve::rhoth_Jis0_Eppmin()
{
    double beta = System.beta;
    double rhopp, Z;
    Z = 1 + exp(-(lambda1e-Eppe)*beta) + exp(-(lambda2e-Eppe)*beta) + exp(-(Emme-Eppe)*beta);
    rhopp = (1 + exp(-(lambda1e-Eppe)*beta))/Z;
    return rhopp;
}

double PH_Evolve::rhoth_Emmmin()
{
    double beta = System.beta;
    double rhopp, Z;
    Z = exp(-(Eppe-Emme)*beta) + exp(-(lambda1e-Emme)*beta) + exp(-(lambda2e-Emme)*beta) + 1;
    rhopp = (exp(-(Eppe-Emme)*beta) + System.wp1*exp(-(lambda1e-Emme)*beta) + System.wp2*exp(-(lambda2e-Emme)*beta))/Z;
    return rhopp;
}

double PH_Evolve::rhoth_Jis0_Emmmin()
{
    double beta = System.beta;
    double rhopp, Z;
    Z = exp(-(Eppe-Emme)*beta) + exp(-(lambda1e-Emme)*beta) + exp(-(lambda2e-Emme)*beta) + 1;
    rhopp = (exp(-(Eppe-Emme)*beta) + exp(-(lambda1e-Emme)*beta))/Z;
    return rhopp;
}

double PH_Evolve::rhoth_lambda1min()
{
    double beta = System.beta;
    double rhopp, Z;
    Z = exp(-(Eppe-lambda1e)*beta) + 1 + exp(-(lambda2e-lambda1e)*beta) + exp(-(Emme-lambda1e)*beta);
    rhopp = (exp(-(Eppe-lambda1e)*beta) + System.wp1 + System.wp2*exp(-(lambda2e-lambda1e)*beta))/Z;
    return rhopp;
}

double PH_Evolve::rhoth_Jis0_lambda1min()
{
    double beta = System.beta;
    double rhopp, Z;
    Z = exp(-(Eppe-lambda1e)*beta) + 1 + exp(-(lambda2e-lambda1e)*beta) + exp(-(Emme-lambda1e)*beta);
    rhopp = (exp(-(Eppe-lambda1e)*beta) + 1)/Z;
    return rhopp;
}

double PH_Evolve::rhoth_lambda2min()
{
    double beta = System.beta;
    double rhopp, Z;
    Z = exp(-(Eppe-lambda2e)*beta) + 1 + exp(-(lambda1e-lambda2e)*beta) + exp(-(Emme-lambda2e)*beta);
    rhopp = (exp(-(Eppe-lambda2e)*beta) + System.wp1*exp(-(lambda1e-lambda2e)*beta) + System.wp2)/Z;
    return rhopp;
}

double PH_Evolve::rhoth_Jis0_lambda2min()
{
    double beta = System.beta;
    double rhopp, Z;
    Z = exp(-(Eppe-lambda2e)*beta) + 1 + exp(-(lambda1e-lambda2e)*beta) + exp(-(Emme-lambda2e)*beta);
    rhopp = (exp(-(Eppe-lambda2e)*beta) + exp(-(lambda1e-lambda2e)*beta))/Z;
    return rhopp;
 }


