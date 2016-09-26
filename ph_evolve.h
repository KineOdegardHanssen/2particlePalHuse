#ifndef PH_EVOLVE_H
#define PH_EVOLVE_H
#include <ph_system.h>

class PH_Evolve
{
public:

    double delta_rho_Epp, delta_rho_Emm, delta_rho_2nd, delta_rho_3rd, rhoQpp, rhoQmm, rhoQ2nd, rhoQ3rd,
           Eppe, Emme, lambda1e, lambda2e, minele, maxele, maxit, tolerance, nm_diff;

    PH_System System;


    //Initialization
    PH_Evolve();
    PH_Evolve(int maxit, double sendinh1, double sendinh2, double sendinJ, double tolerance);



    double condition(double eigenenergy, double beta);
    double condition_derivative(double eigenenergy, double beta);
    double condition_Eppmin(double eigenenergy, double beta);
    double condition_Emmmin(double eigenenergy,double beta);
    double condition_lambda1min(double eigenenergy,double beta);
    double condition_lambda2min(double eigenenergy,double beta);
    double condition_derivative_Eppmin(double eigenenergy, double beta);
    double condition_derivative_Emmmin(double eigenenergy, double beta);
    double condition_derivative_lambda1min(double eigenenergy, double beta);
    double condition_derivative_lambda2min(double eigenenergy, double beta);

    double rhoth();
    double rhoth_Jis0();
    double rhoth_Eppmin();
    double rhoth_Emmmin();
    double rhoth_lambda1min();
    double rhoth_lambda2min();
    double rhoth_Jis0_Eppmin();
    double rhoth_Jis0_Emmmin();
    double rhoth_Jis0_lambda1min();
    double rhoth_Jis0_lambda2min();

    void newtonsmethod(double eigenenergy);

    void calculate_delta_rhos();
    void calculate_delta_rhos_Jis0();
};

#endif // PH_EVOLVE_H