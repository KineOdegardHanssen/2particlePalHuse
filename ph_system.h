#ifndef PH_SYSTEM_H
#define PH_SYSTEM_H
#include <armadillo>
#include <random>

using namespace std;
using namespace arma;

class PH_System
{
public:

    double J, h1, h2, a1, b1, c1, Epp, Emm, lambda1, lambda2, beta, va1y, va2x, wp1, wp2, walpha1, walpha2,
           maxel, minel, vapx, vamy;

    bool Jis0;
    vec elements;

    //Initialization
    PH_System();
    PH_System(double senditinh1, double senditinh2, double senditinJ);

    void constructHamiltonian();
    void constructHamiltonian_Jis0();
    void weights();

};

#endif // PH_SYSTEM_H
