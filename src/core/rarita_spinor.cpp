// Class for the spinor of a massive spin-3/2 particle
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "rarita_spinor.hpp"

// Construct the spin-3/2 wave-function from the stored polarization vector
// and dirac spinor. 
// NOTE: the vector is defined with respect to a particle in the +z direction
// while the spinor adds and explicit theta -> theta + pi because baryons are
// assumed to always be particle 2 
// To batch these we add an explicit M_PI to the angle in the vector components here.
std::complex<double> jpacPhoto::rarita_spinor::component(int i, int mu, int lambda, double s, double theta)
{
    // Double check masses havent changed
    sync_masses();

    double theta_1 = theta + M_PI;
    double theta_2 = theta;

    std::complex<double> v1p, v0, v1m, hp, hm;
    v1p =    _spin_one->component(mu, +1, s, theta_1);
    v0  =   -_spin_one->component(mu,  0, s, theta_1); // terms proportional to v0 get relative minus sign from Jacob-Wick phase
    v1m =    _spin_one->component(mu, -1, s, theta_1);

    hp  =    _spin_half->component(  i, +1, s, theta_2);
    hm  =    _spin_half->component(  i, -1, s, theta_2);

    switch(lambda)
    {
        case  3: return v1p * hp;  
        case  1: return sqrt(1./3.) * v1p * hm + sqrt(2. / 3.) * v0 * hp;
        case -1: return sqrt(1./3.) * v1m * hp + sqrt(2. / 3.) * v0 * hm;
        case -3: return v1m * hm;
        default: 
        {
            std::cout << "Error! Invalid helicity: " << lambda << " passed to rarita_spinor. Returning 0." << std::endl;
            return 0.;
        };
    };

    return 0.;
};
// Adjoint wave-function is similarly constructed from eps^* and bar u
std::complex<double> jpacPhoto::rarita_spinor::adjoint_component(int i, int mu, int lambda, double s, double theta)
{
    // Double check masses havent changed
    sync_masses();

    double theta_1 = theta + M_PI;
    double theta_2 = theta;

    std::complex<double> v1p, v0, v1m, hp, hm;
    v1p =    _spin_one->conjugate_component(mu, +1, s, theta_1);
    v0  =   -_spin_one->conjugate_component(mu,  0, s, theta_1); // terms proportional to v0 get relative minus sign from Jacob-Wick phase
    v1m =    _spin_one->conjugate_component(mu, -1, s, theta_1);

    hp  =    _spin_half->adjoint_component(  i, +1, s, theta_2);
    hm  =    _spin_half->adjoint_component(  i, -1, s, theta_2);

    switch(lambda)
    {
        case  3: return v1p * hp;  
        case  1: return sqrt(1./3.) * v1p * hm + sqrt(2. / 3.) * v0 * hp;
        case -1: return sqrt(1./3.) * v1m * hp + sqrt(2. / 3.) * v0 * hm;
        case -3: return v1m * hm;
        default: 
        {
            std::cout << "Error! Invalid helicity: " << lambda << " passed to rarita_spinor. Returning 0." << std::endl;
            return 0.;
        };
    };

    return 0.;
};
