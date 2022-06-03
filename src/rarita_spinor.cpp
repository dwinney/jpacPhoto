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
    switch(lambda)
    {
        case  3: 
        {
            std::complex<double> v1, hp;

            v1 = _spin_one->component(mu, +1, s, theta + M_PI); 
            hp = _spin_half->component(i, +1, s, theta);

            return v1 * hp;
        };  
        case  1:
        {
            std::complex<double> v1, v0, hp, hm;

            v1 = _spin_one->component(mu, 1, s, theta + M_PI);
            v0 = _spin_one->component(mu, 0, s, theta + M_PI);
            hp = _spin_half->component(i, +1, s, theta);
            hm = _spin_half->component(i, -1, s, theta);

            return sqrt(1./3.) * v1 * hm + sqrt(2. / 3.) * v0 * hp;
        };
        case -1:
        {
            std::complex<double> v1, v0, hp, hm;

            v1 = _spin_one->component(mu, -1, s, theta + M_PI);
            v0 = _spin_one->component(mu,  0, s, theta + M_PI);
            hp = _spin_half->component(i, +1, s, theta);
            hm = _spin_half->component(i, -1, s, theta);

            return sqrt(1./3.) * v1 * hp + sqrt(2. / 3.) * v0 * hm;
        };
        case -3: 
        {
            std::complex<double> v1, hm;
            v1 = _spin_one->component(mu, -1, s, theta + M_PI); 
            hm = _spin_half->component(i, -1, s, theta);

            return v1 * hm;
        };  
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
    switch(lambda)
    {
        case  3: 
        {
            std::complex<double> v1, hp;

            v1 = _spin_one->conjugate_component(mu, +1, s, theta + M_PI); 
            hp = _spin_half->adjoint_component(i, +1, s, theta);

            return v1 * hp;
        };  
        case  1:
        {
            std::complex<double> v1, v0, hp, hm;

            v1 = _spin_one->conjugate_component(mu, 1, s, theta + M_PI);
            v0 = _spin_one->conjugate_component(mu, 0, s, theta + M_PI);
            hp = _spin_half->adjoint_component(i, +1, s, theta);
            hm = _spin_half->adjoint_component(i, -1, s, theta);

            return sqrt(1./3.) * v1 * hm + sqrt(2. / 3.) * v0 * hp;
        };
        case -1:
        {
            std::complex<double> v1, v0, hp, hm;

            v1 = _spin_one->conjugate_component(mu, -1, s, theta + M_PI);
            v0 = _spin_one->conjugate_component(mu,  0, s, theta + M_PI);
            hp = _spin_half->adjoint_component(i, +1, s, theta);
            hm = _spin_half->adjoint_component(i, -1, s, theta);

            return sqrt(1./3.) * v1 * hp + sqrt(2. / 3.) * v0 * hm;
        };
        case -3: 
        {
            std::complex<double> v1, hm;
            v1 = _spin_one->conjugate_component(mu, -1, s, theta + M_PI); 
            hm = _spin_half->adjoint_component(i, -1, s, theta);

            return v1 * hm;
        };  
        default: 
        {
            std::cout << "Error! Invalid helicity: " << lambda << " passed to rarita_spinor. Returning 0." << std::endl;
            return 0.;
        };
    };

    return 0.;
};
