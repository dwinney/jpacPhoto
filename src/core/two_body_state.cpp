// Base class to work with four vectors.
//
// All code can be easily modified to use TLorentzVector instead
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "two_body_state.hpp"

// ---------------------------------------------------------------------------
// The four momentum of the meson
std::complex<double> jpacPhoto::two_body_state::q(int mu, double s, double theta)
{
    switch (mu)
    {
        case 0: return meson_energy(s);
        case 1: return momentum(s) * sin(theta);
        case 2: return 0.;
        case 3: return momentum(s) * cos(theta);
        default: 
        {
            std::cout << "two_body_state: Invalid four vector component! \n";
            return 0.;
        }
    }
};

// ---------------------------------------------------------------------------
// The four momenta of the baryon
std::complex<double> jpacPhoto::two_body_state::p(int mu, double s, double theta)
{
    switch (mu)
    {
        case 0: return baryon_energy(s);
        case 1: return - momentum(s) * sin(theta);
        case 2: return 0.;
        case 3: return - momentum(s) * cos(theta);
        default: 
        {
            std::cout << "two_body_state: Invalid four vector component! \n";
            return 0.;
        }
    }
};
