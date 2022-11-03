// Amplitude assuming a simple effective-range expansion near-threshold
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "scattering_length.hpp"

std::complex<double> jpacPhoto::scattering_length::helicity_amplitude( std::array<int,4> helicities, double s, double t)
{
    // Save inputs 
    update(helicities, s, t);

    // Sum over the partial-wave sum
    std::complex<double> sum = 0.;

    for (int ell = 0; ell <= _ellMax; ell++)
    {
        sum +=  _N * double(2*ell+1)* P_l(ell, cos(_theta)) * f_l(ell);
    };

    return sum / sqrt(6.);
};

double jpacPhoto::scattering_length::P_l(int l, double z)
{
    switch (l) 
    {
        case 0: return 1.;
        case 1: return z;
        case 2: return 0.5*(3.*z*z - 1.);
        case 3: return 0.5*z*(5.*z*z - 3.);
        case 4: return (35.*z*z*z*z - 30.*z*z + 3.)/8.;
        case 5: return z*(63.*z*z*z*z - 70.*z*z + 15.*z)/8.;
        default:
        {
            std::cout << "scattering_length::Pl() : ell value " + std::to_string(l) + " not implemented! Returning 0." << std::endl;
            return 0.;
        }
    };

    return 0.;
};
