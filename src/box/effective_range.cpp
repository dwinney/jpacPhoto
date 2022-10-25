// Amplitude assuming a simple effective-range expansion near-threshold
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "effective_range.hpp"

std::complex<double> jpacPhoto::effective_range::helicity_amplitude( std::array<int,4> helicities, double s, double t)
{
    // Save inputs 
    update(helicities, s, t);

    // Sum over the partial-wave sum
    std::complex<double> sum = 0.;

    for (int ell = 0; ell <= _ellMax; ell++)
    {
        sum += double(2*ell+1)* P_l(ell, cos(_theta)) * f_l(ell);
    };

    return _N * sum / sqrt(6.);
};

double jpacPhoto::effective_range::P_l(int l, double z)
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
            std::cout << "effective_range::Pl() : ell value " + std::to_string(l) + " not implemented! Returning 0." << std::endl;
            return 0.;
        }
    };

    return 0.;
};

std::complex<double> jpacPhoto::effective_range::f_l(int ell)
{
    std::complex<double> K = (ell==0) ? (-1./_a + _r/2.*q2()) : (-1./_a);
    std::complex<double> denominator = 1./K - XI * (rho() * barrier_factor(ell));

    if ( ell == 0 )
    {
        for (int i = 0; i < _extra_couplings.size(); i++ )
        {
            denominator -= XI * _extra_couplings[i] * rho_inelastic(i) * _s;
        };
    };

    return barrier_factor(ell) / denominator;
};