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

    for (int l = 0; l <= _lmax; l++)
    {
        sum += double(2*l+1)* P_l(l, cos(_theta)) * f_l(l);
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
        case 5: return z*(63.*z*z*z*z - 70.*z*z + 15.)/8.;
        default:
        {
            std::cout << "scattering_length::Pl() : ell value " + std::to_string(l) + " not implemented! Returning 0." << std::endl;
            return 0.;
        }
    };

    return 0.;
};

// so that they can be appropriately analytically continued below threshold
std::complex<double> jpacPhoto::scattering_length::barrier_factor(int l, double m1, double m2)
{         
    std::complex<double> pq = sqrt( Kallen(_s * XR, _mB*_mB * XR, _mT*_mT * XR)*Kallen(_s * XR, m1*m1 * XR, m2*m2 * XR) ) / (4.*_s);           
    return pow( pq * XR, double(l) );
};

std::complex<double> jpacPhoto::scattering_length::rho(double m1, double m2, double s)
{
    return sqrt( Kallen(s * XR, m1*m1 * XR, m2*m2 * XR) ) / s;
};

std::complex<double> jpacPhoto::scattering_length::rhoCM(double m1, double m2, double s)
{
    std::complex<double> Rho, xi;
    std::complex<double> result;

    Rho = rho(m1, m2, s);
    xi  = 1. - (m1+m2)*(m1+m2)/s;

    result = ( Rho*log((xi + Rho) / (xi - Rho)) - xi*(m2-m1)/(m2+m1)*log(m2/m1) ) / PI;

    return result;
};