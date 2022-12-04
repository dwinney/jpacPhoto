// Amplitude assuming a simple effective-range expansion near-threshold
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "scattering_length.hpp"

std::complex<double> jpacPhoto::scattering_length::helicity_amplitude( std::array<int,4> helicities, double s, double t)
{
    // This amplitude is helicity independent so we arbitrarily pick one helicty projection to make non-zero
    std::array<int,4> fixed_hels{1, 1, 1, 1};
    if (helicities != fixed_hels) return 0.;

    // Save inputs 
    update(helicities, s, t);

    // Sum over the partial-wave sum
    std::complex<double> sum = 0.;
    for (int l = 0; l <= _lmax; l++)
    {
        sum += double(2*l+1)* legendre(l, cos(_theta)) * A_L(l);
    };

    // Normalization here to get rid of helicity dependence in amplitude::probability_distribution
    // First the 2 removes the factor 1/4 when averaging over initial helicities
    // the 1/sqrt(2) removes the factor of 2 from the parity relation in amplitude::update_cache
    return 2. * sum / sqrt(2.);
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