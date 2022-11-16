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
std::complex<double> jpacPhoto::scattering_length::barrier_factor(int l, double m1, double m2, double scale)
{         
    std::complex<double> pq = sqrt( Kallen(_s * XR, _mB*_mB * XR, _mT*_mT * XR)*Kallen(_s * XR, m1*m1 * XR, m2*m2 * XR) ) / (4.*_s);           
    return pow( pq / scale  * XR, l );
};

 std::complex<double> jpacPhoto::scattering_length::rho(int l)
{
    return barrier_factor(l) * sqrt( Kallen(_s * XR, _mX*_mX * XR, _mR*_mR * XR) ) / _s;
};

std::complex<double> jpacPhoto::scattering_length::rho_inelastic(int l)
{
    std::complex<double> result = 0.;
    for (int i = 0; i < _extra_thresholds.size(); i++)
    {
        double m1 = _extra_thresholds[i][0], m2 = _extra_thresholds[i][1];
        result   += _extra_couplings[i] * barrier_factor(l) * sqrt( Kallen(_s * XR, m1*m1 * XR, m2*m2 * XR) ) / _s;
    }

    return result;
};

// std::complex<double> jpacPhoto::scattering_length::chew_mandelstam(double m1, double m2)
// {
//     double sth     = (m1+m2)*(m1+m2);
//     double stheff  = (m1+m2-_eps)*(m1+m2-_eps);

//     auto rho = [&](double x)
//     {
//         return sqrt( x + IEPS - (m1+m2)*(m1+m2))*sqrt( x + IEPS - (m1-m2)*(m1-m2)) / x;
//     };

//     auto f   = [&](double sp)
//     {
//         return rho(sp) * pow( (sp - stheff), 1.-2.*_n);
//     };

//     std::complex<double> fs;
//     (s > sth) ? (fs = f(s)) : (fs = 0.);

//     auto g = [&] (double sp)
//     {
//         std::complex<double> fsp = f(sp);
//         return (fsp - fs) / (sp - s - IEPS) / ((sp - stheff));
//     };

//     std::complex<double> result;
//     result  = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(g, sth, std::numeric_limits<double>::infinity(), 0, 1.E-6, NULL);
//     if (s > sth) result -= fs  / (s - stheff) * log( ( sth - XR * s - IEPS)/( sth - stheff ) );
//     return result / PI;
// };