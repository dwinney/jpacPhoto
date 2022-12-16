// Simplified version of the 2->2 box diagram with contact S-wave vertices
// i.e. S-wave angular behavior and energy depedence given by scalar bubble integral

// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "bubble_amplitude.hpp"

// ---------------------------------------------------------------------------
// Assume an S-wave contact interaction between the intial and final state

std::complex<double> jpacPhoto::bubble_amplitude::helicity_amplitude(std::array<int, 4> helicities, double s, double t)
{
    // Save energies and helicities
    update(helicities, s, t);

    // Contract indices
    std::complex<double> result = 0.;
    for (int mu = 0; mu < 4; mu++)
    {
        for (int nu = 0; nu < 4; nu++)
        {
            std::complex<double> temp;
            temp  = _covariants->beam_field_tensor(mu, nu);
            temp *= METRIC[mu];
            temp *= _covariants->meson_polarization(mu);
            temp *= METRIC[nu];
            temp *= baryon_current(nu);

            result += temp / s;
        }
    };

    std::complex<double> loop;
    loop = G(s, _m1, _m2);
    if (_twochannel) loop -= _r * G(s, _m3, _m4);

    return _norm  * loop * result;
};

// Alias function for the vector fermion current
std::complex<double> jpacPhoto::bubble_amplitude::baryon_current(int mu)
{
    // Contract dirac-space indices
    std::complex<double> result = 0.;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            std::complex<double> temp;
            temp  = _covariants->recoil_spinor(i);
            temp *= GAMMA[mu][i][j];
            temp *= _covariants->target_spinor(j);

            result += temp;
        }
    }

    return result;
};

// Scalar loop function
// Subtracted at G(th) = 0 
std::complex<double> jpacPhoto::bubble_amplitude::G0(double s, double m1, double m2)
{
    std::complex<double> rho, R, xi;
    std::complex<double> result;

    rho = sqrt( (s + IEPS) - (m1+m2)*(m1+m2))*sqrt( (s - IEPS) - (m1-m2)*(m1-m2)) / s;

    R   = (m1*m1 + m2*m2 - s*(1. + rho)) / (2.*m1*m2) * XR;

    xi  = 1. - (m1+m2)*(m1+m2)/s;

    result = - ( rho*log(R) - xi*(m2-m1)/(m2+m1)*log(m2/m1) ) / PI;

    return result;
};

std::complex<double> jpacPhoto::bubble_amplitude::G(double s, double m1, double m2)
{
    if (2.*_n < 0.5) 
    {
        std::cout << "bubble_amplitude:: Error! Suppression power n < 0.25! Integral diverges... returning 0." << std::endl;
        return 0.;
    };

    double sth     = (m1+m2)*(m1+m2);
    double stheff  = (m1+m2-_eps)*(m1+m2-_eps);

    auto rho = [&](double x)
    {
        return sqrt( x + IEPS - (m1+m2)*(m1+m2))*sqrt( x + IEPS - (m1-m2)*(m1-m2)) / x;
    };

    auto f   = [&](double sp)
    {
        return rho(sp) * pow( (sp - stheff), 1.-2.*_n);
    };

    std::complex<double> fs;
    (s > sth) ? (fs = f(s)) : (fs = 0.);

    auto g = [&] (double sp)
    {
        std::complex<double> fsp = f(sp);
        return (fsp - fs) / (sp - s - IEPS) / ((sp - stheff));
    };

    std::complex<double> result;
    result  = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(g, sth, std::numeric_limits<double>::infinity(), 0, 1.E-6, NULL);
    if (s > sth) result -= fs  / (s - stheff) * log( ( sth - XR * s - IEPS)/( sth - stheff ) );
    return result / PI;
};