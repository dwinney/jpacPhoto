// Container class that calculates specific helicity partial-wave projection
// of a given amplitude
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "helicity_PWA.hpp"

// ---------------------------------------------------------------------------
// Partial-wave projection on the s-channel

std::complex<double> jpacPhoto::helicity_PWA::eval(double s)
{
    if (s < _kinematics->sth()) return 0.;

    // Here we want to make sure that the helicities are defined with respect to the s-channel
    // easiest way to do that currently is to force the covariant evaluation
    _amplitude->force_covariant(true);
    if ( _amplitude->helicity_CM_frame() != S )
    {
        std::cout << "Careful! helicity_PWA() called on amplitude with helicities not defined in s-channel. Results may vary..." << std::endl;
    }

    // Quick check that J is large enough
    // Net helicities
    int lam  = 2 * _helicities[0] - _helicities[1]; // Photon - Target
    int lamp = 2 * _helicities[2] - _helicities[3]; // Meson  - Recoil
    if ( abs(lamp) > _J || abs(lamp) > _J ) return 0.; 

    // Then just calculate the PWA integral
    auto F = [&](double theta)
    {
        double t = _kinematics->t_man(s, theta);

        std::complex<double> integrand;
        integrand  = sin(theta);
        integrand *= wigner_d_half(_J, lam, lamp, theta);
        integrand *= _amplitude->helicity_amplitude(_helicities, s, t);
        return integrand;
    };
    
    std::complex<double> result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(F, 0., PI, 0., 1.E-6, NULL);
    return result / 2.;
};

// ---------------------------------------------------------------------------
// Output the helicity parial wave projection based on the saved interpolation
double jpacPhoto::helicity_PWA::interpolation(double s)
{
    if ( _smax == -1 ) std::cout << "helicity_PWA interpolation interval not set! Returning 0." << std::endl;
    if ( s < _smin || s > _smax ) return 0;
    if (!_interpSaved) update_interpolation();
    return _interp->Eval(s);
};

// Update the inteprolation thats saved
void jpacPhoto::helicity_PWA::update_interpolation()
{
    // Clear any existing interpolation
    if ( _interpSaved ) 
    {
        delete _interp; _interpSaved = false;
    }
    _x.clear(); _fx.clear();

    // Fill up vectors of s values and the pwa
    for (int i = 0; i < _Ninterp; i++)
    {
        double x  = _smin + (_smax - _smin) * double(i) / double(_Ninterp - 1);
        double fx = real_part(x);

        _x.push_back(x); _fx.push_back(fx);
    };  

    _interp = new ROOT::Math::Interpolator(_x, _fx, ROOT::Math::Interpolation::kCSPLINE);
};