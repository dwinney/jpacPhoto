// Axial-vector meson photoproduction proceeding through the primakoff effect
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/primakoff_effect.hpp"

// ---------------------------------------------------------------------------
// Differential cross-sections with all the flux factors
double jpacPhoto::primakoff_effect::differential_xsection(double s, double t)
{
    // update saved energies
    update_kinematics(s, t); // calculate the other kinematics
    
    // output
    long double result = 1.;
    result  = ALPHA * _photonCoupling*_photonCoupling;
    result /= 8. * sqrt(_mA2) * _mX2 * _mX2 * _pGam * t*t;
    result /= (2. * sqrt(_mA2) * _nu - _mQ2);
    result *= W_00();

    // Amplitude depends on LT
    result *= amplitude_squared();
    
    // Convert from GeV^-2 -> nb
    result /= (2.56819E-6); 

    return result;
};

// ---------------------------------------------------------------------------
// Inegrated total cross-section
// IN NANOBARN
double jpacPhoto::primakoff_effect::integrated_xsection(double s)
{
  auto F = [&](double t)
  {
    return differential_xsection(s, t);
  };

  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS61);
  ROOT::Math::Functor1D wF(F);
  ig.SetFunction(wF);

  double t_min = _kinematics->t_man(s, 0.);
  double t_max = _kinematics->t_man(s, 1. * DEG2RAD); // Fall off is extremely fast in t so only integrate over that little bit

  return ig.Integral(t_max, t_min);
};

// ---------------------------------------------------------------------------
// Normalization
void jpacPhoto::primakoff_effect::calculate_norm()
{
    auto F = [&](double r)
    {
        return r * r * charge_distribution(r);
    };

    ROOT::Math::GSLIntegrator ig(   ROOT::Math::IntegrationOneDim::kADAPTIVE,
                                    ROOT::Math::Integration::kGAUSS61);

    ROOT::Math::Functor1D wF(F);
    ig.SetFunction(wF);

    // Integrate [0:inf]
    _rho0 = 1. / ig.IntegralUp(0.);
};

// ---------------------------------------------------------------------------
// Fourier transform the charge_distribution
double jpacPhoto::primakoff_effect::form_factor(double x)
{
    // momentum in the t channel
    double q = sqrt(x * (x - 4. * _mA2)) / (2. * sqrt(_mA2));

    auto dF = [&] (double r)
    {
        return r * sin(q * r) * charge_distribution(r);
    };

    ROOT::Math::GSLIntegrator ig(   ROOT::Math::IntegrationOneDim::kADAPTIVE,
                                    ROOT::Math::Integration::kGAUSS61);

    ROOT::Math::Functor1D wF(dF);
    ig.SetFunction(wF);
    
    return _rho0 * ig.IntegralUp(0.) / q;

};

// ---------------------------------------------------------------------------
// Amplitude
long double jpacPhoto::primakoff_effect::amplitude_squared()
{
    long double result;

    switch (_helProj)
    {
        // Longitudinal photon
        case 0: 
        {
            result = _pX*_pX * _mQ2 * _enX*_enX * _sinX2;
            break;
        };
        // Transverse photon
        case 1:
        {
            double coshalf2 = (1. + _cosX) / 2.;
            double sinhalf2 = (1. - _cosX) / 2.;

            double symC = _pX*_pGam*(_pX + _pGam) + _enX*_nu*(_pGam-_pX) - 2.*_pX*_pGam*_pGam*_cosX;
            double symS = _pX*_pGam*(_pX - _pGam) + _enX*_nu*(_pGam+_pX) - 2.*_pX*_pGam*_pGam*_cosX;

            double temp = pow(_pGam*(_nu*(_mX2+2.*_pX*_pX) - 2.*_enX*_pX*_pGam*_cosX), 2.) / (2. * _mX2);

            result  = pow(coshalf2 * symC, 2.);
            result += pow(sinhalf2 * symS, 2.);
            result += temp * _sinX2;

            break;
        };
        default: result = 0.;
    }

    return result;
};