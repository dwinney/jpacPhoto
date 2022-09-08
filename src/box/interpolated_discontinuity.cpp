// This discontinuity relied on PW projections of the disconinuity saved to file 
// as a interpolation_2D grid
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "interpolated_discontinuity.hpp"

// ---------------------------------------------------------------------------
// Evaluate the discontinutiy by summing partial-wave amplitudes with corresponding d-functions
std::complex<double> jpacPhoto::interpolated_discontinuity::helicity_amplitude(std::array<int,4> helicities, double s, double t)
{
    // Arguemnts of the d-function
    int lam  = 2 * helicities[0] - helicities[1]; // Photon - Target
    int lamp = 2 * helicities[2] - helicities[3]; // Meson  - Recoil

    double theta = _kinematics->theta_s(s, t);

    // For fixed helicities we just multiply the corresponding HPWA's by d-functions and sum over J
    std::complex<double> sum = 0.;
    for (int j = 0; (2*j+1) <= _jmax; j++)
    {
        int J = (2*j+1);
        if ( abs(lamp) > J || abs(lamp) > J ) continue; 
        sum += double(J+1) * wigner_d_half(J, lam, lamp, theta) * helicity_pwa(helicities, J, s);
    };

    return sum;
};

// ---------------------------------------------------------------------------
// Evaluate the discontinutiy by summing helicity amplitudes with corresponding d-functions
std::complex<double> jpacPhoto::interpolated_discontinuity::helicity_pwa(int external_hel_index, int J, double s)
{
    if ((J > _jmax) || (s < _kinematics->sth())) return 0.;
    int  j_index = (J - 1)/2;
    int  hel_index; 
    bool phase;
    if ( external_hel_index < _nAmps/2 ){ hel_index = external_hel_index;              phase = false; }
    else                                { hel_index = _nAmps - external_hel_index - 1; phase = true;  }
    
    auto F = [&] (double x)
    {
        return _hpw_projections[j_index][hel_index]->eval(x, _eta);
    };

    std::complex<double> result = dispersion(F, s);
    if (phase) result *= _kinematics->intrinsic_parity(S);

    return result;
};

// ---------------------------------------------------------------------------
// Take in a function F which is the imaginary part of an amplitude and calculate the dispersion integral
std::complex<double> jpacPhoto::interpolated_discontinuity::dispersion(std::function<double(double)> f, double s)
{
    if (_intermediateThreshold < 0)
    {
        std::cout << "interpolated_discontinutiy: Warning! intermediate threshold not set. Returning 0." << std::endl;
        return 0.;
    }
    if (_cutoff < _intermediateThreshold)
    {
        std::cout << "interpolated_discontinutiy: Warning! Cutoff below threshold. Returning 0." << std::endl;
        return 0.;
    }
    
    // To avoid numerical instabilities we take the principle value by explicitly subtracting the f(s) point and doing the residual integral analyitically
    double fs;
    (s > _intermediateThreshold) ? (fs = f(s)) : (fs = 0.);

    auto g = [&] (double sp)
    {
        return (f(sp) - fs) / (sp - s - IEPS);
    };

    std::complex<double> result;
    result  = (1. / PI) * boost::math::quadrature::gauss_kronrod<double, 31>::integrate(g, _intermediateThreshold, _cutoff, 0, 1.E-6, NULL);
    result += (fs / PI) * (log(_cutoff - XR * s) - log(_intermediateThreshold - XR * s));

    return result;
};