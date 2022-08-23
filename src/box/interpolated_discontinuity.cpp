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
    if (phase) result *= _kinematics->parity_phase(hel_index, S);

    return result;
};

// ---------------------------------------------------------------------------
// Take in a function F which is the imaginary part of an amplitude and calculate the dispersion integral
std::complex<double> jpacPhoto::interpolated_discontinuity::dispersion(std::function<double(double)> f, double s)
{
    if ( _hardCutoff )
    {
        // To avoid numerical instabilities we take the principle value by explicitly subtracting the f(s) point and doing the residual integral analyitically
        double fs = f(s);

        auto g = [&] (double sp)
        {
            return (f(sp) - fs) / (sp - s - IEPS);
        };

        std::complex<double> intpiece, logpiece;
        intpiece = (1. / PI) * boost::math::quadrature::gauss_kronrod<double, 15>::integrate(g, _kinematics->sth() + EPS, _xi, 0, 1.E-6, NULL);
        logpiece = (fs / PI) * (log(_xi - s - 10.*IEPS) - log(_kinematics->sth() + EPS - s - 10.*IEPS));

        std::complex<double> result =  intpiece + logpiece;
    
        return result;
    }
    else
    {
        // To avoid numerical instabilities we take the principle value by explicitly subtracting the f(s) point and doing the residual integral analyitically
        double fs  = f(s);
        double sth = _kinematics->sth();

        auto g = [&] (double sp)
        {
            return ( f(sp) - fs ) / (sp * (sp - s - IEPS));
        };

        std::complex<double> intpiece, logpiece;
        intpiece =   ( s / PI) * boost::math::quadrature::gauss_kronrod<double, 15>::integrate(g, sth + EPS, 100., 0, 1.E-6, NULL);
        logpiece = - (fs / PI) * log(1. - s / sth - IEPS);
        
        std::complex<double> result =  _xi * (/* 1. + */ intpiece + logpiece);

        return result;
    }
};