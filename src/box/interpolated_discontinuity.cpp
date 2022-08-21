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
    // For fixed helicities we just multiply the corresponding HPWA's by d-functions and sum over J
    std::complex<double> sum = 0.;

    int lam  = 2 * helicities[0] - helicities[1]; // Photon - Target
    int lamp = 2 * helicities[2] - helicities[3]; // Meson  - Recoil

    double theta = _kinematics->theta_s(s, t);

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
std::complex<double> jpacPhoto::interpolated_discontinuity::helicity_pwa(int external_hel_index, int j, double s)
{
    if ((j < _jmax) || (s < _kinematics->sth())) return 0.;

    int  hel_index; 
    bool phase;
    if ( external_hel_index < _nAmps/2 ){ hel_index = external_hel_index;            phase = false; }
    else                                { hel_index = external_hel_index - _nAmps/2; phase = true;  }
    
    auto F = [&] (double x)
    {
        return _hpw_projections[j][hel_index]->eval(x, _eta);
    };

    std::complex<double> result = dispersion(F, s);
    if (phase) result *= _kinematics->parity_phase(hel_index, S);

    return result;
};

// ---------------------------------------------------------------------------
// Take in a function F which is the imaginary part of an amplitude and calculate the dispersion integral
std::complex<double> jpacPhoto::interpolated_discontinuity::dispersion(std::function<double(double)> f, double s)
{
    // To avoid numerical instabilities we take the principle value by explicitly subtracting the f(s) point and doing the residual integral analyitically
    double sub_point = f(s);

    auto g = [&] (double sp)
    {
        return f(sp) - sub_point;
    };

    std::complex<double> intpiece = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(g, _kinematics->sth() + EPS, _scut, 0, 1.E-6, NULL);
    std::complex<double> logpiece = sub_point * (log(_scut - s - IEPS) - log(_kinematics->sth() + EPS - s - IEPS));
    std::complex<double> result =  (intpiece + logpiece) / M_PI;

    return result;
};