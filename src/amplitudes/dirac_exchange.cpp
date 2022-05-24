// Spin-1/2 exchange ampltiude from perturbation theory
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/dirac_exchange.hpp"

//------------------------------------------------------------------------------
// Combine everything and contract indices
std::complex<double> jpacPhoto::dirac_exchange::helicity_amplitude(std::array<int, 4> helicities, double s, double t)
{
    // Store the invariant energies to avoid having to pass them around 
    update(helicities, s, t);

    // Exchange mass
    _u = _mB*_mB + _mT*_mT +_mX*_mX + _mR*_mR - _s - _t;

    // Multiply by form factor (this = 1 if no FF is specified)
    return form_factor() * covariant_amplitude();
};

// If we have a form factor this multipled above
double jpacPhoto::dirac_exchange::form_factor()
{
    switch (_useFF)
    {
        // exponential form factor
        case 1: 
        {
            return exp((_u - _kinematics->u_man(_s, 0.)) / _cutoff*_cutoff);
        };

        // monopole form factor
        case 2:
        {
            return (_cutoff*_cutoff - _mEx2) / (_cutoff*_cutoff - _u); 
        };

        default:
        {
            return 1.;
        };
    }
};

//------------------------------------------------------------------------------
// Contract spinor indices using covariant feynman rules
std::complex<double> jpacPhoto::dirac_exchange::covariant_amplitude()
{
    std::complex<double> result = 0.;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            std::complex<double> temp;
            temp  = top_vertex(i);
            temp *= dirac_propagator(i, j);
            temp *= bottom_vertex(j);

            result += temp;
        }
    }
    
    return result;
};

//------------------------------------------------------------------------------
// Beam -- Exchange -- Produced Baryon vertex

std::complex<double> jpacPhoto::dirac_exchange::top_vertex(int i)
{
    std::complex<double> result = 0.;
    for (int k = 0; k < 4; k++)
    {
        std::complex<double> temp;

        temp  = _covariants->recoil_spinor(k);
        temp *= _covariants->slashed_beam_momentum(k, i);

        result += temp;
    }

    return _gG * result;
};

//------------------------------------------------------------------------------
// Target -- Exchange -- Produced Meson vertex

std::complex<double> jpacPhoto::dirac_exchange::bottom_vertex(int j)
{
    std::array<int,2> JP = _kinematics->get_meson_JP();
    int jp = 10 * JP[0] + (1+JP[1])/2;
    
    switch (jp)
    {
        case 10: return vector_coupling(j);
        case  0: return pseudoscalar_coupling(j);
        default: return 0.;        
    }
};

std::complex<double> jpacPhoto::dirac_exchange::vector_coupling(int j)
{
    std::complex<double> result = 0.;

    for (int k = 0; k < 4; k++)
    {
        std::complex<double> temp;
        temp  = _covariants->slashed_meson_polarization(j, k);
        temp *= _covariants->target_spinor(k);

        result += temp;
    }

    return XI * _gN * result;
};

std::complex<double> jpacPhoto::dirac_exchange::pseudoscalar_coupling(int j)
{
    std::complex<double> result = 0.;

    for (int k = 0; k < 4; k++)
    {
        for (int l = 0; l < 4; l++)
        {
            std::complex<double> temp;
            temp  = GAMMA_5[j][k];
            temp *= _covariants->slashed_meson_momentum(k, l);
            temp *= _covariants->target_spinor(l);

            result += temp;
        }
    }

    return XI * _gN * result;
};


//------------------------------------------------------------------------------
std::complex<double> jpacPhoto::dirac_exchange::dirac_propagator(int i, int j)
{
    std::complex<double> result;
    result  = _covariants->slashed_u_momentum(i, j) + (i == j) * _mEx;
    result /= _u - _mEx2;

    return result;
};
