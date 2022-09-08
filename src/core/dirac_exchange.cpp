// Spin-1/2 exchange ampltiude from perturbation theory
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "dirac_exchange.hpp"

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

        // Hyperon form factor 
        case 3:
        {
            return pow(_cutoff, 4.) / (pow(_cutoff, 4.) + pow(_u - _mEx2, 2.));
        };

        default: { return 1.; };
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
    std::array<int,2> JP = _kinematics->get_baryon_JP();
    int jp = 10 * JP[0] + (1+JP[1])/2;
    switch (jp)
    {
        case 11: return halfplus_coupling(i);
        case 10: return halfminus_coupling(i);
        case 30: return threehalfminus_coupling(i);
        default: return 0.;        
    };
}

std::complex<double> jpacPhoto::dirac_exchange::halfplus_coupling(int i)
{   
    std::complex<double> result = 0.;
    for (int k = 0; k < 4; k++)
    {
        std::complex<double> temp;

        temp  = _covariants->recoil_spinor(k);
        temp *= _covariants->slashed_beam_polarization(k, i);

        result += temp;
    }

    if ( _kinematics->is_photon() ) result *= -1.;
    return XI * _gG * result;
};

std::complex<double> jpacPhoto::dirac_exchange::halfminus_coupling(int i)
{   
    std::complex<double> result = 0.;
    for (int k = 0; k < 4; k++)
    {
        for (int l = 0; l < 4; l++)
        {
            std::complex<double> temp, sigmaF;

            temp  = _covariants->recoil_spinor(k);

            sigmaF  = _covariants->slashed_beam_momentum(k, l)     * _covariants->slashed_beam_polarization(l, i);
            sigmaF -= _covariants->slashed_beam_polarization(k, l) * _covariants->slashed_beam_momentum(l, i);

            temp *= sigmaF;

            result += temp;   
        }
    }

    return XI * _gG * result;
};

std::complex<double> jpacPhoto::dirac_exchange::threehalfminus_coupling(int i)
{   
    std::complex<double> result = 0.;
    for (int k = 0; k < 4; k++)
    {
        for (int mu = 0; mu < 4; mu++)
        {
            std::complex<double> temp;
            temp  = _covariants->recoil_spinor(k, mu);
            temp *= METRIC[mu];
            temp *= _covariants->beam_momentum(mu) * _covariants->slashed_beam_polarization(k, i)
                  - _covariants->slashed_beam_momentum(k, i) * _covariants->beam_polarization(mu);

            result += temp;
        }
    }

    return XI * _gG * result;
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

    return - XI * _gN * result;
};

std::complex<double> jpacPhoto::dirac_exchange::pseudoscalar_coupling(int j)
{
    std::complex<double> term1 = 0., term2 = 0.;

    if (abs(_gN) > 1.E-6 );
    {
        for (int k = 0; k < 4; k++)
        {
            std::complex<double> temp;
            temp  = GAMMA_5[j][k];
            temp *= _covariants->target_spinor(k);

            term1 += temp;
        };
    };

    
    if (abs(_gN2) > 1.E-6 );
    {
        for (int k = 0; k < 4; k++)
        {
            for (int l = 0; l < 4; l++)
            {
                std::complex<double> temp;
                temp  = GAMMA_5[j][k];
                temp *= _covariants->slashed_meson_momentum(k, l);
                temp *= _covariants->target_spinor(l);

                term2 += temp;
            };
        };
    };

    return (_gN * term1 + _gN2 * term2);
};


//------------------------------------------------------------------------------
std::complex<double> jpacPhoto::dirac_exchange::dirac_propagator(int i, int j)
{
    std::complex<double> result;
    result  = _covariants->slashed_u_momentum(i, j) - std::complex<double>(i == j) * _mEx;
    result /= _u - _mEx2;
    return - XI * result;
};
