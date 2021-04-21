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
    _s = s; _t = t, _theta = _kinematics->theta_s(s, t);
    _u = _kinematics->u_man(s, _theta);

    std::complex<double> result = 0.;
    result = covariant_amplitude(helicities);
    
    result *= form_factor();

    return result;
};

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
std::complex<double> jpacPhoto::dirac_exchange::covariant_amplitude(std::array<int,4> helicities)
{
    int lam_gam = helicities[0];
    int lam_tar = helicities[1];
    int lam_vec = helicities[2];
    int lam_rec = helicities[3];

    std::complex<double> result = 0.;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            std::complex<double> temp;
            temp  = top_vertex(i, lam_gam, lam_rec);
            temp *= dirac_propagator(i, j);
            temp *= bottom_vertex(j, lam_vec, lam_tar);

            result += temp;
        }
    }
    
    return result;
};

//------------------------------------------------------------------------------
// Photon fermion fermion vertex
// (ubar epsilon-slashed)
std::complex<double> jpacPhoto::dirac_exchange::top_vertex(int i, int lam_gam, int lam_rec)
{
    std::complex<double> result = 0.;
    for (int k = 0; k < 4; k++)
    {
        std::complex<double> temp;
        temp  = _kinematics->_recoil->adjoint_component(k, lam_rec, _s, _theta + PI); // theta_recoil = theta + pi
        temp *= slashed_eps(k, i, lam_gam, _kinematics->_eps_gamma, false, _s, 0.); // theta_gamma = 0

        result += temp;
    }

    return _gGam * result;
};

//------------------------------------------------------------------------------
// Vector fermion fermion vertex
// (epsilon*-slashed u)
std::complex<double> jpacPhoto::dirac_exchange::bottom_vertex(int j, int lam_vec, int lam_tar)
{
    std::complex<double> result = 0.;
    
    // F - F - V coupling
    if (_kinematics->_jp == VECTOR)
    {
        for (int k = 0; k < 4; k++)
        {
            std::complex<double> temp;
            temp  = XI * slashed_eps(j, k, lam_vec, _kinematics->_eps_vec, true, _s, _theta); //theta_vec = theta
            temp *= _kinematics->_target->component(k, lam_tar, _s, PI); // theta_target = pi

            result += temp;
        }
    }

    // F - F - P coupling
    else if (_kinematics->_jp == PSEUDO_SCALAR)
    {
        for (int k = 0; k < 4; k++)
        {
            std::complex<double> temp;
            temp  = XI * GAMMA_5[j][k];
            temp *= _kinematics->_target->component(k, lam_tar, _s, PI); // theta_target = pi

            result += temp;
        }
    }

    return _gVec * result;
};

//------------------------------------------------------------------------------
double jpacPhoto::dirac_exchange::exchange_mass()
{
    double result = 0.;
    for (int mu = 0; mu < 4; mu++)
    {
        std::complex<double> temp;
        temp  = _kinematics->u_exchange_momentum(mu, _s, _theta);
        temp *= METRIC[mu];
        temp *= _kinematics->u_exchange_momentum(mu, _s, _theta);

        result += real(temp);
    }

    return result;
}

std::complex<double> jpacPhoto::dirac_exchange::slashed_exchange_momentum(int i, int j)
{
    std::complex<double> result = 0.;
    for (int mu = 0; mu < 4; mu++)
    {
        std::complex<double> temp;
        temp  = GAMMA[mu][i][j];
        temp *= METRIC[mu];
        temp *= _kinematics->u_exchange_momentum(mu, _s, _theta);

        result += temp;
    }

    return result;
};

//------------------------------------------------------------------------------
// Slashed polarization vectors
std::complex<double> jpacPhoto::dirac_exchange::slashed_eps(int i, int j, double lam, polarization_vector * eps, bool STAR, double s, double theta)
{
    std::complex<double> result = 0.;
    for (int mu = 0; mu < 4; mu++)
    {
        std::complex<double> temp;
        if (STAR == false)
        {
            temp  = eps->component(mu, lam, s, theta);
        }
        else
        {
            temp = eps->conjugate_component(mu, lam, s, theta);
        }
        temp *= METRIC[mu];
        temp *= GAMMA[mu][i][j];

        result += temp;
    }

    return result;
};


//------------------------------------------------------------------------------
std::complex<double> jpacPhoto::dirac_exchange::dirac_propagator(int i, int j)
{
    std::complex<double> result;
    result = slashed_exchange_momentum(i, j);

    if (i == j)
    {
        result += _mEx;
    }

    result /= _u - _mEx2;

    return result;
};
