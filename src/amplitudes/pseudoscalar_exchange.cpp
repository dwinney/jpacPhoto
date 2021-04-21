// Charged axial-vector meson photoproduction proceeding through a pseudoscalar (pion) exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// References:
// [1] arXiv:1503.02125 [hep-ph]
// ---------------------------------------------------------------------------

#include "amplitudes/pseudoscalar_exchange.hpp"

//------------------------------------------------------------------------------
// Combine everything and contract indices
std::complex<double> jpacPhoto::pseudoscalar_exchange::helicity_amplitude(std::array<int, 4> helicities, double s, double t)
{
    int lam_gam = helicities[0];
    int lam_tar = helicities[1];
    int lam_vec = helicities[2];
    int lam_rec = helicities[3];

    // Store the invariant energies to avoid having to pass them around 
    _s = s; _t = t, _theta = _kinematics->theta_s(s, t);

    std::complex<double> result;

    if (_useCovariant == true || _debug == true)
    {
        // Because its a scalar exchange we dont have any loose indices to contract
        result  = top_vertex(lam_gam, lam_vec);
        result *= bottom_vertex(lam_tar, lam_rec);
        result *= scalar_propagator();
    }
    else
    {
        if (lam_vec != lam_gam || lam_tar != lam_rec) 
        {
            return 0.; 
        }
        else
        {
            result  = top_residue(lam_gam, lam_vec);
            result *= bottom_residue(lam_tar, lam_rec);
            result *= scalar_propagator();
        }
    }

    // add form factor if wanted
    result *= form_factor();    
    return result;
};

//------------------------------------------------------------------------------
// Whether to add a form factor to the propagator
double jpacPhoto::pseudoscalar_exchange::form_factor()
{
    switch (_useFormFactor)
    {
        // exponential form factor
        case 1: 
        {
            return exp((_t - _kinematics->t_man(_s, 0.)) / _cutoff*_cutoff);
        };
        // monopole form factor
        case 2:
        {
            return (_cutoff*_cutoff - _mEx2) / (_cutoff*_cutoff - _t); 
        };

        default:
        {
            return 1.;
        };
    }

    return 1.;
};

//------------------------------------------------------------------------------
// ANALYTIC RESIDUES

// Photon resiude 
std::complex<double> jpacPhoto::pseudoscalar_exchange::top_residue(int lam_gam, int lam_vec)
{
    // Spin-flip not allowed with spin-0 exchange
    if (lam_gam != lam_vec) return 0.;

    std::complex<double> result;

    if (_kinematics->_jp == AXIAL_VECTOR)
    {
        result  = 1. / _kinematics->_mX;
    }
    else if (_kinematics->_jp == VECTOR)
    {
        result = 1.;
        (_kinematics->_photon) ? (result *= -4.) : (result *= 1.);
    }

    else if (_kinematics->_jp == PSEUDO_SCALAR)
    {
        result  = 2. * XI / _kinematics->_mB;
    }

    std::complex<double> q_t = sqrt(XR * Kallen(_t, _kinematics->_mX2, _kinematics->_mB2)) / sqrt(4. * _t * XR);
    return 2. * _gGamma * q_t * sqrt(XR * _t) * result;
};

// Nucleon resiude 
std::complex<double> jpacPhoto::pseudoscalar_exchange::bottom_residue(int lam_tar, int lam_rec)
{
    std::complex<double> result;

    if (lam_tar != lam_rec) return 0.;

    result  = _gNN; 
    result *= sqrt(XR * _t - pow((_kinematics->_mT - _kinematics->_mR), 2.));
    result *= lam_tar / 2.;

    return result;
};

//------------------------------------------------------------------------------
// FEYNMAN EVALUATION

// Nucleon vertex
std::complex<double> jpacPhoto::pseudoscalar_exchange::bottom_vertex(int lam_tar, int lam_rec)
{
    std::complex<double> result = 0.;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            // ubar(recoil) * gamma_5 * u(target)
            std::complex<double> temp;
            temp  = _kinematics->_recoil->adjoint_component(i, lam_rec, _s, _theta + PI); // theta_recoil = theta + pi
            temp *= GAMMA_5[i][j];
            temp *= _kinematics->_target->component(j, lam_tar, _s, PI); // theta_target = pi

            result += temp;
        }
    }

    result *= _gNN;

    return result;
};

// Photon vertex
std::complex<double> jpacPhoto::pseudoscalar_exchange::top_vertex(int lam_gam, int lam_vec)
{
    std::complex<double> result = 0.;
    
    // A - V - P
    if (_kinematics->_jp == AXIAL_VECTOR)
    {
         std::complex<double> term1 = 0., term2 = 0.;
        for (int mu = 0; mu < 4; mu++)
        {
            for (int nu = 0; nu < 4; nu++)
            {
                // (eps*_lam . eps_gam)(q_vec . q_gam)
                std::complex<double> temp1;
                temp1  = _kinematics->_eps_vec->conjugate_component(mu, lam_vec, _s, _theta);
                temp1 *= METRIC[mu];
                temp1 *= _kinematics->_eps_gamma->component(mu, lam_gam, _s, 0.);
                temp1 *= _kinematics->_initial_state->q(nu, _s, 0.);
                temp1 *= METRIC[nu];
                temp1 *= _kinematics->_final_state->q(nu, _s, _theta);

                term1 += temp1;

                // (eps*_lam . q_gam)(eps_gam . q_vec)
                std::complex<double> temp2;
                temp2  = _kinematics->_eps_vec->conjugate_component(mu, lam_vec, _s, _theta);
                temp2 *= METRIC[mu];
                temp2 *= _kinematics->_initial_state->q(mu, _s, 0.);
                temp2 *= _kinematics->_eps_gamma->component(nu, lam_gam, _s, 0.);
                temp2 *= METRIC[nu];
                temp2 *= _kinematics->_final_state->q(nu, _s, _theta);

                term2 += temp2;
            }
        }

        result = (term1 - term2) / _kinematics->_mX;
    }

    // V - V - P
    if (_kinematics->_jp == VECTOR)
    {
        // Contract with LeviCivita
        for (int mu = 0; mu < 4; mu ++)
        {
            for (int alpha = 0; alpha < 4; alpha++)
            {
                for (int beta = 0; beta < 4; beta++)
                {
                    for (int gamma = 0; gamma < 4; gamma++)
                    {
                        std::complex<double> temp;
                        temp = levi_civita(mu, alpha, beta, gamma);
                        if (std::abs(temp) < 0.001) continue;
                        temp *= _kinematics->_eps_vec->conjugate_component(mu, lam_vec, _s, _theta);
                        temp *= _kinematics->_eps_gamma->field_tensor(alpha, beta, lam_gam, _s, 0.);
                        temp *= _kinematics->_final_state->q(gamma, _s, _theta) - _kinematics->t_exchange_momentum(gamma, _s, _theta);

                        if (!_kinematics->_photon) temp /= 4.;
                        result += temp;
                    }
                }
            }
        }
    }

    // V - P - P
    if (_kinematics->_jp == PSEUDO_SCALAR)
    {
        if (_kinematics->_mB < 1.E-3) return 0.;

        for (int mu = 0; mu < 4; mu++)
        {
            std::complex<double> temp;
            temp   = _kinematics->_eps_gamma->component(mu, lam_gam, _s, 0.);
            temp  *= METRIC[mu];
            temp  *= _kinematics->_final_state->q(mu, _s, _theta) - _kinematics->t_exchange_momentum(mu, _s, _theta);
            result += XI * temp;
        }
    }

    return _gGamma * result;
};

//------------------------------------------------------------------------------
// Simple pole propagator
std::complex<double> jpacPhoto::pseudoscalar_exchange::scalar_propagator()
{
    if (_reggeized == false)
    {
        return 1. / (_t - _mEx2);
    }
    else
    {
        std::complex<double> alpha_t = _alpha->eval(_t);

        if (std::abs(alpha_t) > 20.) return 0.;

        // Else use the regge propagator
        std::complex<double> result = 1.;
        result  = - _alpha->slope();
        result *= 0.5 * (double(_alpha->_signature) +  exp(-XI * PI * alpha_t));
        result *= cgamma(0. - alpha_t);
        result *= pow(_s, alpha_t);
        return result;
    }
};
