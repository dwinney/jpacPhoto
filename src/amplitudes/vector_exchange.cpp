// Axial-vector meson photoproduction proceeding through a vector meson exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/vector_exchange.hpp"

// ---------------------------------------------------------------------------
// Assemble the helicity amplitude by contracting the lorentz indices
std::complex<double> jpacPhoto::vector_exchange::helicity_amplitude(std::array<int, 4> helicities, double s, double t)
{
    // Save energies and helicities
    update(helicities, s, t);

    _zt = real(_kinematics->z_t(_s, _theta));

    // Output
    std::complex<double> result;

    if ((_useCovariant == true) || (_debug >= 1))
    {
        result = covariant_amplitude(helicities);
    }
    else
    {
        int lam  = _lam_gam - _lam_vec;
        int lamp = (_lam_tar - _lam_rec) / 2.;

        if (abs(lam) == 2) return 0.; // double flip forbidden!

        // Product of residues  
        result  = top_residue(_lam_gam, _lam_vec);
        result *= bottom_residue(_lam_tar, _lam_rec);

        // Pole with d function residue if fixed spin
        if (_ifReggeized == false)
        {
            result *= wigner_d_int_cos(1, lam, lamp, _zt);
            result /= t - _mEx2;
        }
        // or regge propagator if reggeized
        else
        {
            result *= regge_propagator(1, lam, lamp);
        }
    }

    // add form factor if wanted
    result *= form_factor();   
    return result;
;
};

double jpacPhoto::vector_exchange::form_factor()
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
};

// ---------------------------------------------------------------------------
// REGGE EVALUATION
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Analytic residues for Regge form

// Photon coupling
std::complex<double> jpacPhoto::vector_exchange::top_residue(int lam_gam, int lam_vec)
{
    std::complex<double> result = 0.;
    std::complex<double> q_t = sqrt(XR * Kallen(_t, _mX*_mX, _mB*_mB)) / sqrt(4. * _t * XR);

    if (_kinematics->_jp == AXIAL_VECTOR)
    {
        switch (std::abs(lam_gam - lam_vec))
        {
            case 0: result = 1.; break;
            case 1: result = sqrt(XR * _t) / _mX; break;
            default: return 0.;
        };
    }

    else if (_kinematics->_jp == VECTOR)
    {
        switch (std::abs(lam_gam - lam_vec))
        {
            case 0:
            {
                result  = -1. + (1-std::abs(lam_gam)) * _mB / _mX;
                break;
            } 
            case 1: 
            { 
                result = - sqrt(XR * _t) / _mX;
                break;
            }
            default: return 0.;
        };
    }

    else if (_kinematics->_jp == PSEUDO_SCALAR)
    {
        result = sqrt(XR * _t);
        if (_kinematics->is_photon()) result *=  -4.;
    }
    
    return _gGam * q_t * result;
};

// Nucleon - Nucleon - Vector
std::complex<double> jpacPhoto::vector_exchange::bottom_residue(int lam_tar, int lam_rec)
{
    std::complex<double> vector, tensor;
    if (lam_tar == lam_rec)
    {
        vector = (_mT + _mR); 
        tensor = _t / (_mT + _mR);
    }
    else
    {
        vector = sqrt(2.) * sqrt(XR * _t);
        tensor = sqrt(2.) * sqrt(XR * _t);
    }

    std::complex<double> result;
    result  = _gV * vector + _gT * tensor;
    result *= sqrt(XR * _t - pow((_mT - _mR), 2.)) / sqrt(XR * _t);
    result *= double(lam_tar);

    return result;
};

// ---------------------------------------------------------------------------
// Reggeon Propagator
std::complex<double> jpacPhoto::vector_exchange::regge_propagator(int j, int lam, int lamp)
{
    int M = std::max(std::abs(lam), std::abs(lamp));

    if (M > j)
    {
        return 0.;
    }

    std::complex<double> alpha_t = _alpha->eval(_t);

    // the gamma function causes problesm for large t so
    if (std::abs(alpha_t) > 30.)
    {
        return 0.;
    }
    else
    {
        std::complex<double> result;
        result  = wigner_leading_coeff(j, lam, lamp);
        result /= barrier_factor(j, M);
        result *= half_angle_factor(lam, lamp);

        result *= - _alpha->slope();
        result *= 0.5 * (double(_alpha->_signature) + exp(-XI * PI * alpha_t));
        result *= cgamma(1. - alpha_t);
        result *= pow(_s, alpha_t - double(M));

        return result;
    }
};

//------------------------------------------------------------------------------
// Half angle factors
std::complex<double> jpacPhoto::vector_exchange::half_angle_factor(int lam, int lamp)
{
    std::complex<double> sinhalf = sqrt((XR - _zt) / 2.);
    std::complex<double> coshalf = sqrt((XR + _zt) / 2.);

    std::complex<double> result;
    result  = pow(sinhalf, double(std::abs(lam - lamp)));
    result *= pow(coshalf, double(std::abs(lam + lamp)));

    return result;
};

//------------------------------------------------------------------------------
// Angular momentum barrier factor
std::complex<double> jpacPhoto::vector_exchange::barrier_factor(int j, int M)
{
    std::complex<double> q = (_t - _mX*_mX) / sqrt(4. * _t * XR);
    std::complex<double> p = sqrt(XR * _t - 4.* M2_PROTON) / 2.;

    std::complex<double> result = pow(2. * p * q, double(j - M));

    return result;
};

// ---------------------------------------------------------------------------
// FEYNMAN EVALUATION
// ---------------------------------------------------------------------------

std::complex<double> jpacPhoto::vector_exchange::covariant_amplitude(std::array<int, 4> helicities)
{
    int lam_gam = helicities[0];
    int lam_tar = helicities[1];
    int lam_vec = helicities[2];
    int lam_rec = helicities[3];

    std::complex<double> result = 0.;

    // Need to contract the Lorentz indices
    for (int mu = 0; mu < 4; mu++)
    {
        for(int nu = 0; nu < 4; nu++)
        {
            std::complex<double> temp;
            temp  = top_vertex(mu, lam_gam, lam_vec);
            temp *= METRIC[mu];
            temp *= vector_propagator(mu, nu);
            temp *= METRIC[nu];
            temp *= bottom_vertex(nu, lam_tar, lam_rec);

            result += temp;
        }
    }

    return result;
};

// ---------------------------------------------------------------------------
// Photon - Axial Vector - Vector vertex
std::complex<double> jpacPhoto::vector_exchange::top_vertex(int mu, int lam_gam, int lam_vec)
{
    std::complex<double> result = 0.;

    // A-V-V coupling
    if (_kinematics->_jp == AXIAL_VECTOR)
    {
        // Contract with LeviCivita
        for (int alpha = 0; alpha < 4; alpha++)
        {
            for (int beta = 0; beta < 4; beta++)
            {
                for (int gamma = 0; gamma < 4; gamma++)
                {
                    std::complex<double> temp;
                    temp = levi_civita(mu, alpha, beta, gamma);
                    if (std::abs(temp) < 0.001) continue;
                    temp *= METRIC[mu];
                    temp *= _kinematics->_initial_state->q(alpha, _s, 0.);
                    temp *= _kinematics->_eps_gamma->component(beta, lam_gam, _s, 0.);
                    temp *= _kinematics->_eps_vec->component(gamma, lam_vec, _s, _theta);

                    result += temp;
                }
            }
        }
    }

    // V-V-V coupling
    else if (_kinematics->_jp == VECTOR)
    {
        for (int nu = 0; nu < 4; nu++)
        {
            std::complex<double> temp = -1.;
            if (_debug == 2) temp *= METRIC[mu];
            temp *= _kinematics->_eps_gamma->field_tensor(mu, nu, lam_gam, _s, 0.);
            temp *= METRIC[nu];
            temp *= _kinematics->_eps_vec->conjugate_component(nu, lam_vec, _s, _theta);
            result += temp;
        }
    }

    // S-V-V coupling
    else if (_kinematics->_jp == SCALAR)
    {
        for (int nu = 0; nu < 4; nu++)
        {
            std::complex<double> term1, term2;

            // (k . q) eps_gamma^mu
            term1  = _kinematics->t_exchange_momentum(nu, _s, _theta);
            term1 *= METRIC[nu];
            term1 *= _kinematics->_initial_state->q(nu, _s, 0.);
            term1 *= _kinematics->_eps_gamma->component(mu, lam_gam, _s, 0.);

            // (eps_gam . k) q^mu
            term2  = _kinematics->_eps_gamma->component(nu, lam_gam, _s, 0.);
            term2 *= METRIC[nu];
            term2 *= _kinematics->t_exchange_momentum(nu, _s, _theta);
            term2 *= _kinematics->_initial_state->q(mu, _s, 0.);

            result += term1 - term2;
        }

        // Dimensionless coupling requires dividing by the mX
        result /= _mX;
    }

    // P-V-V coupling
    else if (_kinematics->_jp == PSEUDO_SCALAR)
    {
        // Contract with LeviCivita
        for (int alpha = 0; alpha < 4; alpha++)
        {
            for (int beta = 0; beta < 4; beta++)
            {
                for (int gamma = 0; gamma < 4; gamma++)
                {
                    std::complex<double> temp = XI;
                    if (_debug != 2) temp *= METRIC[mu];
                    temp *= levi_civita(mu, alpha, beta, gamma);
                    if (std::abs(temp) < 0.001) continue;

                    // // Explicitly in terms of the Field tensor of the photon
                    // temp *= _kinematics->_eps_gamma->field_tensor(alpha, beta, lam_gam, _s, 0.);
                    // temp *= _kinematics->_final_state->q(gamma, _s, _theta) - _kinematics->t_exchange_momentum(gamma, _s, _theta);

                    // OR reduced to be consistent with jpsi as well
                    temp *= _kinematics->_final_state->q(alpha, _s, _theta);
                    temp *= _kinematics->_eps_gamma->component(beta, lam_gam, _s, 0.);
                    temp *= _kinematics->t_exchange_momentum(gamma, _s, _theta);
                    if (_kinematics->is_photon()) temp *= - 4.;

                    result += temp;
                }
            }
        }
    }

    return result * _gGam;
};

// ---------------------------------------------------------------------------
// Nucleon - Nucleon - Vector vertex
std::complex<double> jpacPhoto::vector_exchange::bottom_vertex(int mu, int lam_tar, int lam_rec)
{
    // Vector coupling piece
    std::complex<double> vector = 0.;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            std::complex<double> temp;
            temp  = _kinematics->_recoil->adjoint_component(i, lam_rec, _s, _theta); // theta_rec = theta + pi
            temp *= GAMMA[mu][i][j];
            temp *= _kinematics->_target->component(j, lam_tar, _s, 0.); // theta_targ = pi

            vector += temp;
        }
    }

    // Tensor coupling piece
    std::complex<double> tensor = 0.;
    if (abs(_gT) > 0.001)
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                std::complex<double> sigma_q_ij = 0.;
                for (int nu = 0; nu < 4; nu++)
                {
                    sigma_q_ij += sigma(mu, nu, i, j) * METRIC[nu] * _kinematics->t_exchange_momentum(nu, _s, _theta) / (2. * M_PROTON);
                }

                std::complex<double> temp;
                temp = _kinematics->_recoil->adjoint_component(i, lam_rec, _s, _theta); // theta_rec = theta + pi
                temp *= sigma_q_ij;
                temp *= _kinematics->_target->component(j, lam_tar, _s, 0.); // theta_targ = pi

                tensor += temp;
            }
        }
    }

    return _gV * vector - _gT * tensor;
};

// ---------------------------------------------------------------------------
// Propagator of a massive spin-one particle
std::complex<double> jpacPhoto::vector_exchange::vector_propagator(int mu, int nu)
{
    // q_mu q_nu / t - g_mu nu
    std::complex<double> result;
    result = _kinematics->t_exchange_momentum(mu, _s, _theta) * _kinematics->t_exchange_momentum(nu, _s, _theta) / _t;

    if (mu == nu)
    {
        result -= METRIC[mu];
    }

    result /= _t - _mEx2;

    return result;
};