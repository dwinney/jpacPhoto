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
    // Save inputs 
    update(helicities, s, t);

    std::complex<double> result;
    if (_useCovariant == true || _debug == 1)
    {
        // If we're using covariants we need to additionally save inputs in covariant_kinematics 
        _covariants->update(helicities, s, t);

        // Because its a scalar exchange we dont have any loose indices to contract
        result  = scalar_propagator();
        result *= top_vertex();
        result *= bottom_vertex();
    }
    else
    {
        _qt = sqrt(XR * Kallen(_t, _mX*_mX, _mB*_mB)) / sqrt(4. * _t * XR);

        if (_lam_vec != _lam_gam || _lam_tar != _lam_rec) 
        {
            return 0.; 
        }
        else
        {
            result  = scalar_propagator();
            result *= top_residue();
            result *= bottom_residue();
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

// To 
std::complex<double> jpacPhoto::pseudoscalar_exchange::top_residue()
{
    // Spin-flip not allowed with spin-0 exchange
    if (_lam_gam != _lam_vec) return 0.;

    std::complex<double> result;

    if (_kinematics->_jp == AXIAL_VECTOR)
    {
        result  = 1. / _mX;
    }
    else if (_kinematics->_jp == VECTOR)
    {
        result = -1.; 
        (_kinematics->is_photon()) ? (result *= -4.) : (result *= 1.);
    }

    else if (_kinematics->_jp == PSEUDO_SCALAR)
    {
        (_kinematics->is_photon()) ? (result =  0.) : (result = 2.*XI/_mB );
    }

    return 2. * _gGamma * _qt * sqrt(XR * _t) * result;
};

// Nucleon resiude 
std::complex<double> jpacPhoto::pseudoscalar_exchange::bottom_residue()
{
    std::complex<double> result;

    if (_lam_tar != _lam_rec) return 0.;

    result  = _gNN; 
    result *= sqrt(XR * _t - pow((_mT - _mR), 2.));
    result /= 2.;

    return result;
};

//------------------------------------------------------------------------------
// FEYNMAN EVALUATION

// Target -- Exchange -- Recoil Vertex
std::complex<double> jpacPhoto::pseudoscalar_exchange::bottom_vertex()
{
    std::complex<double> result = 0.;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            // ubar(recoil) * gamma_5 * u(target)
            std::complex<double> temp;
            temp  = _covariants->recoil_spinor(i);
            temp *= GAMMA_5[i][j];
            temp *= _covariants->target_spinor(j);

            result += temp;
        }
    }

    result *= _gNN;

    return result;
};

// Beam -- Exchange -- Meson vertex
std::complex<double> jpacPhoto::pseudoscalar_exchange::top_vertex()
{
    int jp = 10 * _kinematics->_jp[0] + (1+_kinematics->_jp[1])/2;
    switch (jp)
    {
        case 11: return axialvector_coupling();
        case 10: return vector_coupling();
        case  0: return pseudoscalar_coupling();

        default: return 0.;
    };
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

    return 0.;
};

//------------------------------------------------------------------------------
// These are the different top couplings which depend on the quantum numbers 

std::complex<double> jpacPhoto::pseudoscalar_exchange::axialvector_coupling()
{    
    std::complex<double> term1 = 0., term2 = 0.;
    for (int mu = 0; mu < 4; mu++)
    {
        for (int nu = 0; nu < 4; nu++)
        {
            // (eps*_lam . eps_gam)(q_vec . q_gam)
            std::complex<double> temp1;
            temp1  = _covariants->meson_polarization(mu);
            temp1 *= METRIC[mu];
            temp1 *= _covariants->beam_polarization(mu);
            temp1 *= _covariants->beam_momentum(nu);
            temp1 *= METRIC[nu];
            temp1 *= _covariants->meson_momentum(nu);

            term1 += temp1;

            // (eps*_lam . q_gam)(eps_gam . q_vec)
            std::complex<double> temp2;
            temp2  = _covariants->meson_polarization(mu);
            temp2 *= METRIC[mu];
            temp2 *= _covariants->beam_momentum(mu);
            temp2 *= _covariants->beam_polarization(nu);
            temp2 *= METRIC[nu];
            temp2 *= _covariants->meson_momentum(nu);

            term2 += temp2;
        }
    }

    return _gGamma * (term1 - term2) / _mX;
};

std::complex<double> jpacPhoto::pseudoscalar_exchange::vector_coupling()
{
    std::complex<double> result = 0.;

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
                    temp = - levi_civita(mu, alpha, beta, gamma);
                    if (std::abs(temp) < 0.001) continue;
                    temp *= _covariants->meson_polarization(mu);
                    temp *= _covariants->beam_field_tensor(alpha, beta);
                    temp *= _covariants->meson_momentum(gamma);

                    result += temp;
                }
            }
        }
    }
    
    if (!_kinematics->is_photon()) result /= -4.;
    return _gGamma * result;
};

std::complex<double> jpacPhoto::pseudoscalar_exchange::pseudoscalar_coupling()
{
    if (_kinematics->is_photon()) return 0.;

    std::complex<double> result = 0.;
    for (int mu = 0; mu < 4; mu++)
    {
        std::complex<double> temp;
        temp  = _covariants->beam_polarization(mu);
        temp *= METRIC[mu];
        temp *= _covariants->beam_momentum(mu);
        result += -XI * temp;
    }
    return _gGamma * result;
};