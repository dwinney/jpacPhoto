// Charged axial-vector meson photoproduction proceeding through a pseudoscalar (pion) exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// References:
// [1] arXiv:1503.02125 [hep-ph]
// ---------------------------------------------------------------------------

#include "pseudoscalar_exchange.hpp"

//------------------------------------------------------------------------------
// Combine everything and contract indices
std::complex<double> jpacPhoto::pseudoscalar_exchange::helicity_amplitude(std::array<int, 4> helicities, double s, double t)
{
    // Save inputs 
    update(helicities, s, t);

    std::complex<double> result;
    if (_useCovariant && !_reggeized)
    {
        // Because its a scalar exchange we dont have any loose indices to contract
        result  = scalar_propagator();
        result *= top_vertex();
        result *= bottom_vertex();
    }
    else
    {
        _qt = sqrt(XR * Kallen(_t, _mX*_mX, _mB*_mB)) / sqrt(4. * _t * XR);
        _pt = sqrt(XR * Kallen(_t, _mR*_mR, _mT*_mT)) / sqrt(4. * _t * XR);

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
            return exp((_t - _kinematics->t_man(_s, 0.)) / (_cutoff*_cutoff));
        };
        // monopole form factor
        case 2:
        {
            return (_cutoff*_cutoff - _mEx2) / (_cutoff*_cutoff - _t); 
        };
        case 3:
        {
            double u = _kinematics->u_man(_s, _theta);

            auto Fx = [&] (double x, double m)
            {
                return pow(_cutoff, 4.) / (pow(_cutoff, 4.) + pow(x - m*m, 2.));
            };

            return 1. - (1. - Fx(_s, M_PROTON)) * (1. - Fx(_t, M_PION)) * (1. - Fx(u, M_LAMBDA));
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
    
    std::array<int,2> JP = _kinematics->get_meson_JP();
    int jp = 10 * JP[0] + (1+JP[1])/2;

    switch (jp)
    {
        case (11): result = 1./_mX; break;
        case (10): result = _kinematics->is_photon() * -4. + !_kinematics->is_photon() * -1.; break;
        case ( 0): result = _kinematics->is_photon() * 2.*XI/_mB; break;
        default: return 0.;
    };

    return 2. * _gT * _qt * sqrt(XR * _t) * result;
};

// Nucleon resiude 
std::complex<double> jpacPhoto::pseudoscalar_exchange::bottom_residue()
{
    if (_lam_tar != _lam_rec) return 0.
    ;

    std::array<int,2> JP = _kinematics->get_baryon_JP();
    int jp = 10 * JP[0] + (1+JP[1])/2;

    std::complex<double> result;
    switch (jp)
    {
        case (11): result = sqrt(XR * _t - pow((_mT - _mR), 2.)) / 2.; break;
        case (31): 
        {
            result  = sqrt(2./3.) * _pt * sqrt(XR * _t) / _mR;
            result *= sqrt(XR * _t - pow((_mT + _mR), 2.)) / 2.;
            break;
        };
        default: return 0.;
    }

    return _gB * result;
};

//------------------------------------------------------------------------------
// FEYNMAN EVALUATION

std::complex<double> jpacPhoto::pseudoscalar_exchange::bottom_vertex()
{
    std::array<int,2> JP = _kinematics->get_baryon_JP();
    int jp = 10 * JP[0] + (1+JP[1])/2;

    switch (jp)
    {
        case 11: return halfplus_coupling();
        case 31: return threehalvesplus_coupling();
        default: return 0.;
    };
};


// Beam -- Exchange -- Meson vertex
std::complex<double> jpacPhoto::pseudoscalar_exchange::top_vertex()
{
    std::array<int,2> JP = _kinematics->get_meson_JP();
    int jp = 10 * JP[0] + (1+JP[1])/2;
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
// These are the different bottom couplings which depend on the quantum numbers 

// Target -- Exchange -- Recoil Vertex
std::complex<double> jpacPhoto::pseudoscalar_exchange::halfplus_coupling()
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

    result *= _gB;

    return result;
};

// Target -- Exchange -- Recoil Vertex
std::complex<double> jpacPhoto::pseudoscalar_exchange::threehalvesplus_coupling()
{
    std::complex<double> result = 0.;
    for (int i = 0; i < 4; i++)
    {
        for (int mu = 0; mu < 4; mu ++)
        {
            std::complex<double> temp;
            temp  = _covariants->recoil_spinor(i, mu);
            temp *= METRIC[mu];
            temp *= _covariants->t_momentum(mu);
            temp *= _covariants->target_spinor(i);
            result += temp;
        }
    }
    result *= _gB;

    return result;
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

    return _gT * (term1 - term2) / _mX;
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
                    temp *= _covariants->meson_momentum(gamma) - _covariants->t_momentum(gamma);

                    result += temp;
                }
            }
        }
    }
    
    if (!_kinematics->is_photon()) result /= -4.;
    return _gT * result;
};

std::complex<double> jpacPhoto::pseudoscalar_exchange::pseudoscalar_coupling()
{
    std::complex<double> result = 0.;
    for (int mu = 0; mu < 4; mu++)
    {
        std::complex<double> temp;
        temp  = _covariants->beam_polarization(mu);
        temp *= METRIC[mu];
        temp *= _covariants->t_momentum(mu) - _covariants->meson_momentum(mu);
        result += - temp; 
    }

    return _gT * result;
};