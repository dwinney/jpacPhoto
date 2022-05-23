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
    update(helicities, s, t + EPS);

    // Output
    std::complex<double> result;

    if ((_useCovariant == true))
    {
        result = covariant_amplitude();
    }
    else
    {
        result = analytic_amplitude();
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
        case 1: return exp((_t - _kinematics->t_man(_s, 0.)) / _cutoff*_cutoff);
        // monopole form factor
        case 2: return (_cutoff*_cutoff - _mEx2) / (_cutoff*_cutoff - _t); 
        // No form-factor
        default: return 1.;
    }
};

// ---------------------------------------------------------------------------
// REGGE EVALUATION
// ---------------------------------------------------------------------------

std::complex<double> jpacPhoto::vector_exchange::analytic_amplitude()
{
    _lam  =  _lam_gam - _lam_vec;
    _lamp = (_lam_tar - _lam_rec) / 2.;
    _M    = std::max(std::abs(_lam), std::abs(_lamp));

    _zt = std::real(_kinematics->z_t(_s, _theta));
    _qt = sqrt(XR * Kallen(_t, _mX*_mX, _mB*_mB)) / sqrt(4. * _t * XR);

    if (abs(_lam) == 2) return 0.; // double flip forbidden!

    // Product of residues  
    std::complex<double> result;
    result  = top_residue();
    result *= bottom_residue();

    // Pole with d function residue if fixed spin
    if (_ifReggeized == false)
    {
        result *= wigner_d_int_cos(1, _lam, _lamp, _zt);
        result /= _t - _mEx2;
    }
    // or regge propagator if reggeized
    else
    {
        result *= regge_propagator();
    }  

    return result;
};
// ---------------------------------------------------------------------------
// Analytic residues for Regge form

// ---------------------------------------------------------------------------
// Top residue

std::complex<double> jpacPhoto::vector_exchange::top_residue()
{
    if (abs(_lam) >= 2) return 0.;

    int jp = 10*_kinematics->_jp[0] + (1+_kinematics->_jp[1])/2;
    switch (jp)
    {
        case 11: return axialvector_residue();
        case 10: return vector_residue();
        case  1: return pseudoscalar_residue();
        default: return 0.;
    };
};

std::complex<double> jpacPhoto::vector_exchange::axialvector_residue()
{
    std::complex<double> result;
    (abs(_lam) == 0) ? (result = 1.) : (result = sqrt(XR * _t) / _mX);
    return _gGam * _qt * result;
};

std::complex<double> jpacPhoto::vector_exchange::vector_residue()
{
    std::complex<double> result;
    (abs(_lam) == 0) ? (result = -1. + (1- abs(_lam_gam)) *_mB/_mX) : (result = - sqrt(XR * _t) / _mX);
    return _gGam * _qt * result;
};

std::complex<double> jpacPhoto::vector_exchange::pseudoscalar_residue()
{
    std::complex<double> result =  sqrt(XR * _t);
    if (_kinematics->is_photon()) result *=  -4.;
    return _gGam * _qt * result;
};

// ---------------------------------------------------------------------------
// Bottom 

// Nucleon - Nucleon - Vector
std::complex<double> jpacPhoto::vector_exchange::bottom_residue()
{
    std::complex<double> vector, tensor;
    if (abs(_lamp) == 0)
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
    result *= double(_lam_tar);

    return result;
};

// ---------------------------------------------------------------------------
// Reggeon Propagator
std::complex<double> jpacPhoto::vector_exchange::regge_propagator()
{
    if (_M > 1) return 0.;

    std::complex<double> alpha_t = _alpha->eval(_t);

    // the gamma function causes problesm for large t so
    if (std::abs(alpha_t) > 30.)
    {
        return 0.;
    }
    else
    {
        std::complex<double> result;
        result  = wigner_leading_coeff(1, _lam, _lamp);
        result /= barrier_factor();
        result *= half_angle_factor();

        result *= - _alpha->slope();
        result *= 0.5 * (double(_alpha->_signature) + exp(-XI * PI * alpha_t));
        result *= cgamma(1. - alpha_t);
        result *= pow(_s, alpha_t - double(_M));

        return result;
    }
};

//------------------------------------------------------------------------------
// Half angle factors
std::complex<double> jpacPhoto::vector_exchange::half_angle_factor()
{
    std::complex<double> sinhalf = sqrt((XR - _zt) / 2.);
    std::complex<double> coshalf = sqrt((XR + _zt) / 2.);

    std::complex<double> result;
    result  = pow(sinhalf, double(std::abs(_lam - _lamp)));
    result *= pow(coshalf, double(std::abs(_lam + _lamp)));

    return result;
};

//------------------------------------------------------------------------------
// Angular momentum barrier factor
std::complex<double> jpacPhoto::vector_exchange::barrier_factor()
{
    std::complex<double> q = (_t - _mX*_mX) / sqrt(4. * _t * XR);
    std::complex<double> p = sqrt(XR * _t - 4.* M2_PROTON) / 2.;

    std::complex<double> result = pow(2. * p * q, double(1 - _M));

    return result;
};

// ---------------------------------------------------------------------------
// FEYNMAN EVALUATION
// ---------------------------------------------------------------------------

std::complex<double> jpacPhoto::vector_exchange::covariant_amplitude()
{
    std::complex<double> result = 0.;

    // Need to contract the Lorentz indices
    for (int mu = 0; mu < 4; mu++)
    {
        for(int nu = 0; nu < 4; nu++)
        {
            std::complex<double> temp;
            temp  = top_vertex(mu);
            temp *= METRIC[mu];
            temp *= vector_propagator(mu, nu);
            temp *= METRIC[nu];
            temp *= bottom_vertex(nu);

            result += temp;
        }
    }

    return result;
};

// ---------------------------------------------------------------------------
// Photon - Axial Vector - Vector vertex
std::complex<double> jpacPhoto::vector_exchange::top_vertex(int mu)
{
        int jp = 10 * _kinematics->_jp[0] + (1+_kinematics->_jp[1])/2;
    switch (jp)
    {
        case 11: return axialvector_coupling(mu);
        case 10: return vector_coupling(mu);
        case  1: return scalar_coupling(mu);
        case  0: return pseudoscalar_coupling(mu);

        default: return 0.;
    };
}

std::complex<double> jpacPhoto::vector_exchange::axialvector_coupling(int mu)
{
    std::complex<double> result = 0.;

    //Contract with LeviCivita
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
                temp *= _covariants->beam_momentum(alpha);
                temp *= _covariants->beam_polarization(beta);
                temp *= _covariants->meson_polarization(gamma);
                result += temp;
            }
        }
    }

    return _gGam * result;
};

std::complex<double> jpacPhoto::vector_exchange::vector_coupling(int mu)
{
    std::complex<double> result = 0.;

    for (int nu = 0; nu < 4; nu++)
    {
        std::complex<double> temp = -1.;
        temp *= _covariants->beam_field_tensor(mu, nu);
        temp *= METRIC[nu];
        temp *= _covariants->meson_polarization(nu);

        result += temp;
    }
    return _gGam * result;
};

std::complex<double> jpacPhoto::vector_exchange::pseudoscalar_coupling(int mu)
{  
    std::complex<double> result = 0.;

    // Contract with LeviCivita
    for (int alpha = 0; alpha < 4; alpha++)
    {
        for (int beta = 0; beta < 4; beta++)
        {
            for (int gamma = 0; gamma < 4; gamma++)
            {
                std::complex<double> temp = XI;

                temp *= levi_civita(mu, alpha, beta, gamma);
                if (std::abs(temp) < EPS) continue;

                temp *= _covariants->meson_momentum(alpha);
                temp *= _covariants->beam_polarization(beta);
                temp *= _covariants->t_momentum(gamma);

                result += temp;
            }
        }
    }

    if (_kinematics->is_photon()) result *= - 4.;
    return _gGam * result;
};

std::complex<double> jpacPhoto::vector_exchange::scalar_coupling(int mu)
{
    std::complex<double> result = 0.;

    for (int nu = 0; nu < 4; nu++)
    {
        std::complex<double> term1, term2;

        // (k . q) eps_gamma^mu
        term1  = _covariants->t_momentum(nu);
        term1 *= METRIC[nu];
        term1 *= _covariants->beam_momentum(nu);
        term1 *= _covariants->beam_polarization(mu);


        // (eps_gam . k) q^mu
        term2  = _covariants->beam_polarization(nu);
        term2 *= METRIC[nu]; 
        term2 *= _covariants->t_momentum(nu);
        term2 *= _covariants->beam_momentum(mu);

        result += term1 - term2;
    }

    // Dimensionless coupling requires dividing by the mX
    return _gGam * result / _mX;
};

// ---------------------------------------------------------------------------
// Nucleon - Nucleon - Vector vertex
std::complex<double> jpacPhoto::vector_exchange::bottom_vertex(int mu)
{
    // Vector coupling piece
    std::complex<double> vector = 0.;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            std::complex<double> temp;
            temp  = _covariants->recoil_spinor(i);
            temp *= GAMMA[mu][i][j];
            temp *= _covariants->target_spinor(j);
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
                    sigma_q_ij += sigma(mu, nu, i, j) * METRIC[nu] * _covariants->t_momentum(nu) / (_mT + _mR);
                }

                std::complex<double> temp;
                temp  = _covariants->recoil_spinor(i);
                temp *= sigma_q_ij;
                temp *= _covariants->target_spinor(j);
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
    result = _covariants->t_momentum(mu) * _covariants->t_momentum(nu) / _t;

    if (mu == nu)
    {
        result -= METRIC[mu];
    }

    result /= _t - _mEx2;

    return result;
};