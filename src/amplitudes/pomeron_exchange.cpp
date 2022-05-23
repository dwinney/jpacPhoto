// Vector meson photoproduction dynamics proceeding through a pomeron exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/pomeron_exchange.hpp"

// ---------------------------------------------------------------------------
// Given a set of helicities for each particle, assemble the helicity amplitude by contracting Lorentz indicies
std::complex<double> jpacPhoto::pomeron_exchange::helicity_amplitude(std::array<int, 4> helicities, double s, double t)
{
    // Save energies and helicities
    update(helicities, s, t);

    // We use covariants to eval the amplitude so update that too
    _covariants->update(helicities, s, t);
 
    std::complex<double> result = 0.;
    // IF using helicity conserving delta fuction model
    if (_model == 1)
    {
        (_lam_gam == _lam_vec && _lam_rec == _lam_tar) ? (result = regge_factor()) : (result = 0.);
        return result;
    }

    // else contract indices
    for (int mu = 0; mu < 4; mu++)
    {
        std::complex<double> temp = 1.;
        temp *= regge_factor();
        temp *= top_vertex(mu);
        temp *= METRIC[mu];
        temp *= bottom_vertex(mu);

        result += temp;
    }

    return result;
};

// ---------------------------------------------------------------------------
// Bottom vertex coupling the target and recoil proton spinors to the vector pomeron
std::complex<double> jpacPhoto::pomeron_exchange::bottom_vertex(int mu)
{
    std::complex<double> result = 0.;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            std::complex<double> temp;
            temp  = _covariants->recoil_spinor(i);
            temp *= GAMMA[mu][i][j];
            temp *= _covariants->target_spinor(j);

            result += temp;
        }
    }

    return result;
};

// ---------------------------------------------------------------------------
// Top vertex coupling the photon, pomeron, and vector meson.
std::complex<double> jpacPhoto::pomeron_exchange::top_vertex(int mu)
{
    std::complex<double> result = 0.;

    std::complex<double> sum1 = 0., sum2 = 0.;
    for (int nu = 0; nu < 4; nu++)
    {
        std::complex<double> temp1, temp2;

        // (q . eps_vec^*) eps_gam^mu
        temp1  = _covariants->beam_momentum(nu);
        temp1 *= METRIC[nu];
        temp1 *= _covariants->meson_polarization(nu);
        sum1  += temp1 * _covariants->beam_polarization(mu);

        // (eps_vec^* . eps_gam) q^mu
        temp2  = _covariants->beam_polarization(nu);
        temp2 *= METRIC[nu];
        temp2 *= _covariants->meson_polarization(nu);
        sum2  += temp2 * _covariants->beam_momentum(mu);
    }

    result = -sum1 + sum2;

    return result;
};

// ---------------------------------------------------------------------------
// Usual Regge power law behavior, s^alpha(t) with an exponential fall from the forward direction
std::complex<double> jpacPhoto::pomeron_exchange::regge_factor()
{
    if (_s < _kinematics->sth())
    {
        std::cout << " \n pomeron_exchange: Trying to evaluate below threshold (sqrt(s) = " << sqrt(_s) << ")! Quitting... \n";
        exit(0);
    }

    std::complex<double> result = 0.;
    
    switch (_model)
    {
        case 0:
        {
            double t_min = _kinematics->t_man(_s, 0.); // t_min = t(theta = 0)
            result  = exp(_b0 * (_t - t_min));
            result *= pow(_s - _kinematics->sth(), _traj->eval(_t));
            result *= XI * _norm * E;
            result /= _s;
            break;
        }
        case 1:
        {
            double t_min = _kinematics->t_man(_s, 0.); // t_min = t(theta = 0)
            result  = exp(_b0 * (_t - t_min));
            result *= pow(_s - _kinematics->sth(), _traj->eval(_t));
            result *= XI * _norm * E;
            break;
        }
        case 2:
        {
            double mX2 = _mX*_mX;
            double th  = pow((_mT + _mR), 2.);

            double beta_0 = 2.;           // Pomeron - light quark coupling
            double beta_c = _norm;        // Pomeron - charm quark coupling
            double mu2 = _b0 * _b0;       // cutoff parameter 
            double etaprime = real(_traj->slope());

            std::complex<double> F_t;
            F_t  = 3. * beta_0;
            F_t *= (th - 2.8* _t);
            F_t /= (th - _t) *  pow((1. - (_t / 0.7)) , 2.);

            std::complex<double> G_p = -XI;
            G_p  *= pow(XR * etaprime * _s, _traj->eval(_t) - 1.);

            result  = - XI * 8. * beta_c * mu2 * G_p * F_t;
            result *= 2. * E * F_JPSI / M_JPSI; // Explicitly only for the jpsi... 
            result /= (mX2 - _t) * (2.*mu2 + mX2 - _t);
            break;
        }
        default: return 0.;
    }

    return result;
};
