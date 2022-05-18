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
    int lam_gam = helicities[0];
    int lam_targ = helicities[1];
    int lam_vec = helicities[2];
    int lam_rec = helicities[3];

    // Save energies 
    _s = s; _t = t; _theta = _kinematics->theta_s(s, t);
 
    std::complex<double> result = 0.;

    // IF using helicity conserving delta fuction model
    if (_model == 1)
    {
        (lam_gam == lam_vec && lam_rec == lam_targ) ? (result = regge_factor()) : (result = 0.);
        return result;
    }

    // else contract indices
    for (int mu = 0; mu < 4; mu++)
    {
        std::complex<double> temp = 1.;
        temp *= regge_factor();
        temp *= top_vertex(mu, lam_gam, lam_vec);
        temp *= METRIC[mu];
        temp *= bottom_vertex(mu, lam_targ, lam_rec);

        result += temp;
    }

    return result;
};

// ---------------------------------------------------------------------------
// Bottom vertex coupling the target and recoil proton spinors to the vector pomeron
std::complex<double> jpacPhoto::pomeron_exchange::bottom_vertex(int mu, int lam_targ, int lam_rec)
{
    std::complex<double> result = 0.;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            std::complex<double> temp;
            // Recoil oriented an angle theta + pi
            temp = _kinematics->_recoil->adjoint_component(i, lam_rec, _s, _theta);

            // vector coupling
            temp *= GAMMA[mu][i][j];

            // target oriented in negative z direction
            temp *= _kinematics->_target->component(j, lam_targ, _s, 0.);

            result += temp;
        }
    }

    return result;
};

// ---------------------------------------------------------------------------
// Top vertex coupling the photon, pomeron, and vector meson.
std::complex<double> jpacPhoto::pomeron_exchange::top_vertex(int mu, int lam_gam, int lam_vec)
{
    std::complex<double> result = 0.;

    if (_model == 0)
    {
        std::complex<double> sum1 = 0., sum2 = 0.;
        for (int nu = 0; nu < 4; nu++)
        {
            std::complex<double> temp1, temp2;

            // (q . eps_vec^*) eps_gam^mu
            temp1  = _kinematics->_initial_state->q(nu, _s, 0.);
            temp1 *= METRIC[nu];
            temp1 *= _kinematics->_eps_vec->conjugate_component(nu, lam_vec, _s, _theta);
            sum1  += _kinematics->_eps_gamma->component(mu, lam_gam, _s, 0.) * temp1;

            // (eps_vec^* . eps_gam) q^mu
            temp2  = _kinematics->_eps_gamma->component(nu, lam_gam, _s, 0.);
            temp2 *= METRIC[nu];
            temp2 *= _kinematics->_eps_vec->conjugate_component(nu, lam_vec, _s, _theta);
            sum2  += _kinematics->_initial_state->q(mu, _s, 0.) * temp2;
        }

        result = -sum1 + sum2;
    }
    else if (_model == 2)
    {
        std::complex<double> sum1 = 0., sum2 = 0.;
        for (int nu = 0; nu < 4; nu++)
        {
            std::complex<double> temp1, temp2;

            // -2 * (q . eps_vec^*) eps_gam^mu
            temp1  = _kinematics->_initial_state->q(nu, _s, 0.);
            temp1 *= METRIC[nu];
            temp1 *= _kinematics->_eps_vec->conjugate_component(nu, lam_vec, _s, _theta);
            sum1  += -2. * _kinematics->_eps_gamma->component(mu, lam_gam, _s, 0.) * temp1;

            // (eps_vec . eps_gam) (q + q')^mu
            temp2  = _kinematics->_eps_vec->conjugate_component(nu, lam_vec, _s, _theta);
            temp2 *= METRIC[nu];
            temp2 *= _kinematics->_eps_gamma->component(nu, lam_gam, _s, 0.);
            sum2  += (_kinematics->_initial_state->q(mu, _s, 0.) + _kinematics->_final_state->q(mu, _s, _theta)) * temp2;
        }
      
        result = (sum1 + sum2);
    };

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
            double mX2 = _kinematics->_mX2;
            double th  = pow((_kinematics->_mT + _kinematics->_mR), 2.);

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
