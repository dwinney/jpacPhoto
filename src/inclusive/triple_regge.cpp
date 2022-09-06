// Form of the triple Regge interaction using JPAC's parameterization
// i.e. using the t dependence from properly normalized Regge propagators
// and M2 dependence from the total hadronic cross-section of the bottom vertex.
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "triple_regge.hpp"

// ---------------------------------------------------------------------------
// Parse the passed amplitude_name string from an exclusive amplitude
// to make sure the appropriate coupling function is used for the top vertex
void jpacPhoto::triple_regge::initialize(std::string amp_name)
{
    // Pion exchange on top vertices
    if (amp_name == "pseudoscalar_exchange")
    {
        // Axial-vector - photon - psuedoscalar coupling
        _coupling = [&](double t)
        {
            double ff =  exp(_b * (t - _kinematics->TMINfromM2( M2_PROTON )));
            return   ff * (_g / _kinematics->_mX) * (t - _kinematics->_mX2) / 2.;
        };

        // Default: pi- exchange with the PDG parameterization (no resonances)
        _sigma_tot = new PDG_parameterization(M_PION, M_PROTON, {-1., 1., 9.56, 1.767, 18.75});
        _exchange_spin = 0;

        if (_useRegge)
        {
            _alpha_prime = std::real(_trajectory->slope());
            _cutoff = exp(_b / _alpha_prime) / _alpha_prime;
        } 
    }
    else 
    {
        _coupling = [&](double x)
        {
            return 0.;
        };

        _sigma_tot = new zero_xsection();
    };
};

// ---------------------------------------------------------------------------
// Change the default sigma_total from initialize to a user-selected one
void jpacPhoto::triple_regge::set_sigma_total(sigma_option opt)
{
    // We need to make sure to free up _sigma_tot first
    delete _sigma_tot;
    _sigma_tot = get_total_xsection(opt);
};

// ---------------------------------------------------------------------------
// Evaluate the invariant amplitude
double jpacPhoto::triple_regge::d3sigma_d3p(double s, double t, double mm)
{
    // Make sure to pass the CM energy to the kinematics
    _kinematics->_s = s;

    // Things tend to blow up at exactly x = 1
    if (_useTX && (abs(mm - 1) < 0.001)) return 0.;

    double s_piece, phase_space, M2;
    if (_useTX)
    {
        // We are using x and t as our variables
        s_piece = (1. - mm);
        M2 = _kinematics->M2fromTX(t, mm);
        phase_space = s_piece;
    }
    else
    {
        // We are using M2 and t as our variables
        s_piece = mm / s;
        M2 = mm;
        phase_space = sqrt(Kallen(M2, t, M2_PROTON)  / Kallen(s, 0., M2_PROTON) );
    }
    
    // Coupling squared
    double coupling2   = _coupling(t) * _coupling(t);

    // The form of the exchange propagators depends on if we want a regge form or not
    double exchange_propagator2;
    if (_useRegge)
    {
        double alpha      = std::real(_trajectory->eval(t));

        // First check t isnt too big to make the gamma function blow up
        if (std::abs(t) - _cutoff > 0.) return 0.;

        std::complex<double> signature_factor = (1. + double(_trajectory->_signature) * exp(- XI * M_PI * alpha)) / 2.; 
        double t_piece = std::norm(_alpha_prime * signature_factor * cgamma(double(_trajectory->_minJ) - alpha));

        exchange_propagator2 = t_piece * pow(s_piece, -2. * alpha);
    }
    else
    {
        double pole          = 1. / (_exchange_mass2 - t);                         // Simple pole 
        exchange_propagator2 = pole * pole * pow(s_piece, -2. * _exchange_spin);   // Squared
    };

    // Get missing mass to put into the sigma_tot
    double sigma_tot =  _sigma_tot->eval(M2, t) * 1.E6;

    if (_debug == 1)
    {
        debug("s", s);
        debug("t", t);
        debug("M2", M2);
        debug("phase_space", phase_space);
        debug("coup2", coupling2);
        debug("prop2", exchange_propagator2);
        debug("sigma_tot", sigma_tot * 1.E-6);
    };

    return sigma_tot * coupling2 * exchange_propagator2 * phase_space / (16. * pow(M_PI, 3.));
};  