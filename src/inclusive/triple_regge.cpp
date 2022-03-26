// Form of the triple Regge interaction using JPAC's parameterization
// i.e. using the t dependence from properly normalized Regge propagators
// and M2 dependence from the total hadronic cross-section of the bottom vertex.
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "inclusive/triple_regge.hpp"

// Evaluate the invariant amplitude
double jpacPhoto::triple_regge::d3sigma_d3p(double s, double t, double mm)
{
    if (_sigma_tot == NULL)
    {
        std::cout << "ERROR! No sigma_tot set! Returning 0! \n";
        return 0.;
    };
    if (_useRegge == true && _trajectory == NULL)
    {
        std::cout << "ERROR! No regge_trajectory set! Returning 0! \n";
        return 0.;
    };
    if (_couplingSet == false)
    {
        std::cout << "ERROR! No coupling set! Returning 0! \n";
        return 0.;   
    };

    // Coupling squared
    double coupling2   = _coupling(t) * _coupling(t);

    // Form factor with tprime corresponding to the exclusive limit
    double formfactor2 = exp(2. * _b * (t - _kinematics->TMINfromM2( _kinematics->_minM2 )));

    // phase_space factor (depends on whether mm is x or M2)
    double s_piece;
    (_useTX) ? (s_piece = (1. - mm)) : (s_piece = mm / s);

    double exchange_propagator2;
    if (_useRegge)
    {
        double alpha = std::real(_trajectory->eval(t));
        double alphaPrime = real(_trajectory->slope());

        // First check t isnt too big to make the gamma function blow up
        if ( _b + alphaPrime - alphaPrime * log(- alphaPrime * t) < 0.) return 0.;

        std::complex<double> signature_factor = (1. + double(_trajectory->_signature) * exp(- XI * M_PI * alpha)) / 2.; 
        double t_piece = std::norm(alphaPrime * signature_factor * gamma(- alpha));

        exchange_propagator2 = t_piece * pow(s_piece, -2. * alpha);
    }
    else
    {
        double pole          = 1. / (_exchange_mass2 - t);  // Simple pole 
        exchange_propagator2 = pole * pole;                // Squared
    };

    double sigma_tot;
    (_useTX) ? (sigma_tot = _sigma_tot->eval(s * (1. - mm))) : (sigma_tot = _sigma_tot->eval(mm));

    return sigma_tot * coupling2 * formfactor2 * exchange_propagator2 * s_piece / pow(4.* M_PI, 3.);
};  