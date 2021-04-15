// Vector production via a one-loop box diagram
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "one_loop/box_amplitude.hpp"

// ---------------------------------------------------------------------------
// Evaluate the entire loop diagram in terms of external helicites and invariant mass / momentum transfer of the gamma p -> jpsi p process
std::complex<double> jpacPhoto::box_amplitude::helicity_amplitude(std::array<int,4> helicities, double s, double t)
{
    // Store the invariant energies to avoid having to pass them around 
    _s = s; _t = t, _theta = _kinematics->theta_s(s, t);

    std::complex<double> result = 0.;
    // Pass external values to the discontinuity
    for (int n = 0; n < _discs.size(); n++)
    {
        _discs[n]->set_externals(helicities, _theta);

        double sub =  _discs[n]->eval(s);
        auto F = [&](double sp)
        {
            std::complex<double> result = (_discs[n]->eval(sp) - sub) / (sp - s - IEPS);
            return result;
        };

        std::complex<double> intpiece = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(F, _discs[n]->_threshold + EPS, _s_cut, 0, 1.E-6, NULL);
        std::complex<double> logpiece = sub * (log(_s_cut - s - IEPS) - log(_discs[n]->_threshold + EPS - s - IEPS));
        std::complex<double> result_n =  (intpiece + logpiece) / M_PI;
        result += result_n;
    }

    return result;
};

// ---------------------------------------------------------------------------
// Override the usual integrated_xsection to use a gauss-legendre integrator since t behavior is smooth but extremely slow
double jpacPhoto::box_amplitude::integrated_xsection(double s)
{
    auto F = [&](double t)
    {
        double result = differential_xsection(s, t);
        debug(t, result);
        return result;
    };

    ROOT::Math::GaussLegendreIntegrator ig(10);
    ROOT::Math::Functor1D wF(F);
    ig.SetFunction(wF);

    double t_min = _kinematics->t_man(s, 0.);
    double t_max = _kinematics->t_man(s, PI);

    return ig.Integral(t_max, t_min);
};