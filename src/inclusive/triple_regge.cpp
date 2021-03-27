// Form of the invariant cross-section from a triple regge interaction.
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "inclusive/triple_regge.hpp"

double jpacPhoto::triple_regge::invariant_xsection(double s, double t, double M2)
{
    std:double result = 0.;

    for (int i = 0; i < _termsFF.size(); i++)
    {
        result += _termsFF[i]->eval(s, t, M2);
    }
    for (int i = 0; i < _termsJPAC.size(); i++)
    {
        result += _termsJPAC[i]->eval(s, t, M2);
    }

    return result; // in nanobarn!
};

double jpacPhoto::triple_regge::integrated_xsection(double s)
{
    auto dSigma = [&](const double * in)
    {
        double M2 = _kinematics->M2(s, in[1]);
        double t  = _kinematics->t_man(s, in[0], M2);
        
        double result = invariant_xsection(s, t, M2);
        debug(t, in[0], in[1]);
        return result;
    };

    // Integrate over costheta = [-1, 1] and x = [0., 1]   
    double min[2] = {-1., 0.1};
    double max[2] = { 1., 1.};
    
    ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kDEFAULT);
    ROOT::Math::Functor wF(dSigma, 2);
    ig.SetFunction(wF, 2);

    double result = ig.Integral(min, max);

    return result;
};