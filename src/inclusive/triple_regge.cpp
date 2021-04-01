// Form of the invariant cross-section from a triple regge interaction.
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "inclusive/triple_regge.hpp"

// ---------------------------------------------------------------------------
// The fully differential invariant cross-section 
// E dsigma/dtdM2
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

// ---------------------------------------------------------------------------
// Singly-differential cross-sections in terms of (t, M2)
double jpacPhoto::triple_regge::dsigma_dt(double s, double t)
{
    auto dSigma = [&](double M2)
    {
        double jacobian = 2. * sqrt(s) * _kinematics->pGamma(s) / PI;
        double result = invariant_xsection(s, t, M2);
        return jacobian * result;
    };

    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS61);
    ROOT::Math::Functor1D wF(dSigma);
    ig.SetFunction(wF);

    double xmin = _kinematics->M2_bounds(-1, s, t);
    double xmax = _kinematics->M2_bounds(+1, s, t);

    return ig.Integral(xmin, xmax);
};

double jpacPhoto::triple_regge::dsigma_dM2(double s, double M2)
{
    auto dSigma = [&](double t)
    {
        double jacobian = 2. * sqrt(s) * _kinematics->pGamma(s) / PI;
        double result = invariant_xsection(s, t, M2);
        return result;
    };

    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS61);
    ROOT::Math::Functor1D wF(dSigma);
    ig.SetFunction(wF);

    double xmin = _kinematics->t_bounds(-1, s, M2);
    double xmax = _kinematics->t_bounds(+1, s, M2);

    return ig.Integral(xmin, xmax);
};

// ---------------------------------------------------------------------------
// Singly-differential cross-sections in terms of (x, pT2)

// Integrate over x
double jpacPhoto::triple_regge::dsigma_dpT2(double s, double pT2)
{
    auto dSigma = [&](double x)
    {
        double M2 = _kinematics->M2_from_xpT2(s, x, pT2);
        double t  = _kinematics->t_from_xpT2(s, x, pT2);

        double jacobian = M_PI / sqrt(x*x + (_kinematics->_mX2 + pT2) / pow(_kinematics->pX_max(s), 2.));
        double result = jacobian * invariant_xsection(s, t, M2);
        
        return result;
    };

    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS61);
    ROOT::Math::Functor1D wF(dSigma);
    ig.SetFunction(wF);

    double xmin = 0.;
    double xmax = sqrt(1. - pT2 / pow(_kinematics->pX_max(s), 2.));

    return ig.Integral(xmin, xmax);
};

// Integrate over pT2
double jpacPhoto::triple_regge::dsigma_dx(double s, double x)
{
    auto dSigma = [&](double pT2)
    {
        double M2 = _kinematics->M2_from_xpT2(s, x, pT2);
        double t  = _kinematics->t_from_xpT2(s, x, pT2);

        double jacobian = M_PI / sqrt(x*x + (_kinematics->_mX2 + pT2) / pow(_kinematics->pX_max(s), 2.));
        double result = jacobian * invariant_xsection(s, t, M2);
        
        return result;
    };

    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS61);
    ROOT::Math::Functor1D wF(dSigma);
    ig.SetFunction(wF);
    
    double xmin = 0.;
    double xmax = pow(_kinematics->pX_max(s), 2.) * (1. - x*x);

    return ig.Integral(xmin, xmax);
};

// double jpacPhoto::triple_regge::integrated_xsection(double s)
// {
//     auto dSigma = [&](const double * in)
//     {
//         double x   = in[0], pT2 = in[1];

//         double M2 = _kinematics->M2(s, x, pT2);
//         double t  = _kinematics->t_man(s, x, pT2);

//         double jacobian = M_PI / sqrt(x*x + (_kinematics->_mX2 + pT2) / pow(_kinematics->pX_max(s), 2.));
//         double result = jacobian * invariant_xsection(s, t, M2);
        
//         return result;
//     };

//     // Integrate over costheta = [-1, 1] and x = [0., 1]   
//     double min[2] = { 0.,                  0.   };
//     double max[2] = { pow(pX_max(s), 2.),  M_PI };
    
//     ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kDEFAULT);
//     ROOT::Math::Functor wF(dSigma, 2);
//     ig.SetFunction(wF, 2);

//     double result = ig.Integral(min, max);

//     return result;
// };