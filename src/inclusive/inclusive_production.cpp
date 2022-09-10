// Abstract class for the invariant cross-section from a triple regge interaction.
// Contains inclusve_kinematics objects as well as dynamical objects with either
// 'JPAC' or 'Field & Fox' models.
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "inclusive_production.hpp"

// ---------------------------------------------------------------------------
// Differential cross-sections
// OUTPUT IN NANOBARN


double jpacPhoto::inclusive_production::dsigma_dtdM2(double s, double t, double M2)
{
    if (sqrt(s) <= _kinematics->_mX + _kinematics->_mT) return 0.;

    // Pass the total energy to the kinematics object
    _kinematics->_s = s; 

    double mx;
    (_useTX) ? (mx = _kinematics->XfromTM2(t, M2)) : (mx = M2); 

    return _kinematics->jacobianTM2(t, M2) * d3sigma_d3p(s, t, mx);
}

double jpacPhoto::inclusive_production::dsigma_dtdx(double s, double t, double x)
{
    if (sqrt(s) <= _kinematics->_mX + _kinematics->_mT) return 0.;

    // Pass the total energy to the kinematics object
    _kinematics->_s = s; 

    double mx;
    (_useTX) ? (mx = x) : (mx = _kinematics->M2fromTX(t, x)); 

    return _kinematics->jacobianTX(t, x) * d3sigma_d3p(s, t, mx);
}

double jpacPhoto::inclusive_production::dsigma_dxdy2(double s, double x, double y2)
{
    if (sqrt(s) <= _kinematics->_mX + _kinematics->_mT) return 0.;

    // Pass the total energy to the kinematics object
    _kinematics->_s = s; 

    double t = _kinematics->TfromXY2(x ,y2);
    double mx;
    (_useTX) ? (mx = x) : (mx = _kinematics->M2fromXY2(t, x)); 

    return _kinematics->jacobianXY2(x, y2) * d3sigma_d3p(s, t, mx);
}

// ---------------------------------------------------------------------------
// Singly-differential cross-sections in terms of (t, M2)
double jpacPhoto::inclusive_production::dsigma_dt(double s, double t)
{
    if (sqrt(s) <= _kinematics->_mX + _kinematics->_mT) return 0.;

    // Pass the total energy to the kinematics object
    _kinematics->_s = s;

    // How we integrate depends on what variables are being used
    double result = 0.;
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS15);

    // Integrate the differential cross-section
    auto dSigma = [&](double M2)
    {
        return dsigma_dtdM2(s, t, M2);
    };
    ROOT::Math::Functor1D wF(dSigma);
    ig.SetFunction(wF);

    result = ig.Integral(_kinematics->M2MINfromT(t), _kinematics->M2MAXfromT(t));

    return result;
};

double jpacPhoto::inclusive_production::dsigma_dM2(double s, double M2)
{
    if (sqrt(s) <= _kinematics->_mX + sqrt(_kinematics->_minM2)) return 0.;
    if (!_useTX && (sqrt(M2) >= sqrt(s) - _kinematics->_mX)) return 0.;

    // Pass the total energy to the kinematics object
    _kinematics->_s = s;

    // How we integrate depends on what variables are being used
    double result = 0.;
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS15);

    // integrate over usual t bounds
    auto dSigma = [&](double t)
    {
        return dsigma_dtdM2(s, t, M2);
    };
    ROOT::Math::Functor1D wF(dSigma);
    ig.SetFunction(wF);

    result = ig.Integral(_kinematics->TMAXfromM2(M2), _kinematics->TMINfromM2(M2));

    return result;
};

// ---------------------------------------------------------------------------
// Singly-differential cross-sections with respect to fixed (x, y2)
double jpacPhoto::inclusive_production::dsigma_dy2(double s, double y2)
{
    if (sqrt(s) <= _kinematics->_mX + _kinematics->_mT) return 0.;

    // Pass the total energy to the kinematics object
    _kinematics->_s = s;

    // How we integrate depends on what variables are being used
    double result = 0.;
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS15);

    // Argument 3 is x and we dont need to convert
    auto dSigma = [&](double x)
    {
        return dsigma_dxdy2(s, x, y2);
    };
    ROOT::Math::Functor1D wF(dSigma);
    ig.SetFunction(wF);

    result = ig.Integral(0., sqrt(1. - y2));

    return result;
};

// Cross-section with fixed x
double jpacPhoto::inclusive_production::dsigma_dx(double s, double x)
{
    if (sqrt(s) <= _kinematics->_mX + sqrt(_kinematics->_minM2)) return 0.;

    // Pass the total energy to the kinematics object
    _kinematics->_s = s;

    // How we integrate depends on what variables are being used
    double result = 0.;
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS15);

    auto dSigma = [&](double y2)
    {
        return dsigma_dxdy2(s, x, y2);
    };
    ROOT::Math::Functor1D wF(dSigma);
    ig.SetFunction(wF);

    result = ig.Integral(0., 1. - x*x);
   
    return result;
};

// ---------------------------------------------------------------------------
// OUTPUT IN NANOBARN
double jpacPhoto::inclusive_production::integrated_xsection(double s)
{
    if (sqrt(s) <= _kinematics->_mX + sqrt(_kinematics->_minM2)) return 0.;

    // Pass the total energy to the kinematics object
    _kinematics->_s = s;

    double result = 0.;

    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS15);
    
    if (_useTX == true)
    {
        // Assume argument 3 is M2
        auto dSigma = [&](double x)
        {
            return dsigma_dx(s,x);
        };
        ROOT::Math::Functor1D wF(dSigma);
        ig.SetFunction(wF);

        // Integrate over x
        result  = ig.Integral(0., 1.);
    }
    else
    {
        // Assume argument 3 is M2
        auto dSigma = [&](double m2)
        {
            return dsigma_dM2(s, m2);
        };
        ROOT::Math::Functor1D wF(dSigma);
        ig.SetFunction(wF);

        // Integrate over x
        result = ig.Integral(_kinematics->_minM2, pow(sqrt(s) - _kinematics->_mX, 2));
    }

    return result;
};