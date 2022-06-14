// Abstract class for the invariant cross-section from a triple regge interaction.
// Contains inclusve_kinematics objects as well as dynamical objects with either
// 'JPAC' or 'Field & Fox' models.
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "inclusive/inclusive_production.hpp"

// ---------------------------------------------------------------------------
// Singly-differential cross-sections in terms of (t, M2)
double jpacPhoto::inclusive_production::dsigma_dt(double s, double t)
{
    if (sqrt(s) <= _kinematics->_mX + _kinematics->_mT) return 0.;

    // Pass the total energy to the kinematics object
    _kinematics->_s = s;

    // How we integrate depends on what variables are being used
    double result = 0.;
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS61);

    if (_useTX == false)
    {
        // Assume argument 3 is M2
        auto dSigma = [&](double M2)
        {
            return _kinematics->jacobianTM2(t, M2) * d3sigma_d3p(s, t, M2);
        };
        ROOT::Math::Functor1D wF(dSigma);
        ig.SetFunction(wF);

        result = ig.Integral(_kinematics->M2MINfromT(t), _kinematics->M2MAXfromT(t));
    }
    else
    {
        // TODO: add integration over x at fixed t
        result = 0.;
    };

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
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS61);

    if (_useTX == false)
    {
        // Assume argument 3 is M2
        // integrate over usual t bounds
        auto dSigma = [&](double t)
        {
            return _kinematics->jacobianTM2(t, M2) * d3sigma_d3p(s, t, M2);
        };
        ROOT::Math::Functor1D wF(dSigma);
        ig.SetFunction(wF);

        result = ig.Integral(_kinematics->TMAXfromM2(M2), _kinematics->TMINfromM2(M2));
    }
    else
    {
        // Argument 3 is x
        // Convert the M2 taken as argument to an x in the high-energy limit
        double x = 1. - M2 / s;

        auto dSigma = [&](double t)
        {
            return _kinematics->jacobianTX(t, x) * d3sigma_d3p(s, t, x);
        };
        ROOT::Math::Functor1D wF(dSigma);
        ig.SetFunction(wF);

        result = ig.Integral(_kinematics->TMINfromX(x), _kinematics->TMAXfromX(x));
    };

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
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS61);

    if (_useTX == false)
    {
        // Assume argument 3 is M2
        auto dSigma = [&](double x)
        {
            // convert x and y2 to fixed M2
            double M2 = _kinematics->M2fromXY2(x, y2);
            double t  = _kinematics->TfromXY2(x, y2);
            return _kinematics->jacobianXY2(t, M2) * d3sigma_d3p(s, t, M2);
        };
        ROOT::Math::Functor1D wF(dSigma);
        ig.SetFunction(wF);

        result = ig.Integral(0., sqrt(1. - y2));
    }
    else
    {
        // Argument 3 is x and we dont need to convert
        auto dSigma = [&](double x)
        {
            double t  = _kinematics->TfromXY2(x, y2);
            return _kinematics->jacobianXY2(x, y2) * d3sigma_d3p(s, t, x);
        };
        ROOT::Math::Functor1D wF(dSigma);
        ig.SetFunction(wF);

        result = ig.Integral(0., sqrt(1. - y2));
    };

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
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS61);
    
    if (_useTX == false)
    {
        // Assume argument 3 is M2
        auto dSigma = [&](double y2)
        {
            // convert x and y2 to fixed M2
            double M2 = _kinematics->M2fromXY2(x, y2);
            double t  = _kinematics->TfromXY2(x, y2);
            return _kinematics->jacobianXY2(x, y2) * d3sigma_d3p(s, t, M2);
        };
        ROOT::Math::Functor1D wF(dSigma);
        ig.SetFunction(wF);

        result = ig.Integral(0., 1. - x*x);
    }
    else
    {
        // Argument 3 is x and we dont need to convert
        auto dSigma = [&](double y2)
        {
            double t  = _kinematics->TfromXY2(x, y2);
            return _kinematics->jacobianXY2(x, y2) * d3sigma_d3p(s, t, x);
        };
        ROOT::Math::Functor1D wF(dSigma);
        ig.SetFunction(wF);

        result = ig.Integral(0., 1. - x*x);
    };

    return result;
};

// ---------------------------------------------------------------------------
double jpacPhoto::inclusive_production::integrated_xsection(double s)
{
    if (sqrt(s) <= _kinematics->_mX + sqrt(_kinematics->_minM2)) return 0.;

    // Pass the total energy to the kinematics object
    _kinematics->_s = s;

    double result = 0.;

    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS61);
    
    if (_useTX == true)
    {
        // Assume argument 3 is M2
        auto dSigma = [&](double x)
        {
            double sigma = dsigma_dx(s,x);
            return sigma;
        };
        ROOT::Math::Functor1D wF(dSigma);
        ig.SetFunction(wF);

        // Integrate over x
        result  = ig.Integral(0.1, 1.);
    }
    else
    {
        // Assume argument 3 is M2
        auto dSigma = [&](double m2)
        {
            double sigma = dsigma_dM2(s,m2);
            return sigma;
        };
        ROOT::Math::Functor1D wF(dSigma);
        ig.SetFunction(wF);

        // Integrate over x
        result = ig.Integral(_kinematics->_minM2, pow(sqrt(s) - _kinematics->_mX, 2));
    }

    return result * 1.E3;
};