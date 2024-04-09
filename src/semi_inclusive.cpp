// Extention of the amplitude structure but for inclusive reactions,
// here things are much easier since everything needs to be implements
// at the cross-section level not the amplitude level
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------

#include "semi_inclusive.hpp"

namespace jpacPhoto
{
    //----------------------------------------------------------------------------------------------
    // Method for summing terms together
    
    bool are_compatible(semi_inclusive a, semi_inclusive b)
    {
        bool same_kinem = (a->get_kinematics() == b->get_kinematics());
        
        std::string error_msg = "Attempted to add incompatible semi-inclusives: " + a->id() + " and " + b->id() + "!";
        if (!same_kinem) return error(error_msg + " (Contain different kinematics objects!)", false);

        return same_kinem;
    };  

    bool are_compatible(semi_inclusive a, amplitude b)
    {
        bool same_kinem = (a->get_kinematics() == b->get_kinematics());
        
        std::string error_msg = "Attempted to add incompatible objects: " + a->id() + " and " + b->id() + "!";
        if (!same_kinem) return error(error_msg + " (Contain different kinematics objects!)", false);
        
        return same_kinem;
    };  
    bool are_compatible(amplitude a, semi_inclusive b){ return are_compatible(b, a); };

    // If summing two semi-inclusives
    semi_inclusive operator+(semi_inclusive a, semi_inclusive b)
    {
        if (!are_compatible(a, b)) return nullptr;

        // New id will be the "sum" of the individual id's
        std::string id = a->id() + " + " + b->id(); 
        kinematics kinem = a->get_kinematics();

        // Get any potential sub diagrams
        std::vector<semi_inclusive> infrom_a = (a->is_sum()) ? a->_inclusives : std::vector<semi_inclusive>{{a}};
        std::vector<semi_inclusive> infrom_b = (b->is_sum()) ? b->_inclusives : std::vector<semi_inclusive>{{b}};
        infrom_a.insert(infrom_a.end(), infrom_b.begin(), infrom_b.end());

        std::vector<amplitude>      exfrom_a = (a->is_sum()) ? a->_exclusives : std::vector<amplitude>();
        std::vector<amplitude>      exfrom_b = (b->is_sum()) ? b->_exclusives : std::vector<amplitude>();
        exfrom_a.insert(exfrom_a.end(), exfrom_b.begin(), exfrom_b.end());

        return std::make_shared<raw_semi_inclusive>(key(), kinem, infrom_a, exfrom_a, id);
    };

    // If summing one semi-inclusie and an exclusive amplitude
    semi_inclusive operator+(semi_inclusive a, amplitude b)
    {
        if (!are_compatible(a, b)) return nullptr;

        // New id will be the "sum" of the individual id's
        std::string id = a->id(); 
        kinematics kinem = a->get_kinematics();

        // Get any potential sub diagrams
        std::vector<semi_inclusive> infrom_a = (a->is_sum()) ? a->_inclusives : std::vector<semi_inclusive>{{a}};
        
        std::vector<amplitude>      exfrom_a = (a->is_sum()) ? a->_exclusives : std::vector<amplitude>();
        std::vector<amplitude>      exfrom_b = extract_subamplitudes(b);
        exfrom_a.insert(exfrom_a.end(), exfrom_b.begin(), exfrom_b.end());

        return std::make_shared<raw_semi_inclusive>(key(), kinem, infrom_a, exfrom_a, id);
    };

    void operator+=(semi_inclusive a, amplitude b)
    {
        std::string error_msg = "Attempted to add incompatible objects: " + a->id() + " and " + b->id() + "!";
        if (!are_compatible(a, b)) return error("semi_inclusive::+=", error_msg + " (Contain different kinematics objects!)");
        
        a->_exclusives.push_back(b);
    };
    
    //----------------------------------------------------------------------------------------------
    // 
    double raw_semi_inclusive::invariant_xsection(double s, double t, double mm)
    {
        if (_inclusives.size() == 0) return std::nan("");

        double sum = 0;
        if (sqrt(s) >= sqrt(_mX2) + sqrt(minimum_M2()))
            for (semi_inclusive term : _inclusives) sum += term->invariant_xsection(s, t, mm);

        // Add any exclusive pole if applicable 
        if (_exclusives.size() > 0 && !use_TX() && are_equal(sqrt(mm), _kinematics->get_recoil_mass()))
        {
            amplitude amp = (_exclusives.size() > 1) ? _exclusives[0] + _exclusives[1] : _exclusives[0];
            for (int i = 2; i < _exclusives.size(); i++) amp += _exclusives[i];
            sum += amp->differential_xsection(s, t); 
        };

        return sum;
    };

    //----------------------------------------------------------------------------------------------
    // Parameter handling

    bool raw_semi_inclusive::correct_size(std::vector<double> pars)
    {
        if (pars.size() != _N_pars) return error(id()+"::set_parameters", "Number of parameters passed not the expected size!", false);
        return true;
    };

    void raw_semi_inclusive::set_parameters( std::vector<double> x )
    {
        // Check the size is alright
        if (!correct_size(x)) return;

        // Allocate them (amplitude specific)
        allocate_parameters(x);
    };

    // Momentum of the photon
    double raw_semi_inclusive::qGamma(double s)
    {
        return (s - _mT2) / sqrt(4*s);
    };

    // Max momentum of the produced particle
    double raw_semi_inclusive::pMax(double s)
    {
        return sqrt(Kallen(s, _mX2, minimum_M2())) / (2*sqrt(s));
    };

    // ---------------------------------------------------------------------------
    // POLAR VARIABLES (r, cos)
    
    // Energy of the particle X
    double raw_semi_inclusive::EfromRCOS(double s, double r, double cos)
    {
        return sqrt(_mX2 + pMax(s)*pMax(s)*r*r);
    };

    // Jacobian in polar coordinates
    double raw_semi_inclusive::jacobianRCOS(double s, double r, double cos)
    {
        return(2*PI*r*r*pow(pMax(s), 3.)) /  EfromRCOS(s, r, cos);
    };

    // ---------------------------------------------------------------------------
    // CARTESIAN VARIABLES (x, y)

    // Energy of X in (x, y)
    double raw_semi_inclusive::EfromXY(double s, double x, double y)
    {
        return sqrt(_mX2 + pMax(s)*pMax(s) * (x*x + y*y));
    };

    // Jacobian in (x, y)
    double raw_semi_inclusive::jacobianXY(double s, double x, double y)
    {
        return (2.* M_PI * y * pow(pMax(s), 3)) / EfromXY(s, x, y);
    };

    // Energy of X in (x, y2)
    double raw_semi_inclusive::EfromXY2(double s, double x, double y2)
    {
        return sqrt(_mX2 + pMax(s)*pMax(s) * (x*x + y2));
    };

    // Jacobian in (x, y2)
    double raw_semi_inclusive::jacobianXY2(double s, double x, double y2)
    {
        return (M_PI * pow(pMax(s), 3)) / EfromXY2(s, x, y2);
    };

    double raw_semi_inclusive::XfromRCOS(double s, double r, double cos)
    {
        return r * cos;
    };

    double raw_semi_inclusive::XfromTM2(double s, double t, double M2)
    {
        double r   = RfromTM2(s, t, M2);
        double cos = COSfromTM2(s, t, M2);
        return XfromRCOS(s, r, cos);
    };

    // ---------------------------------------------------------------------------
    // INVARIANT variables (t, M2)

    // momentum of produced at a fixed missing mass
    double raw_semi_inclusive::pXfromM2(double s, double M2)
    {
        return sqrt(Kallen(s, _mX2, M2)) / (2*sqrt(s));
    };  

    double raw_semi_inclusive::COSfromTM2(double s, double t, double M2)
    {
        double u = _mT2 + _mX2 + M2 - s - t;
        double num  = s * (t - u) - _mT2*(_mX2 - M2);

        return num / (4*_s*qGamma(s)*pXfromM2(s, M2));
    };

    double raw_semi_inclusive::RfromTM2(double s, double t, double M2)
    {
        return pXfromM2(s, M2) / pMax(s);
    };

    // Jacobian in (t, M2)
    double raw_semi_inclusive::jacobianTM2(double s, double t, double M2)
    {
        return PI/ (2*sqrt(s)*qGamma(s));
    };

    // t from cos and M2 
    double raw_semi_inclusive::TfromM2COS(double s, double M2, double cos)
    {
        return _mX2 - (s - _mT2) * (s - M2 + _mX2) / (2*s) + 2*qGamma(s)*pXfromM2(s, M2)*cos;
    };

    double raw_semi_inclusive::M2fromRCOS(double s, double r, double cos)
    {
        return _mX2 + s - 2*sqrt(_mX2*s + pMax(s)*pMax(s)*r*r*s);
    };

    double raw_semi_inclusive::TfromRCOS(double s, double r, double cos)
    {
        double M2 = M2fromRCOS(s, r, cos);
        return _mX2 - (s - _mT2) * (s - M2 + _mX2) / (2*s) + 2*qGamma(s)*pMax(s)*r*cos;
    };

    // Bounds of integration in t at fixed M2
    double raw_semi_inclusive::TMINfromM2(double s, double M2) { return TfromM2COS(s, M2,  1); };
    double raw_semi_inclusive::TMAXfromM2(double s, double M2) { return TfromM2COS(s, M2, -1); };

    // Bounds of integration in M2 at fixed t
    double raw_semi_inclusive::M2MINfromT(double s, double t) { return minimum_M2(); };
    double raw_semi_inclusive::M2MAXfromT(double s, double t) 
    {
        return (_mT2 + _mX2 - s - t) * (_mT2 * _mX2 - s * t) / (_mT2 - s) / (_mX2 - t);
    };

    // ---------------------------------------------------------------------------
    // MIXED variables (t, x)

    // Missing mass from the cartesian XY
    double raw_semi_inclusive::M2fromXY(double s, double x, double y)
    {
        return s + _mX2 - 2*sqrt(s)*EfromXY(s, x, y);
    }; 
    double raw_semi_inclusive::M2fromXY2(double s, double x, double y2)
    {
        return s + _mX2 - 2*sqrt(s)*EfromXY(s, x, sqrt(y2));
    }; 

    // Similarly, momentum transfer from XY
    double raw_semi_inclusive::TfromXY(double s, double x, double y)
    {
        return _mX2 - 2*qGamma(s)*(EfromXY(s, x, y) - pMax(s)*x);
    };
    double raw_semi_inclusive::TfromXY2(double s, double x, double y2){ return TfromXY(s, x, sqrt(y2)); };

    // Bounds of integration for T at fixed X
    double raw_semi_inclusive::TMINfromX(double s, double x) { return TMINfromM2( s, M2fromXY(s, x, sqrt(1.-x*x))); };
    double raw_semi_inclusive::TMAXfromX(double s, double x) { return TMAXfromM2( s, M2fromXY(s, x, 0)); };

    // Jacobian for mixed variables
    double raw_semi_inclusive::jacobianTX(double s, double t, double x)
    {
        double inv = qGamma(s) / M_PI / pMax(s);
        return 1./inv;
    };

    // Also useful is M2 as a function of X and T
    double raw_semi_inclusive::M2fromTX(double s, double t, double x)
    {
        double lami = Kallen(s, 0., _mT2);
        double lamf = Kallen(s, _mX2, minimum_M2());
        double num = _mT2 * _mX2 + _mT2 * s + _mX2 * s - s*s - 2*s*t + sqrt(lami*lamf)*x;
        return num / (_mT2 - s);
    };

    // ---------------------------------------------------------------------------
    // Differential cross-sections
    // OUTPUT IN NANOBARN

    double raw_semi_inclusive::dsigma_dtdM2(double s, double t, double M2)
    {
        double mx = (use_TX()) ? XfromTM2(s, t, M2) : M2; 
        return jacobianTM2(s, t, M2) * invariant_xsection(s, t, mx);
    };

    double raw_semi_inclusive::dsigma_dtdx(double s, double t, double x)
    {
        double mx = (use_TX()) ? x : M2fromTX(s, t,x); 
        return jacobianTX(s, t, x) * invariant_xsection(s, t, mx);
    };

    double raw_semi_inclusive::dsigma_dxdy2(double s, double x, double y2)
    {
        double t  = TfromXY2(s, x ,y2);
        double mx = (use_TX()) ? x : M2fromXY2(s, t, x); 
        return jacobianXY2(s, x, y2) * invariant_xsection(s, t, mx);
    }

    // ---------------------------------------------------------------------------
    // Singly-differential cross-sections in terms of (t, M2)

    double raw_semi_inclusive::dsigma_dt(double s, double t)
    {
        gErrorIgnoreLevel = 6001;

        ROOT::Math::GSLIntegrator ig( ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS15);
        ROOT::Math::Functor1D wF( [&](double M2){ return dsigma_dtdM2(s, t, M2); } );
        ig.SetFunction(wF);
        double total =  ig.Integral(M2MINfromT(s, t), M2MAXfromT(s, t));

        if (_exclusives.size() != 0)
        {
            for (amplitude amp : _exclusives) total += amp->differential_xsection(s, t);
        };

        return total;
    };

    double raw_semi_inclusive::dsigma_dM2(double s, double M2)
    {
        gErrorIgnoreLevel = 6001;

        auto dSigma = [&](double t)
        {
            return dsigma_dtdM2(s, t, M2);
        };
        ROOT::Math::GSLIntegrator ig( ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS15);
        ROOT::Math::Functor1D wF(dSigma);
        ig.SetFunction(wF);

        return ig.Integral(TMAXfromM2(s, M2), TMINfromM2(s, M2));
    };

    // ---------------------------------------------------------------------------
    // Singly-differential cross-sections with respect to fixed (x, y2)
    // TODO: The exclusive reaction pole might not be well implemented here

    double raw_semi_inclusive::dsigma_dy2(double s, double y2)
    {
        gErrorIgnoreLevel = 6001;

        auto dSigma = [&](double x)
        {
            return dsigma_dxdy2(s, x, y2);
        };
        ROOT::Math::GSLIntegrator ig( ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS15);
        ROOT::Math::Functor1D wF(dSigma);
        ig.SetFunction(wF);

        return ig.Integral(0, sqrt(1 - y2));
    };

    double raw_semi_inclusive::dsigma_dx(double s, double x)
    {
        gErrorIgnoreLevel = 6001;

        auto dSigma = [&](double y2)
        {
            return dsigma_dxdy2(s, x, y2);
        };
        ROOT::Math::GSLIntegrator ig( ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS15);
        ROOT::Math::Functor1D wF(dSigma);
        ig.SetFunction(wF);

        return ig.Integral(0, 1 - x*x);
    };

    // ---------------------------------------------------------------------------
    // OUTPUT IN NANOBARN
    double raw_semi_inclusive::integrated_xsection(double s, double cos_min)
    {
        gErrorIgnoreLevel = 6001;

        double total = 0;
        if (sqrt(s) >= sqrt(_mX2) + sqrt(minimum_M2()))
        {
            auto dSigma = [&](const double * rcos)
            {
                double r = rcos[0], cos = rcos[1];
                double t  = TfromRCOS(s, r, cos);
                double mx = (use_TX()) ? XfromRCOS(s, r, cos) : M2fromRCOS(s, r, cos);
                return jacobianRCOS(s, r, cos) * invariant_xsection(s, t, mx); 
            };
            ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::Type::kVEGAS, 1E-5, 1E-5);
            double min[2] = {0., cos_min}, max[2] = {1., 1.};
            
            total = ig.Integral(dSigma, 2, min, max);
        };
        
        if (_exclusives.size() > 0)
        {
            amplitude amp = (_exclusives.size() > 1) ? _exclusives[0] + _exclusives[1] : _exclusives[0];
            for (int i = 2; i < _exclusives.size(); i++) amp += _exclusives[i];
            total += amp->integrated_xsection(s); 
        };

        return total;
    };
};