// Form of the terms following JPAC's normalization
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef TRIPLE_JPAC
#define TRIPLE_JPAC

#include <complex>
#include <tuple>
#include <functional>
#include "misc_math.hpp"
#include "inclusive_kinematics.hpp"
#include "regge_trajectory.hpp"
#include "sigma_tot.hpp"

namespace jpacPhoto
{
    class jpacTripleRegge
    {
        public:

        // Fully general constructor
        jpacTripleRegge(inclusive_kinematics * xkinem, regge_trajectory * trajectory, 
                        const std::function<double(double)>& coupling, const sigma_tot * sigmatot)
        : _kinematics(xkinem), _trajectory(trajectory), _coupling(coupling), _sigmatot(sigmatot)
        {};
        
        inline double eval(double s, double t, double M2)
        {
            double result;
            result  = _sigmatot->operator()(s);
            result *= pow((s / M2), 2.*_trajectory->eval(t) - 1.); 
            result *= norm(xi(t));
            result *= _coupling(t) * _coupling(t);
            
            double normalization = _trajectory->slope() / (16. * PI*PI*PI);
            
            return normalization * result;
        };

        protected: 
        
        inclusive_kinematics * _kinematics;
        regge_trajectory * _trajectory;

        std::function<double(double)> _coupling;
        const sigma_tot * _sigmatot;

        constexpr static double _scale = 1.0;

        inline std::complex<double> xi(double t)
        {
            std::complex<double> signature_factor = 1., gamma_factor = 1.;
            signature_factor = 0.5 * (1. + double(_trajectory->_signature) * exp(XI * PI * _trajectory->eval(t)));
            gamma_factor = cgamma(double(_trajectory->_minJ) - _trajectory->eval(t));

            return signature_factor * gamma_factor;
        };
    };
};

#endif