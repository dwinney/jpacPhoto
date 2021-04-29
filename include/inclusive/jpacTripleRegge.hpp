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
            double tmin = _kinematics->t_bounds(+1, s, M2);
            double tprime = t - tmin;

            std::complex<double> regge_prop;
            regge_prop  = pow((s / M2), _trajectory->eval(t)); 

            switch (_debug)
            {
                case 0: 
                {
                    regge_prop *= xi(t);
                    regge_prop *= _coupling(t);
                    regge_prop *= cgamma(double(_trajectory->_minJ) - _trajectory->eval(t));
                    break;
                };
                case 1:
                {
                    regge_prop *= xi(tprime);
                    regge_prop *= _coupling(tprime);
                    regge_prop *= cgamma(double(_trajectory->_minJ) - _trajectory->eval(tprime));
                    break;
                };
                case 2:
                {
                    regge_prop *= xi(tprime);
                    regge_prop *= _coupling(t);
                    regge_prop *= cgamma(double(_trajectory->_minJ) - _trajectory->eval(tprime));
                    break;
                };
                case 3:
                {
                    regge_prop *= xi(t);
                    regge_prop *= _coupling(t);
                    regge_prop *= cgamma(double(_trajectory->_minJ) - _trajectory->eval(tprime));
                    break;
                };
                default: 
                {
                    regge_prop *= xi(t);
                    regge_prop *= _coupling(t);
                    regge_prop *= cgamma(double(_trajectory->_minJ) - _trajectory->eval(t));
                    break;
                };
            };
            
            double normalization = _trajectory->slope() / (16. * PI*PI*PI);
            double phase_space = sqrt(Kallen(M2, _kinematics->_mT2, t) / Kallen(s, _kinematics->_mT2, 0.));

            double result = normalization * phase_space * norm(regge_prop) * _sigmatot->operator()(s);
            return result;
        };

        int _debug = 0.;
        
        protected: 
        
        inclusive_kinematics * _kinematics;
        regge_trajectory * _trajectory;

        std::function<double(double)> _coupling;
        const sigma_tot * _sigmatot;

        constexpr static double _scale = 1.0;

        inline std::complex<double> xi(double t)
        {
            std::complex<double> signature_factor;
            signature_factor = 0.5 * (1. + double(_trajectory->_signature) * exp(XI * PI * _trajectory->eval(t)));

            return signature_factor;
        };
    };
};

#endif