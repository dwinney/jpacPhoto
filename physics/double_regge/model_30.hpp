//  Model 3.0 according to the overleaf
//
// --------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
//               University of Bonn (HISKP)
// Email:        winney@hiskp.uni-bonn.de
// --------------------------------------------------------------------------------

#ifndef MODEL_30_HPP
#define MODEL_30_HPP

#include "constants.hpp"
#include "cgamma.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"

namespace jpacPhoto
{
    class model_30 : public two_meson::raw_amplitude
    {
        public: 

        // Basic constructor
        model_30(two_meson::amplitude_key key, two_meson::kinematics xkinem, std::string id = "phase_space")
        : two_meson::raw_amplitude(key, xkinem, id)
        {
            initialize(10);
        };

        inline complex helicity_amplitude(std::array<int,3> helicities, double s, double t, double s12, double thetaGJ, double phiGJ)
        {
            store(helicities, s, t, s12, thetaGJ, phiGJ);

            // Here labels 1 and 2 are 1 == Fast, 2 == Slow and NOT particle 1 and 2 as labeled in amplitude::kinematics
            // So depending on what _option is selected we filter the appropriate kinematic variables

            double alpha2 = _az2 + _t*_ap2, s1 = s12; // This always the same
            double alpha1, s2;
            if (_option == kFast1) // Particle 1 is fast 
            {
                alpha1 = _az1 + _t1*_ap1;
                s2 = _s2;
            }
            else // swap _t1 <-> _t2 and _s1 <-> _s2
            {
                alpha1 = _az1 + _t2*_ap1;
                s2 = _s1;
            }

            double V1 = _v1[0] + _v1[1]*s1/_s0 + _v1[2]*s2/_s0;
            double V2 = _v2[0] + _v2[1]*s1/_s0 + _v2[2]*s2/_s0;
            
            complex term1 = gamma(-alpha1) * pow(s/_s0, alpha1) * pow(s2/_s0, alpha2 - alpha1) * xi(-1, alpha1) * xi(+1, alpha2 - alpha1) * V1;
            complex term2 = gamma(-alpha2) * pow(s/_s0, alpha2) * pow(s1/_s0, alpha1 - alpha2) * xi(-1, alpha2) * xi(+1, alpha1 - alpha2) * V2;

            return term1 + term2;
        };
        
        inline void allocate_parameters(std::vector<double> pars)
        {
            _az1   = pars[0]; _ap1 = pars[1];
            _az2   = pars[2]; _ap2 = pars[3];
            _v1[0] = pars[4]; _v1[1] = pars[5]; _v1[2] = pars[6];
            _v2[0] = pars[7]; _v2[1] = pars[8]; _v2[2] = pars[9];
        }

        // Use amplitude::set_option to choose which particle is the fast one
        static const int kFast1 = 0;
        static const int kFast2 = 1;

        protected:
        
        // Trajectories
        double _az1, _az2, _ap1, _ap2;

        // Vertex factors
        std::array<double,3> _v1, _v2;

        // Scale parameter
        double _s0 = 1.; // GeV^2

        // Signature factor
        inline complex xi(int tau, double alpha)
        {
            return (tau + exp(-I*PI*alpha)) / 2.;
        };
    };
};

#endif