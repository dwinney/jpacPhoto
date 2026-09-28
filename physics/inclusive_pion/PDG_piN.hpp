// Implementation of the HPR1R2 parameterization of hadronic total cross sections
// from the PDG 
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------

#ifndef PDG_PIN_HPP
#define PDG_PIN_HPP

#include "constants.hpp"

namespace jpacPhoto
{
    class PDG_piN
    {
        public: 

        PDG_piN()
        {};

        // Only available for on-shell beams so no q2 dependence
        double operator()(int iso, double s)
        {
            double sab = pow(M_PION + M_PROTON + _M, 2);
            return _delta*(_H*pow(log(s/sab), 2) + _P) + _R1*pow(s/sab, -_eta1) - iso*_R2*pow(s/sab, -_eta2);
        };

        private:

        // Process dependent constants
        double _delta   = 1;
        double _R1 = 9.56, _R2 = 1.767, _P = 18.75;

        // Process independent constants        
        double _M = 2.1206, _H = 0.2720, _eta1 = 0.4473, _eta2 = 0.5486;
    };
};

#endif