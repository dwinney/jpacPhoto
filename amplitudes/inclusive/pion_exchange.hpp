//
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------

#ifndef PION_EXCHANGE_HPP       
#define PION_EXCHANGE_HPP

#include "constants.hpp"
#include "inclusive_process.hpp"
#include "total_xsection.hpp"
#include "sigma_tot/JPAC_piN.hpp"
#include "sigma_tot/PDG.hpp"

namespace jpacPhoto
{
    namespace fixed_spin
    {
        class pion_exchange : public raw_inclusive_process
        {
            public: 

            pion_exchange(inclusive_key key, double mX, int pm, std::string id = "")
            : raw_inclusive_process(key, mX, id),
              _sigma(new_total_xsection<JPAC_piN>(pm))
            {
                set_N_pars(2);
            };

            // Minimum mass is the proton 
            inline double minimum_M2(){ return pow(M_PROTON + M_PION, 2); };

            // Only free parameters are the top photocoupling and the form factor cutoff
            inline void allocate_parameters(std::vector<double> pars)
            {
                _g      = pars[0];
                _LamPi2 = pars[1]*pars[1];
            };

            // The invariant cross section used S T and M2 as independent vareiables
            inline double invariant_xsection(double s, double t, double M2)
            {
                double sigmatot    = _sigma->evaluate(M2, t) * 1E6; // in nb
                double phase_space = sqrt(Kallen(M2, t, M2_PROTON)/Kallen(s, 0., M2_PROTON));

                return pow(coupling(t)/(M2_PION - t), 2) * sigmatot * phase_space / (16*PI*PI*PI);
            };

            protected:

            inline double coupling(double t)
            {
                return exp((t - TMINfromM2(M2_PROTON))/_LamPi2) * (_g/sqrt(_mX2))*(_mX2 - t)/2;
            };

            private:
            
            total_xsection _sigma;
            double _g      = 0; // Top coupling
            double _LamPi2 = 0; // Exponential cut-off
        };
    };
};

#endif