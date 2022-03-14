// Phenomenological expressions for the total cross-sections.
// We use a generic class callable by double sigma_tot(double) to select different
// parameterizations or reactions
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef SIGMA_TOT
#define SIGMA_TOT

#include "constants.hpp"

#include <vector>
#include <array>
#include <fstream>
#include <sstream>

namespace jpacPhoto
{
    
    // ---------------------------------------------------------------------------
    // Generic template class

    class sigma_tot
    {
        public:

        // Default constructor 
        sigma_tot(){};
        
        // Only thing that is needed is a way to evaluate the cross-section
        virtual double operator()(double s) = 0;
    };

    // ---------------------------------------------------------------------------
    // Coded up PDG parameterization for the asymptotic cross-sections

    class sigma_tot_PDG : public sigma_tot
    {
        public:

        // Constructor with masses and a filename to look for data
        sigma_tot_PDG(double m1, double m2, std::array<double, 5> pdgparams)
        : _mBeam(m1), _mTarget(m2), _threshold((m1+m2)*(m1+m2))
        {
            _iso   = pdgparams[0];
            _delta = pdgparams[1];
            _R1    = pdgparams[2];
            _R2    = pdgparams[3];
            _P     = pdgparams[4];  
        };

        double operator()(double s)
        {
            double result = 0.;
            if (s < _threshold + 10.*EPS) 
            {
                return 0.;
            } 
            else 
            {
                result = PDG_parameterization(s);
            } 
            return result *= 1.E6; // convet mb -> nb
        };

        private:

        double _mBeam, _mTarget;
        double _threshold;

        // Else use the fits from the PDG 
        inline double PDG_parameterization(double s)
        {
            // These params are universal for all collisions
            double M = 2.1206, H = 0.2720, eta1 = 0.4473, eta2 = 0.5486;

            double sab;
            if (_mBeam < 1.E-4) sab = pow(M_RHO + _mTarget + M, 2);
            else sab = pow(_mBeam + _mTarget + M, 2);

            return _delta*(H*pow(log(s/sab), 2.) + _P) + _R1*pow(s/sab, -eta1) - _iso*_R2*pow(s/sab, -eta2);
        };

        int _iso = 0;
        double _delta = 1.;
        double _R1, _R2, _P;
    };

        // Pi+ Proton 
    sigma_tot_PDG sigma_tot_pipp(M_PION, M_PROTON,  {+1., 1., 9.56, 1.767, 18.75});

    // Pi- Proton
   sigma_tot_PDG sigma_tot_pimp(M_PION, M_PROTON,  {-1., 1., 9.56, 1.767, 18.75});

    // Gamma Proton
   sigma_tot_PDG sigma_tot_gamp(0.,     M_PROTON,  { 0., 3.065E-3, 0.0139, 0., 34.41});
};

#endif