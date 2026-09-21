// Simple s-channel Breit-Wigner amplitude for narrow baryon resonances decaying to Vector + proton
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef BARYON_BW_HPP
#define BARYON_BW_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"

namespace jpacPhoto
{
    class baryon_BW : public raw_amplitude
    {
        public: 

        // Each resonance is specified by: 
        // quantum numbers  J = |jp|/2, P = sign(jp)
        // pole parameters  mass = mw[0], width = mw[1]
        baryon_BW(key k, kinematics xkinem, int jp)
        : raw_amplitude(k, xkinem, id), 
          _resJ(abs(jp)), _resP(sgn(jp)),
        { 
            initialize(5); // 2 free parameters
            
            // J^P dependent quantities
            _naturality = _resP*pow(-1, (_resJ-1)/2);
            switch (jp)
            {
                case  1: {_lmin = 0; _pt = 2./3.; break;}
                case -1: {_lmin = 1; _pt = 3./5.; break;}
                case  3: {_lmin = 1; _pt = 3./5.; break;}
                case -3: {_lmin = 0; _pt = 2./3.; break;}
                case  5: {_lmin = 1; _pt = 3./5.; break;}
                case -5: {_lmin = 2; _pt = 1./3.; break;}
                default: warning("baryon_resonance", 
                                 "requested spin-parity combination not available!");
            };
        }; 

        inline complex helicity_amplitude(std::array<int,4> helicities, double s, double t)
        {
            // Save inputes
            store(helicities, s, t);

            int lam_i = 2 * helicities[0] - helicities[1];
            int lam_f = 2 * helicities[2] - helicities[3];

            // momentum for decay couplings
            _pibar = real(xkinem->initial_momentum(_mass*_mass));
            _pfbar = real(xkinem->final_momentum(_mass*_mass));

            std::complex<double> residue = 1.;
            residue  = photo_coupling(lam_i);
            residue *= hadronic_coupling(lam_f);
            residue *= threshold_factor(1.5);

            residue *= wigner_d_half(_resJ, lam_i, lam_f, _theta);
            residue /= (s-_mass*_mass+I*_mass*_width);

            return residue;
        };

        inline helicity_frame         native_helicity_frame(){ return S_CHANNEL; };
        inline std::vector<quantum_numbers> allowed_mesons() { return { VECTOR }; };
        inline std::vector<quantum_numbers> allowed_baryons(){ return { HALFPLUS }; };

        protected:

        // Three free parameters:
        // [0] the decay branching fraction to the final state
        // [1] decay constant of the final state vector F_V
        // [2] ratio of lambda = 1/2 to 3/2 photocoupling [0,1]
        inline void allocate_parameters(std::vector<double> pars)
        {
            _mass   = pars[0];
            _width  = pars[1];
            _xBR    = pars[2];
            _fv     = pars[3];
            _photoR = pars[4];
        };
        
        private:

        // Resonance quantum numbers 
        int _resJ = 0, _resP = +1, _naturality = +1, _lmin = 0;

        // Resonance pole parameters
        double _mass, _width;

        // Momenta at the resonance mass
        double _pibar, _pfbar;

        // Couplings which are free parameters
        double _xBR = 0, _photoR = 0, _fv = 0, _pt = 0;

        // Photoexcitation helicity amplitude for the process gamma p -> R
        inline complex photo_coupling(int lam_i)
        {
            // For spin-1/2 exchange no double flip
            if (abs(lam_i) > _resJ) return 0.;

            // A_1/2 or A_3/2 depending on ratio R_photo
            double a;
            (abs(lam_i) == 1) ? (a = _photoR) : (a = sqrt(1. - _photoR * _photoR));

            // Electromagnetic decay width given by VMD assumption
            complex emGamma = (_xBR*_width) * pow(_fv/_mX,2);
            emGamma *= pow(_pibar/_pfbar, 2*_lmin+1) * _pt;

            // Photo-coupling overall size of |A_1/2|^2 + |A_3/2|^2 is restriced from VMD
            complex A_lam = emGamma*PI*_mass*(_resJ + 1) / (2.*_mT*_pibar*_pibar);
            A_lam = csqrt(A_lam);

            complex result = sqrt(_s)*_pibar/_mass;
            result *= csqrt(8.*_mT*_mass/_kinematics->initial_momentum(_s));
            result *= A_lam * a;
            result *= sqrt(4. * PI * ALPHA);
            return result;
        };

        // Hadronic decay helicity amplitude for the R -> J/psi p process
        inline complex hadronic_coupling(int lam_f)
        {
            // Hadronic coupling constant g, given in terms of branching ratio xBR
            complex g;
            g  = 8. * PI * _xBR * _width;
            g *=  _mass*_mass *(_resJ + 1) / 6.;
            g /= pow(_pfbar, 2 * _lmin + 1);
            g = csqrt(XR * g);

            complex gpsi;
            gpsi = g * pow(_kinematics->final_momentum(_s), _lmin);
            if (lam_f < 0) gpsi *= _naturality;
            return gpsi;
        };

        // Ad-hoc threshold factor to kill the resonance at threshold
        inline double threshold_factor(double beta)
        {
            double result = pow((_s - _kinematics->sth()) / _s, beta);
            result /= pow((_mass*_mass - _kinematics->sth()) / (_mass*_mass), beta);

            return result;
        };        
    };
};

#endif 