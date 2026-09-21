// Methods which implement crossing matrix numerically
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2026)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Instituto de Ciencias Nucleares (ICN)
// Email:        daniel.winney@nucleares.unam.mx
// ------------------------------------------------------------------------------

#ifndef CROSSING_HPP
#define CROSSING_HPP

#include "kinematics.hpp"
#include "amplitude.hpp"

namespace jpacPhoto 
{
    class crossed_to_schannel: public raw_amplitude
    {
        public: 
        
        crossed_to_schannel(key key, amplitude to_cross)
        : raw_amplitude(key, to_cross->get_kinematics(), to_cross->id()),
        _amplitude(to_cross)
        {
            initialize(to_cross->N_pars());

            // We'll follow the particle labels in Martin & Spearman Ch. 7.2.2
            _ma2 = std::norm(_mT); // Target
            _mb2 = std::norm(_mB); // Beam
            _mc2 = std::norm(_mR); // Produced baryon
            _md2 = std::norm(_mX); // Produced meson
            _delta = _mb2 - _md2 - _ma2 - _mc2;
        };

        complex helicity_amplitude(std::array<int,4> helicities, double s, double t);

        // We now set the helicity_frame to S_CHANNEL
        inline helicity_frame native_helicity_frame(){ return S_CHANNEL; };

        // Everything else just gets piped through to the underlying amplitude
        inline std::vector<quantum_numbers> allowed_mesons() { return _amplitude->allowed_mesons();  };
        inline std::vector<quantum_numbers> allowed_baryons(){ return _amplitude->allowed_baryons(); };
        inline void set_option( int opt ){ _amplitude->set_option(opt); };
        inline std::vector<std::string> parameter_labels(){ return _amplitude->parameter_labels(); };

        protected:
        
        inline void allocate_parameters(std::vector<double> x){ _amplitude->allocate_parameters(x); };

        private:

        // Amplitude which is being crossed to the s-channel
        amplitude _amplitude;

        // Particle masses
        double _ma2, _mb2, _mc2, _md2, _delta;

        // Cosine of the crossing angles for each particule
        complex wigner_cos_a(double s, double t);
        complex wigner_cos_b(double s, double t);
        complex wigner_cos_c(double s, double t);
        complex wigner_cos_d(double s, double t);

        // Cache for wigner angles
        complex _cosa, _cosb, _cosc, _cosd;
    };
};

#endif 