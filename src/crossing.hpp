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
#include "utilities.hpp"

namespace jpacPhoto 
{
    // This is a templated class to allow crossing to any channel but right now
    // we only use the S_CHANNEL, more crossing in the future
    template<helicity_frame OUT_CHANNEL>
    class crossed_to: public raw_amplitude
    {
        public: 
        
        crossed_to(key key, amplitude to_cross)
        : raw_amplitude(key, to_cross->get_kinematics(), to_cross->id()),
         _amplitude(to_cross)
        {
            initialize(to_cross->N_pars());
        };

        inline complex helicity_amplitude(std::array<int,4> helicities, double s, double t)
        {
            store(helicities, s,t);

            if (_amplitude->native_helicity_frame() == OUT_CHANNEL)
            {
                return _amplitude->helicity_amplitude(helicities, s, t);
            };

            if (_amplitude->native_helicity_frame() == helicity_frame::U_CHANNEL)
            {   
                return NaN<complex>();
            };

            if (!are_equal(_cached_s, s, _cache_tolerance) || !are_equal(_cached_t, t, _cache_tolerance))
            {
                _cached_amps = _amplitude->get_cache(s, t);
                _coswa = -cos_wigner_a(s,t); _coswb = -cos_wigner_b(s,t);
                _coswc = -cos_wigner_c(s,t); _coswd = -cos_wigner_d(s,t);
            };

            // Relabel helicities
            int Ja   = 1, Jb = 1, Jc = _kinematics->get_baryon_JP()[0], Jd = _kinematics->get_meson_JP()[0];
            int lama = helicities[1], lamb = helicities[0], lamc = helicities[3], lamd = helicities[2];
            complex sum = 0.;

            auto all_helicities = _amplitude->get_kinematics()->helicities();
            for (int i = 0; i < all_helicities.size(); i++)
            {
                auto help = all_helicities[i];
                int lampa = help[1], lampb = help[0], lampc = help[3], lampd = help[2];
                
                complex crossing_matrix;
                switch (OUT_CHANNEL)
                {
                    case (helicity_frame::S_CHANNEL):
                    {
                        crossing_matrix  = wigner_d_half_cos(Ja, lama, lampa, _coswa);
                        crossing_matrix *= wigner_d_int_cos( Jb, lamb, lampb, _coswb);
                        crossing_matrix *= wigner_d_half_cos(Jc, lamc, lampc, _coswc);
                        crossing_matrix *= wigner_d_int_cos( Jd, lamd, lampd, _coswd);
                        break;
                    };
                    case (helicity_frame::T_CHANNEL):
                    {
                        crossing_matrix  = wigner_d_half_cos(Ja, lampa, lama, _coswa);
                        crossing_matrix *= wigner_d_int_cos( Jb, lampb, lamb, _coswb);
                        crossing_matrix *= wigner_d_half_cos(Jc, lampc, lamc, _coswc);
                        crossing_matrix *= wigner_d_int_cos( Jd, lampd, lamd, _coswd);
                        break;
                    }
                };
                sum  += crossing_matrix*_cached_amps[i];
            };
            return sum;
        };

        // We now set the helicity_frame to S_CHANNEL
        inline helicity_frame native_helicity_frame(){ return OUT_CHANNEL; };

        // Everything else just gets piped through to the underlying amplitude
        inline std::vector<quantum_numbers> allowed_mesons() { return _amplitude->allowed_mesons();  };
        inline std::vector<quantum_numbers> allowed_baryons(){ return _amplitude->allowed_baryons(); };
        inline void set_option( int opt ){ _amplitude->set_option(opt); };
        inline std::vector<std::string> parameter_labels(){ return _amplitude->parameter_labels(); };

        protected:
        
        inline void allocate_parameters(std::vector<double> x){ _amplitude->allocate_parameters(x); };

        inline void store(std::array<int,4> helicities, double s, double t)
        {
            raw_amplitude::store(helicities, s, t);

            // We'll follow the particle labels in Martin & Spearman Ch. 7.2.2
            _ma2 = std::norm(_mT); // Target
            _mb2 = std::norm(_mB); // Beam
            _mc2 = std::norm(_mR); // Produced baryon
            _md2 = std::norm(_mX); // Produced meson
            _delta = _ma2 -_mb2 - _mc2 + _md2;
        };

        private:

        // Amplitude which is being crossed to the s-channel
        amplitude _amplitude;

        // Particle masses
        double _ma2, _mb2, _mc2, _md2, _delta;

        // Cosine of the crossing angles for each particule
        inline complex cos_wigner_a(double s, double t)
        {
            complex Pab = csqrt(kallen(_s, _ma2, _mb2));
            complex Tac = csqrt(kallen(_t, _ma2, _mc2));
            complex num = +(_s+_ma2-_mb2)*(_t+_ma2-_mc2)-2*_ma2*_delta;
            return num/Pab/Tac;
        };

        inline complex cos_wigner_b(double s, double t)
        {
            complex Pab = csqrt(kallen(_s, _ma2, _mb2));
            complex Tbd = csqrt(kallen(_t, _mb2, _md2));
            complex num = -(_s+_mb2-_ma2)*(_t+_mb2-_md2)-2*_mb2*_delta;
            return num/Pab/Tbd;
        };

        inline complex cos_wigner_c(double s, double t)
        {
            complex Pcd = csqrt(kallen(_s, _mc2, _md2));
            complex Tac = csqrt(kallen(_t, _ma2, _mc2));
            complex num = -(_s+_mc2-_md2)*(_t+_mc2-_ma2)-2*_mc2*_delta;
            return num/Pcd/Tac;
        };

        inline complex cos_wigner_d(double s, double t)
        {
            complex Pcd = csqrt(kallen(_s, _mc2, _md2));
            complex Tbd = csqrt(kallen(_t, _mb2, _md2));
            complex num = +(_s+_md2-_mc2)*(_t+_md2-_mb2)-2*_md2*_delta;
            return num/Pcd/Tbd;
        };

        // Cache for wigner angles
        complex _coswa, _coswb, _coswc, _coswd;
        std::vector<complex> _cached_amps;
    };

    template<helicity_frame C>
    amplitude cross_to(amplitude to_cross)
    {
        auto crossed_amp = std::make_shared<crossed_to<C>>(key(), to_cross);
        return crossed_amp;
    };
};

#endif 