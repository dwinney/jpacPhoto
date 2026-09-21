// Methods which implement crossing matrix numerically
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2026)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Instituto de Ciencias Nucleares (ICN)
// Email:        daniel.winney@nucleares.unam.mx
// ------------------------------------------------------------------------------

#include "crossing.hpp"

namespace jpacPhoto 
{
    // --------------------------------------------------------------------------
    // Cosine of the crossing angles for each particle
    complex crossed_to_schannel::wigner_cos_a(double s, double t)
    {
        complex Pab = csqrt(kallen(_s, _ma2, _mb2));
        complex Tac = csqrt(kallen(_t, _ma2, _mc2));
        complex num = -(_s+_ma2-_mb2)*(_t+_ma2-_mc2)-2*_ma2*_delta;
        return num/Pab/Tac;
    };

    complex crossed_to_schannel::wigner_cos_b(double s, double t)
    {
        complex Pab = csqrt(kallen(_s, _ma2, _mb2));
        complex Tbd = csqrt(kallen(_t, _mb2, _md2));
        complex num = +(_s+_mb2-_ma2)*(_t+_mb2-_md2)-2*_mb2*_delta;
        return num/Pab/Tbd;
    };

    complex crossed_to_schannel::wigner_cos_c(double s, double t)
    {
        complex Pcd = csqrt(kallen(_s, _mc2, _md2));
        complex Tac = csqrt(kallen(_t, _ma2, _mc2));
        complex num = +(_s+_mc2-_mb2)*(_t+_mc2-_ma2)-2*_mc2*_delta;
        return num/Pcd/Tac;
    };

    complex crossed_to_schannel::wigner_cos_d(double s, double t)
    {
        complex Pcd = csqrt(kallen(_s, _mc2, _md2));
        complex Tbd = csqrt(kallen(_t, _mb2, _md2));
        complex num = -(_s+_md2-_mc2)*(_t+_md2-_mb2)-2*_md2*_delta;
        return num/Pcd/Tbd;
    };

    // --------------------------------------------------------------------------

    complex crossed_to_schannel::helicity_amplitude(std::array<int,4> helicities, double s, double t)
    {
        if (_amplitude->native_helicity_frame() == helicity_frame::S_CHANNEL)
        {
            return _amplitude->helicity_amplitude(helicities, s, t);
        };

        if (_amplitude->native_helicity_frame() != helicity_frame::T_CHANNEL)
        {   
            return NaN<complex>();
        };

        if (!are_equal(_cached_s, s, _cache_tolerance) || !are_equal(_cached_t, t, _cache_tolerance))
        {
            _cosa = wigner_cos_a(s,t); _cosb = wigner_cos_b(s,t);
            _cosc = wigner_cos_c(s,t); _cosd = wigner_cos_d(s,t);
        };
        
        // Relabel helicities
        int lama = helicities[1], lamb = helicities[0], lamc = helicities[3], lamd = helicities[2];

        complex sum = 0.;
        for (auto help : _amplitude->get_kinematics()->helicities())
        {
            int lampa = help[1], lampb = help[0], lampc = help[3], lampd = help[2];
            // sum += 
        };

    };
};
