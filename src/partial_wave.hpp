// Extension of the raw_amplitude which supports individual partial waves in the s-channel
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef PARTIAL_WAVE_HPP
#define PARTIAL_WAVE_HPP

#include "utilities.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"
#include "crossing.hpp"

#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <memory>

namespace jpacPhoto
{
    // Foward declare the PW amplitude
    class raw_partial_wave;

    // Similar to amplitude we only ever want partial_waves to be pointers
    using partial_wave = std::shared_ptr<raw_partial_wave>;

    template<class A>
    inline partial_wave new_partial_wave(kinematics xkinem, uint J)
    {
        auto amp = std::make_shared<A>(key(), xkinem, J);
        return std::static_pointer_cast<raw_partial_wave>(amp);
    };

    template<class A, class B>
    inline partial_wave new_partial_wave(kinematics xkinem, uint J, B extra)
    {
        auto amp = std::make_shared<A>(key(), xkinem, J, extra);
        return std::static_pointer_cast<raw_partial_wave>(amp);
    };

    // Summing two partial waves defaults to a "full" amplitude
    inline amplitude operator+(partial_wave a, partial_wave b)
    {
        amplitude wa = std::static_pointer_cast<raw_amplitude>(a);
        amplitude wb = std::static_pointer_cast<raw_amplitude>(b);
        
        return wa + wb;
    };

    // ---------------------------------------------------------------------------
    // Raw_amplitude class

    class raw_partial_wave : public raw_amplitude
    {
        public:

        // This constructor should be used for any user defined derived classes
        raw_partial_wave(key key, kinematics xkinem,  int J, std::string id = "partial_wave")
        : raw_amplitude(key, xkinem, id), _J(J)
        {
            set_N_pars(0);
        };

        // Return the J-th term to the full amplitude by multiplying by angular function
        // These may be overloaded with an explicit model for the PWA
        complex helicity_amplitude(std::array<int,4> helicities, double s, double t)
        {
            // s-channel scattering angle
            store(helicities, s, t);

            switch (this->native_helicity_frame())
            {
                case helicity_frame::S_CHANNEL: 
                {
                    // Net helicities
                    int lam  = 2 * helicities[0] - helicities[1]; // Photon - Target
                    int lamp = 2 * helicities[2] - helicities[3]; // Meson  - Recoil
                    
                    return (_J + 1) * wigner_d_half(_J, lam, lamp, _theta) * this->partial_wave(helicities, s);
                }
                case helicity_frame::HELICITY_INDEPENDENT:
                {
                    return (2*_J+1) * legendre(_J, cos(_theta)) * this->partial_wave(s);
                };
                default: return NaN<complex>();
            };
        };

        // This is the relevent method that gets over-ridden, which tells us how to compute the specific J-projected amplitude
        virtual complex partial_wave(std::array<int,4> helicities, double s) = 0;
        // For helicity-independent amplitudes its useful to have this shorthand where the helicity dependence is ignored
        complex partial_wave(double s){ return partial_wave(_kinematics->helicities(0), _s); };

        // Output the J quantum number
        inline int J(){ return _J; };
        
        protected:

        // Produce a string of a parameter name which appends the J quantum number to it, i.e. "name[J]"
        inline std::string J_label(std::string name){ return name + "[" + std::to_string(_J) + "]"; };

        // Fixed spin identifier.
        // This may either be the whole-spin orbital angular momentum L
        // or half-integer total spin J
        int _J    = 0;
    };
};

#endif