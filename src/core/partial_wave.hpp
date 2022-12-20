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

#include "constants.hpp"
#include "helicities.hpp"
#include "kinematics.hpp"
#include "amplitude_options.hpp"
#include "amplitude.hpp"

#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <memory>

namespace jpacPhoto
{
    // Foward declare the PW amplitude
    class raw_partial_wave;

    // Similar to amplitude we only ever want partial_waves to be pointers
    using partial_wave = std::shared_ptr<raw_partial_wave>;

    // Summing two partial waves defaults to a "full" amplitude
    amplitude operator+(partial_wave a, partial_wave b);

    // ---------------------------------------------------------------------------
    // Methods to make partial_Waves from existing amplitudes

    // "Constructor" function which projects an existing amplitude
    amplitude project(int J, amplitude to_project, std::string id = "");

    // ---------------------------------------------------------------------------
    // Raw_amplitude class

    class raw_partial_wave : public raw_amplitude
    {
        public:

        // This constructor should be used for any user defined derived classes
        raw_partial_wave(amplitude_key key, kinematics xkinem, int J, std::string name, std::string id)
        : raw_amplitude(key, xkinem, name, id), 
          _J(J)
        {};

        // This constructor is specifically for use with the project() function
        raw_partial_wave(amplitude_key key, int J, amplitude to_project, std::string id)
        : raw_amplitude(key, to_project->_kinematics, "partial_wave", id), 
          _J(J), _amplitude(to_project)
        {
            set_N_pars(0);
        };

        // These are always assumed to be s-channel helicities so this is fixed
        helicity_channel native_helicity_frame()
        {
            return helicity_channel::S;
        };

        virtual std::vector<std::array<int,2>> allowed_meson_JP(){  return _amplitude->allowed_meson_JP(); };
        virtual std::vector<std::array<int,2>> allowed_baryon_JP(){ return _amplitude->allowed_baryon_JP(); };
        
        // Return the J-th term to the full amplitude by multiplying by angular function
        // These may be overloaded with an explicit model for the PWA
        virtual complex helicity_amplitude(std::array<int,4> helicities, double s, double t);

        // By default we calculate the partial-wave projection integral numerically from 
        // the saved amplitude 
        virtual complex evaluate(std::array<int,4> helicities, double s);

        protected:

        virtual inline void allocate_parameters(int x)
        {
            warning("partial_wave::allocate_parameters", 
                    "Partial-waves created using project() cannot change parameters! Change them in the amplitude being projected");
            return;
        };

        // By default partial_waves carry a pointer to a full amplitude and which they project
        amplitude _amplitude = nullptr;

        // Fixed spin identifier.
        // This may either be the whole-spin orbital angular momentum L
        // or half-integer total spin J
        int _J    = 0;
    };
};

#endif