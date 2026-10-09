// Theres a few classes which I want to prohibit being created on the stack
// We use this struct to lock out the public contructors
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef KEY_HPP
#define KEY_HPP

#include <memory>

namespace jpacPhoto
{
    class raw_kinematics;
    class raw_amplitude;
    class raw_partial_wave;
    class raw_semi_inclusive;
    enum helicity_frame : uint;

    class key
    {
        private:

        // Private constructor only accessable via friend methods below
        key(){};

        // Kinematics factories
        friend std::shared_ptr<raw_kinematics> new_kinematics(double, double, double, double);
        friend std::shared_ptr<raw_kinematics> new_kinematics(double, double);

        // Amplitude factories        
        template<class A>             friend std::shared_ptr<raw_amplitude> new_amplitude(std::shared_ptr<raw_kinematics>);
        template<class A, typename B> friend std::shared_ptr<raw_amplitude> new_amplitude(std::shared_ptr<raw_kinematics>, B);
        friend std::shared_ptr<raw_amplitude> operator+(std::shared_ptr<raw_amplitude>, std::shared_ptr<raw_amplitude>);
        template<helicity_frame C> friend std::shared_ptr<raw_amplitude> cross_to(std::shared_ptr<raw_amplitude>);

        // Partial wave factories
        template<class A>             friend std::shared_ptr<raw_partial_wave> new_partial_wave(std::shared_ptr<raw_kinematics>, uint);
        template<class A, typename B> friend std::shared_ptr<raw_partial_wave> new_partial_wave(std::shared_ptr<raw_kinematics>, uint, B);
        friend std::shared_ptr<raw_partial_wave> project(uint, std::shared_ptr<raw_amplitude>);

        // Inclusive factories     
        template<class A>          friend std::shared_ptr<raw_semi_inclusive> new_semi_inclusive(std::shared_ptr<raw_kinematics> kin);
        template<class A, class B> friend std::shared_ptr<raw_semi_inclusive> new_semi_inclusive(std::shared_ptr<raw_kinematics> kin, B extra);
        friend std::shared_ptr<raw_semi_inclusive> operator+(std::shared_ptr<raw_semi_inclusive> a, std::shared_ptr<raw_semi_inclusive> b);
        friend std::shared_ptr<raw_semi_inclusive> operator+(std::shared_ptr<raw_semi_inclusive> a, std::shared_ptr<raw_amplitude> b);
    };
};

#endif