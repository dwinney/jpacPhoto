// Template file to give an idea how to specify a AmpTools interface struct
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------

#ifndef TEMPLATE_HPP
#define TEMPALTE_HPP

#include "constants.hpp"
#include "amplitude.hpp"
#include "kinematics.hpp"

#include "some_jpacPhoto_amplitude.hpp"
#include "some_other_amplitude.hpp"

namespace jpacPhoto
{
    struct my_process
    {
        // Here we should set up the kinematics and amplitude we want to inject into AmpTools
        // parameters accessable by jpacPhoto::amplitude::set_parameters() will be set by 
        // AmpTools' config file interface
        inline static one_meson::amplitude initialize_amplitude()
        {
            using namespace jpacPhoto::one_meson;

            kinematics kin = new_kinematics( M_MESON, M_BARYON );
            kin->set_meson_JP( 1, +1);
            kin->set_baryon_JP(3, +1);

            amplitude amp1 = new_amplitude<some_jpacPhoto_amplitude>(kin, "first amp" );
            amplitude amp2 = new_amplitude<some_other_amplitude>(    kin, "second amp");
            return amp1 + amp2;
        };

       
        // Here we specify what we will use as the intensity
        // In the simple example of no subsequent decays and no polarization information
        // we can use the differential cross section directly
        // In more complicated processes, here we calculate whatever we want
        inline static double intensity( double s, double t, one_meson::amplitude amp) 
        { 
            return amp->differential_xsection(s, t);               
        };
    };
};

#endif