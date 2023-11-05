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
    struct my_decay_process
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

        // Here we need to specify what variables we need to calculate the intensity
        // We always assume you dont need the full four-vectors to calculate the final
        // intensity since we have an analytic formula

        // For a 2->2 (unpolarized) process we need only 2 variables: s and t
        // We can generalize this to have 3 and include a polarization angle
        enum UsersVars { kVar_s = 0, kVar_t = 1 };
        inline static unsigned int N_variables() { return 2; }

        // Take in the 4-vectors and calculate the variables we need
        inline static void  calculate_variables( GDouble** pKin, GDouble* userVars )
        {
            TLorentzVector meson( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0]);
            TLorentzVector baryon(pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0]);

            double s = (meson + baryon).M2();
            userVars[kVar_s] = s;

            double Egam = jpacPhoto::E_beam(sqrt(s));
            TLorentzVector beam( 0, 0, Egam, Egam);
            double t = (beam - meson).M2();
            userVars[kVar_t] = t;
        };

        // Here we specify what we will use as the intensity
        // In the simple example of no subsequent decays we can use the differential
        // cross section directly.
        // In more complicated processes, here we calculate the decay chains
        inline static double intensity( GDouble* userVars, one_meson::amplitude amp) 
        { 
            return amp->differential_xsection(userVars[kVar_s], userVars[kVar_t]);               
        };
    };
};

#endif