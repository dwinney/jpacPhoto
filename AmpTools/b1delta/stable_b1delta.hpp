// Production of a stable B1 Delta final state via pion and rho meson exchanges
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------

#ifndef B1Delta_FINALSTATE_HPP
#define B1Delta_FINALSTATE_HPP

#include "constants.hpp"
#include "regge/unnatural_exchange.hpp"
#include "regge/natural_exchange.hpp"
#include <string>

namespace jpacPhoto
{
    struct stable_b1delta
    {
        inline static one_meson::amplitude initialize_amplitude()
        {
            using namespace jpacPhoto::one_meson;

            kinematics kB1Delta = new_kinematics( M_B1, M_DELTA );
            kB1Delta->set_meson_JP( 1, +1);
            kB1Delta->set_baryon_JP(3, +1);

            amplitude pi  = new_amplitude<regge::unnatural_exchange>(kB1Delta, "pi  exchange");
            amplitude rho = new_amplitude<regge::natural_exchange>(  kB1Delta, "rho exchange");
            return pi + rho;
        };

        // For a 2->2 process we only want the mandelstam invariants s and t
        enum UsersVars { kVar_s = 0, kVar_t = 1 };
        inline static unsigned int N_variables() { return 2; }
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

        // Use the differential cross-section as unpolarized intensity
        inline static double intensity( GDouble* userVars, jpacPhoto::one_meson::amplitude amp) 
        { 
            return amp->differential_xsection(userVars[kVar_s], userVars[kVar_t]);               
        };
    };
};

#endif