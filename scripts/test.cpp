// Photoproduction cross-sections for axial vector Z-mesons in charm and bottom sectors
// Reproduces fig. 2 of [1] in pdf form.
//
// OUTPUT: Z_mesons.pdf
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------
// REFERENCES:
// [1] arXiv:2008.01001 [hep-ph]
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "kinematics.hpp"
#include "plotter.hpp"
#include "crossing.hpp"

#include "analytic/pseudoscalar_exchange.hpp"
#include "covariant/pseudoscalar_exchange.hpp"
#include "regge/pseudoscalar_exchange.hpp"

void test()
{
    using namespace jpacPhoto;

    // ---------------------------------------------------------------------------
    // Couplings and constants
    // ---------------------------------------------------------------------------

    // Bottom vertex coupling (pi - nucleon - nucleon)
    double g_piNN = sqrt(2) * sqrt(4*PI*13.81); 

    // Cutoff for exponential form factor
    double lambda_pi = .900;  // MeV 

    // Zc(3900) couplings 
    double gc_jpsi  = 1.91; // psi coupling before VMD scaling
    double gc_gamma = E * F_JPSI * gc_jpsi / M_JPSI;

    // ---------------------------------------------------------------------------
    // Kinematics
    // ---------------------------------------------------------------------------

    kinematics kZc  = new_kinematics(M_ZC3900);
    kZc->set_meson_JP(AXIALVECTOR);

    // ---------------------------------------------------------------------------

    // Covariant amplitudes always defined in the s-channel
    amplitude  x  = new_amplitude<covariant::pseudoscalar_exchange>(kin);
    // Analytic expressions usually defined in the t-channel
    amplitude  y  = new_amplitude<analytic::pseudoscalar_exchange>(kin);

    // We want to sum these so we can do:
    amplitude z = x + y; // but this will produce an error because helicities in different channels
    amplitude z = x + cross_to<helicity_frame::S_CHANNEL>(y); // Both defined in the s-channel
    amplitude z = cross_to<helicity_frame::T_CHANNEL>(x) + y; // Both defined in the t-channel



    double s = 36, t = -0.5;  
    auto hels = kZc->helicities();
    auto xamps = x->get_cache(s,t);
    auto zamps = z->get_cache(s,t);
    for (int i = 0; i < 12; i++)
    {
        print(print_helicities(hels[i]), xamps[i], zamps[i]);
    };

    print(x->differential_xsection(s,t));
    print(y->differential_xsection(s,t));
    print(z->differential_xsection(s,t));
};