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

#include "covariant/vector_exchange.hpp"
#include "analytic/vector_exchange.hpp"

void test()
{
    using namespace jpacPhoto;

    // ---------------------------------------------------------------------------
    // Couplings and constants
    // ---------------------------------------------------------------------------

    // Nucleon couplings 
    double gV_omega = 16.,    gT_omega = 0.;
    double gV_rho   = 2.4,    gT_rho   = 14.6;
    double gV_phi   = -6.2,   gT_phi   = 2.1;
    double gV_psi   = 1.6E-3, gT_psi   = 0.;
    
    // Photon couplings
    double gChi_omega   = 5.2E-4;
    double gChi_rho     = 9.2E-4;
    double gChi_phi     = 4.2E-4;
    double gChi_psi     = 1.;
    double gX_omega     = 8.2E-3;
    double gX_rho       = 3.6E-3;
    
    // Form factor cutoffs
    double LamOmega = 1.2;
    double LamRho   = 1.4; 

    // ---------------------------------------------------------------------------
    // Kinematics
    // ---------------------------------------------------------------------------

    kinematics kin  = new_kinematics(M_ZC3900);
    kin->set_meson_JP(AXIALVECTOR);

    // ---------------------------------------------------------------------------

    amplitude x = new_amplitude<covariant::vector_exchange>(kin);
    x->set_parameters({M_OMEGA, gChi_omega, gV_omega, gT_omega, LamOmega});

    amplitude y = new_amplitude<analytic::vector_exchange>(kin);
    y->set_parameters({M_OMEGA, gChi_omega, gV_omega, gT_omega, LamOmega});

    amplitude   z = cross_to<helicity_frame::S_CHANNEL>(y);

    double s = 36, t = -0.5;  
    auto hels = kin->helicities();
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