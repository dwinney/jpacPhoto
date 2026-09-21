// Photoproduction cross-sections for X(3872) and axial-vector charmonium mesons
// Reproduces fig. 3 of [1] in pdf form.
//
// OUTPUT: X_mesons.pdf
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
#include "crossing.hpp"
#include "plotter.hpp"

#include "covariant/vector_exchange.hpp"
#include "analytic/vector_exchange.hpp"
#include "regge/vector_exchange.hpp"

void test()
{
    using namespace jpacPhoto;

    // ---------------------------------------------------------------------------
    // Couplings and constants
    // ---------------------------------------------------------------------------'

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
    
    // Form factor cutoffs
    double LamOmega = 1.2;
    double LamRho   = 1.4; 

    // ---------------------------------------------------------------------------
    // Kinematics
    // ---------------------------------------------------------------------------
    kinematics kChiC1  = new_kinematics( M_CHIC1 );
    kChiC1->set_meson_JP(AXIALVECTOR);

    // ---------------------------------------------------------------------------
    // Low Energy amplitudes (Fixed-spin vector exchange)
    // ---------------------------------------------------------------------------

    // chi_c1
    amplitude x = new_amplitude<covariant::vector_exchange>(kChiC1);
    x->set_parameters({M_RHO, gChi_rho, gV_rho, gT_rho, LamRho});

    amplitude y = new_amplitude<analytic::vector_exchange>(kChiC1);
    y->set_parameters({M_RHO, gChi_rho, gV_rho, gT_rho, LamRho});

    amplitude z = cross_to<helicity_frame::S_CHANNEL>(y);

    double s = 36., t = -0.3;

    auto xamps = x->get_cache(s, t);
    auto zamps = z->get_cache(s, t);

    for (int i = 0; i < 12; i++)
    {
        print(i, xamps[i], zamps[i]);
    };

    line();
    print(x->differential_xsection(s,t));
    print(y->differential_xsection(s,t));
    print(z->differential_xsection(s,t));
};