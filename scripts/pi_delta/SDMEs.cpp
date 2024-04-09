// Model for high energy photoproduction of the pi Delta system via natural
// and unnatural exchanges.
//
// Adapted from models and codes originally written by J. Nys and V. Mathieu in [1]
// Plots model versus new data from GlueX
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Univeritat Bonn (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// References:
// [1] - https://arxiv.org/abs/1710.09394
// ------------------------------------------------------------------------------


#include "kinematics.hpp"
#include "amplitude.hpp"
#include "regge/natural_exchange.hpp"
#include "regge/unnatural_exchange.hpp"
#include "piDelta/plots.hpp"
#include "piDelta/data.hpp"

#include "plotter.hpp"

void SDMEs()
{
    using namespace jpacPhoto;

    //--------------------------------------------------------------------------
    // Couplings and things

    int kPole = regge::natural_exchange::kPole;
    int kCut  = regge::natural_exchange::kCut;

    //--------------------------------------------
    // Pole model 

    // Degenerate natural trajectory
    double aPN = 0.9, a0N  = 1 - aPN*M_RHO*M_RHO; 
    double aPU = 0.7, a0U  =   - aPU*M_PION*M_PION; 
    
    // Top couplings
    double ga2 = 0.204268,  grho = 0.109677;
    double gpi = 0.428507,  gb1  = 0.0976007;

    // Bottom couplings
    double g0 = 19.16;                        // gpiDelta
    double g1 = -48.2, g2 = -52.4, g3 = 40.2; // grhoDelta

    // Suppression factors
    double bN = -0.4161, bU = 0.06043;

    // Pion absorption constant
    double cpi = 1.0587;
    
    //-----------------------------------------
    // Cut model 

    double bUc =  0.135,    bNc =  2.122,  cpic =  1.041;
    double g1c = -370.777, g2c = -242.416, g3c  = -138.9583;
    
    //--------------------------------------------------------------------------
    // Kinematics for pi Delta

    kinematics kpi = new_kinematics(M_PION, M_DELTA);
    kpi->set_meson_JP( PSEUDOSCALAR );
    kpi->set_baryon_JP( THREEPLUS );

    //--------------------------------------------------------------------------
    // Natural exchanges

    amplitude rho = new_amplitude<regge::natural_exchange>(kpi, -1, "#rho");
    rho->set_parameters({a0N, aPN,  grho, g1, g2, g3, bN});
    
    amplitude a2  = new_amplitude<regge::natural_exchange>(kpi, +1, "a_{2}");
    a2->set_parameters( {a0N, aPN, -ga2,  g1, g2, g3, bN});

    amplitude N = a2 + rho;
    N->set_id("Natural");

    //--------------------------------------------------------------------------
    // Unnatural exchanges

    amplitude pi  = new_amplitude<regge::unnatural_exchange>(kpi, +1, "#pi");
    pi->set_parameters( {a0U, aPU, -gpi,  g0, bU, cpi});

    amplitude b1  = new_amplitude<regge::unnatural_exchange>(kpi, -1, "b_{1}");
    b1->set_parameters( {a0U, aPU,  gb1,  g0, bU, cpi});

    amplitude U = pi + b1;
    U->set_id("Unnatural");

    //--------------------------------------------------------------------------
    // Total sum

    amplitude total = N + U;
    total->set_id("Sum");

    plotter plotter;

    plot pdif = piDelta::plot_differential(plotter);
    pdif.add_curve({1E-3, 1.}, [&](double mt){ return total->differential_xsection(s_cm(8.0), -mt) * 1E-3; }, "Pole model");

    plot pBA = piDelta::plot_beam_asymmetry(plotter);
    pBA.add_curve({1E-3, 1.}, [&](double mt){ return total->beam_asymmetry_4pi(s_cm(8.0), -mt); }, "Pole model");

    std::vector<plot> pSDMEs = piDelta::plot_SDMEs(plotter);
    for (int i = 1; i <= 9; i++)
    {
        std::array<int,3> ind = piDelta::SDME_indices(i);
        if (ind[0] == 2) pSDMEs[i-1].add_curve({1E-3, 1.}, [&](double mt){ return std::imag(total->bSDME_H(ind[0], ind[1], ind[2], s_cm(8.5), -mt)); }, "Pole model");
        else             pSDMEs[i-1].add_curve({1E-3, 1.}, [&](double mt){ return std::real(total->bSDME_H(ind[0], ind[1], ind[2], s_cm(8.5), -mt)); }, "Pole model");
    };

    plotter.combine({2,1}, {pdif, pBA}, "diff_and_BA.pdf");
    plotter.combine({3,3},  pSDMEs,     "SDMEs.pdf");
};