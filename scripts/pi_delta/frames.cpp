// Model for high energy photoproduction of the pi Delta system via natural
// and unnatural exchanges.
//
// Adapted from models and codes originally written by J. Nys and V. Mathieu 
// Reproduces fig. 2 in [1].
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
#include "plotter.hpp"

void frames()
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
    amplitude a2  = new_amplitude<regge::natural_exchange>(kpi, +1, "a_{2}");
    amplitude N = a2 + rho;
    N->set_id("Natural");

    //--------------------------------------------------------------------------
    // Unnatural exchanges

    amplitude pi  = new_amplitude<regge::unnatural_exchange>(kpi, +1, "#pi");
    amplitude b1  = new_amplitude<regge::unnatural_exchange>(kpi, -1, "b_{1}");
    amplitude U = pi + b1;
    U->set_id("Unnatural");

    //--------------------------------------------------------------------------
    // Total sum

    amplitude total = N + U;
    total->set_id("Sum");

    //--------------------------------------------------------------------------
    // Isospin states
    
    // pi- Delta++
    auto pimDpp = [&](int x)
    {
        rho->set_option(x); a2->set_option(x);
        if (x == regge::natural_exchange::kPole)
        {
            rho->set_parameters({a0N, aPN,  grho, g1, g2, g3, bN});
            a2->set_parameters( {a0N, aPN, -ga2,  g1, g2, g3, bN});
            pi->set_parameters( {a0U, aPU, -gpi,  g0, bU, cpi});
            b1->set_parameters( {a0U, aPU,  gb1,  g0, bU, cpi});
        }
        else 
        {
            rho->set_parameters({a0N, aPN,  grho, g1c, g2c, g3c, bNc});
            a2->set_parameters( {a0N, aPN, -ga2,  g1c, g2c, g3c, bNc});
            pi->set_parameters( {a0U, aPU, -gpi,  g0, bUc, cpic});
            b1->set_parameters( {a0U, aPU,  gb1,  g0, bUc, cpic});
        }
    };

    // pi+ Delta0
    auto pipD0 = [&](int x)
    {
        rho->set_option(x); a2->set_option(x);
        if (x == regge::natural_exchange::kPole)
        {
            rho->set_parameters({a0N, aPN,  grho, g1, g2, g3, bN});
            a2->set_parameters( {a0N, aPN,  ga2 , g1, g2, g3, bN});
            pi->set_parameters( {a0U, aPU,  gpi , g0, bU, cpi});
            b1->set_parameters( {a0U, aPU,  gb1 , g0, bU, cpi});
        }
        else 
        {
            rho->set_parameters({a0N, aPN,  grho, g1c, g2c, g3c, bNc});
            a2->set_parameters( {a0N, aPN,  ga2 , g1c, g2c, g3c, bNc});
            pi->set_parameters( {a0U, aPU,  gpi , g0,  bUc, cpic});
            b1->set_parameters( {a0U, aPU,  gb1 , g0,  bUc, cpic});
        }
    };

    //--------------------------------------------------------------------------
    // Plots
    
    // Print  SDMEs section vs sqrt -t
    auto SDME_GJ = [](int alpha, int m, int mp, double Egam, amplitude amp)
    {
        double s = W_cm(Egam)*W_cm(Egam);
        return [amp, alpha, m, mp, s](double mt)
        { 
            double x = (alpha == 2) ? std::imag(amp->bSDME_GJ(alpha, m, mp, s, -mt)) : std::real(amp->bSDME_GJ(alpha, m, mp, s, -mt)); 
            return /* amp->differential_xsection(s, -mt) * */ x /* * 1E-3 */; 
        };
    };
    
    auto SDME_H = [](int alpha, int m, int mp, double Egam, amplitude amp)
    {
        double s = W_cm(Egam)*W_cm(Egam);
        return [amp, alpha, m, mp, s](double mt)
        { 
            double x = (alpha == 2) ? std::imag(amp->bSDME_H(alpha, m, mp, s, -mt)) : std::real(amp->bSDME_H(alpha, m, mp, s, -mt)); 
            return /* amp->differential_xsection(s, -mt) * */ x /* * 1E-3 */; 
        };
    };
    

    std::array<double,2> ybounds = {-0.5, 0.5};

    plotter plotter;
    plot p1 = plotter.new_plot();
    p1.set_ranges({0,1}, ybounds);
    p1.set_legend(0.2, 0.2);
    p1.add_header("#pi^{#minus} #Delta^{#plus#plus}");
    p1.set_labels("-t  [GeV^{2}]", "#rho_{11}^{0}");
    pimDpp(kPole);  p1.add_curve(  {0.01, 1.0}, SDME_H(0, 1, 1, 8, total),  "H");
    pimDpp(kPole);  p1.add_curve(  {0.01, 1.0}, SDME_GJ(0, 1, 1, 8, total), "GJ");

    plot p2 = plotter.new_plot();
    p2.set_legend(0.2, 0.2);
    p2.set_ranges({0,1},ybounds);
    p2.add_header("#pi^{#minus} #Delta^{#plus#plus}");
    p2.set_labels("-t  [GeV^{2}]", "Re #rho_{31}^{0}");
    pimDpp(kPole);  p2.add_curve(  {0.01, 1.0}, SDME_H(0, 3, 1, 8, total),  "H");
    pimDpp(kPole);  p2.add_curve(  {0.01, 1.0}, SDME_GJ(0, 3, 1, 8, total), "GJ");

    plot p3 = plotter.new_plot();
    p3.set_legend(0.2, 0.2);
    p3.set_ranges({0,1}, ybounds);
    p3.add_header("#pi^{#minus} #Delta^{#plus#plus}");
    p3.set_labels("-t  [GeV^{2}]", "Re #rho_{3-1}^{0}");
    pimDpp(kPole);  p3.add_curve(  {0.01, 1.0}, SDME_H(0, 3, -1, 8, total),  "H");
    pimDpp(kPole);  p3.add_curve(  {0.01, 1.0}, SDME_GJ(0, 3, -1, 8, total), "GJ");

    plot p4 = plotter.new_plot();
    p4.set_legend(0.2, 0.2);
    p4.set_ranges({0,1}, ybounds);
    p4.add_header("#pi^{#minus} #Delta^{#plus#plus}");
    p4.set_labels("-t  [GeV^{2}]", "#rho_{11}^{1}");
    pimDpp(kPole);  p4.add_curve(  {0.01, 1.0}, SDME_H(1, 1, 1, 8, total),  "H");
    pimDpp(kPole);  p4.add_curve(  {0.01, 1.0}, SDME_GJ(1, 1, 1, 8, total), "GJ");

    plot p5 = plotter.new_plot();
    p5.set_legend(0.2, 0.2);
    p5.set_ranges({0,1}, ybounds);
    p5.add_header("#pi^{#minus} #Delta^{#plus#plus}");
    p5.set_labels("-t  [GeV^{2}]", "#rho_{33}^{1}");
    pimDpp(kPole);  p5.add_curve(  {0.01, 1.0}, SDME_H(1, 3, 3, 8, total),  "H");
    pimDpp(kPole);  p5.add_curve(  {0.01, 1.0}, SDME_GJ(1, 3, 3, 8, total), "GJ");

    plot p6 = plotter.new_plot();
    p6.set_legend(0.2, 0.2);
    p6.set_ranges({0,1}, ybounds);
    p6.add_header("#pi^{#minus} #Delta^{#plus#plus}");
    p6.set_labels("-t  [GeV^{2}]", "Re #rho_{31}^{1}");
    pimDpp(kPole);  p6.add_curve(  {0.01, 1.0}, SDME_H(1, 3, 1, 8, total),  "H");
    pimDpp(kPole);  p6.add_curve(  {0.01, 1.0}, SDME_GJ(1, 3, 1, 8, total), "GJ");

    plot p7 = plotter.new_plot();
    p7.set_legend(0.2, 0.2);
    p7.set_ranges({0,1}, ybounds);
    p7.add_header("#pi^{#minus} #Delta^{#plus#plus}");
    p7.set_labels("-t  [GeV^{2}]", "Re #rho_{3-1}^{1}");
    pimDpp(kPole);  p7.add_curve(  {0.01, 1.0}, SDME_H(1, 3, -1, 8, total),  "H");
    pimDpp(kPole);  p7.add_curve(  {0.01, 1.0}, SDME_GJ(1, 3, -1, 8, total), "GJ");

    plot p8 = plotter.new_plot();
    p8.set_legend(0.2, 0.2);
    p8.set_ranges({0,1}, ybounds);
    p8.add_header("#pi^{#minus} #Delta^{#plus#plus}");
    p8.set_labels("-t  [GeV^{2}]", "Im #rho_{31}^{2}");
    pimDpp(kPole);  p8.add_curve(  {0.01, 1.0}, SDME_H(2, 3, 1, 8, total),  "H");
    pimDpp(kPole);  p8.add_curve(  {0.01, 1.0}, SDME_GJ(2, 3, 1, 8, total), "GJ");

    plot p9 = plotter.new_plot();
    p9.set_legend(0.2, 0.2);
    p9.set_ranges({0,1}, ybounds);
    p9.add_header("#pi^{#minus} #Delta^{#plus#plus}");
    p9.set_labels("-t  [GeV^{2}]", "Im #rho_{3-1}^{2}");
    pimDpp(kPole);  p9.add_curve(  {0.01, 1.0}, SDME_H(2, 3, -1, 8, total),  "H");
    pimDpp(kPole);  p9.add_curve(  {0.01, 1.0}, SDME_GJ(2, 3, -1, 8, total), "GJ");
    
    plotter.combine({3,3}, {p1, p2, p3, p4, p5, p6, p7, p8, p9}, "sdmes_old.pdf");

    return;
};