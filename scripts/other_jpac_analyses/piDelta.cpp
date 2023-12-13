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

void piDelta()
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
    amplitude b1  = new_amplitude<regge::unnatural_exchange>(kpi,  -1, "b_{1}");
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
    
    // Print differential cross section vs sqrt -t'
    // in mub!
    auto dsig = [kpi](double Egam, amplitude amp)
    {
        double s = W_cm(Egam)*W_cm(Egam);
        return [amp, kpi, s](double sqrtmt)
        {
            double tmin = kpi->t_min(s);
            double t = - sqrtmt*sqrtmt + tmin;
            return amp->differential_xsection(s, t) * 1E-3;
        };
    };

    auto dsigPol = [kpi](double Egam, int x, amplitude amp)
    {
        double s = W_cm(Egam)*W_cm(Egam);
        return [amp, kpi, x, s](double sqrtmt)
        {
            double tmin = kpi->t_min(s);
            double t = - sqrtmt*sqrtmt + tmin;
            return amp->polarized_differential_xsection(x, s, t) * 1E-3;
        };
    };

    // Beam asymmetry vs sqrt -t'
    auto BeamAsym = [kpi](double Egam, amplitude amp)
    {
        double s = W_cm(Egam)*W_cm(Egam);
        return [amp, kpi, s](double sqrtmt)
        {
            double tmin = kpi->t_min(s);
            double t = - sqrtmt*sqrtmt + tmin;
            return amp->beam_asymmetry_4pi(s, t);
        };
    };
    
    plotter plotter;
    plot p1 = plotter.new_plot();
    p1.set_logscale(false, true);
    p1.set_ranges({0., 1.4}, {1E-3, 1E2});
    p1.set_legend(0.2, 0.2);
    p1.add_header("#pi^{#minus} #Delta^{#plus#plus}");
    p1.set_labels("#sqrt{ #minus (#it{t} - #it{t}_{min})}    [GeV]", "d#sigma/d#it{t}   [#mub / GeV^{2}]");
    pimDpp(kCut);  p1.add_curve( {0, 1.4}, dsig(5,  total),  "5 GeV");
    pimDpp(kPole); p1.add_dashed({0, 1.4}, dsig(5,  total));
    pimDpp(kCut);  p1.add_curve( {0, 1.4}, dsig(8,  total),  "8 GeV");
    pimDpp(kPole); p1.add_dashed({0, 1.4}, dsig(8,  total));
    pimDpp(kCut);  p1.add_curve( {0, 1.4}, dsig(11, total), "11 GeV");
    pimDpp(kPole); p1.add_dashed({0, 1.4}, dsig(11, total));
    pimDpp(kCut);  p1.add_curve( {0, 1.4}, dsig(16, total), "16 GeV");
    pimDpp(kPole); p1.add_dashed({0, 1.4}, dsig(16,  total));
    p1.save("dxs.pdf");

    plot p2 = plotter.new_plot();
    p2.set_logscale(false, true);
    p2.set_ranges({0., 1.4}, {1E-2, 10});
    p2.set_legend(0.2, 0.2);
    p2.add_header("#pi^{#plus} #Delta^{0}");
    p2.set_labels("#sqrt{ #minus (#it{t} - #it{t}_{min})}    [GeV]", "d#sigma/d#it{t}   [#mub / GeV^{2}]");
    pipD0(kCut);  p2.add_curve( {0, 1.4}, dsig(9,  total),  "9 GeV");
    pipD0(kPole); p2.add_dashed({0, 1.4}, dsig(9,  total));

    plot p3 = plotter.new_plot();
    p3.set_ranges({0, 1}, {-1, 1});
    p3.set_legend(0.2, 0.8);
    p3.add_header("16 GeV");
    p3.set_labels("#sqrt{ #minus (#it{t} - #it{t}_{min})}    [GeV]", "#Sigma_{4#pi}");
    pipD0(kCut);   p3.add_curve( {0, 1.0}, BeamAsym(16,  total),  "#pi^{#plus} #Delta^{0}");
    pipD0(kPole);  p3.add_dashed({0, 1.0}, BeamAsym(16,  total));
    pimDpp(kCut);  p3.add_curve( {0, 1.0}, BeamAsym(16,  total),  "#pi^{#minus} #Delta^{#plus#plus}");
    pimDpp(kPole); p3.add_dashed({0, 1.0}, BeamAsym(16,  total));


    plot p4 = plotter.new_plot();
    p4.set_ranges({0, 1}, {-1, 1});
    p4.set_legend(0.2, 0.8);
    p4.add_header("9 GeV");
    p4.set_labels("#sqrt{ #minus (#it{t} - #it{t}_{min})}    [GeV]", "#Sigma_{4#pi}");
    pipD0(kCut);   p4.add_curve( {0, 1.0}, BeamAsym(9,  total),  "#pi^{#plus} #Delta^{0}");
    pipD0(kPole);  p4.add_dashed({0, 1.0}, BeamAsym(9,  total));
    pimDpp(kCut);  p4.add_curve( {0, 1.0}, BeamAsym(9,  total),  "#pi^{#minus} #Delta^{#plus#plus}");
    pimDpp(kPole); p4.add_dashed({0, 1.0}, BeamAsym(9,  total));

    plot p5 = plotter.new_plot();
    p5.set_logscale(false, true);
    p5.set_legend(0.2, 0.8);
    p5.add_header("16 GeV");
    p5.set_ranges({0, 1}, {4E-3, 10});
    p5.set_labels("#sqrt{ #minus (#it{t} - #it{t}_{min})}    [GeV]", "d#sigma/d#it{t}   [#mub / GeV^{2}]");
    pipD0(kCut);   p5.add_curve( {0, 1}, dsigPol(16, +1, total), "#pi^{#plus} #Delta^{0}");
    pipD0(kPole);  p5.add_dashed({0, 1}, dsigPol(16, +1, total));
    pimDpp(kCut);  p5.add_curve( {0, 1}, dsigPol(16, +1, total), "#pi^{#minus} #Delta^{#plus#plus}");
    pimDpp(kPole); p5.add_dashed({0, 1}, dsigPol(16, +1, total));

    plot p6 = plotter.new_plot();
    p6.set_logscale(false, true);
    p6.set_legend(0.2, 0.8);
    p6.add_header("16 GeV");
    p6.set_ranges({0, 1}, {4E-3, 10});
    p6.set_labels("#sqrt{ #minus (#it{t} - #it{t}_{min})}    [GeV]", "d#sigma/d#it{t}   [#mub / GeV^{2}]");
    pipD0(kCut);   p6.add_curve( {0, 1}, dsigPol(16, -1, total), "#pi^{#plus} #Delta^{0}");
    pipD0(kPole);  p6.add_dashed({0, 1}, dsigPol(16, -1, total));
    pimDpp(kCut);  p6.add_curve( {0, 1}, dsigPol(16, -1, total), "#pi^{#minus} #Delta^{#plus#plus}");
    pimDpp(kPole); p6.add_dashed({0, 1}, dsigPol(16, -1, total));

    plotter.combine({3, 2}, {p1, p6, p3, p2, p5, p4}, "piDelta.pdf");
    return;
};