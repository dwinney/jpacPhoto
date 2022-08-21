#include "triple_regge.hpp"
#include "inclusive_kinematics.hpp"
#include "regge_trajectory.hpp"

#include "jpacUtils.hpp"
#include "jpacGraph1D.hpp"

#include <iostream>
#include <cstring>

#include "Math/GSLIntegrator.h"
#include "Math/IntegrationTypes.h"
#include "Math/Functor.h"

using namespace jpacPhoto;

void deltapp()
{
    // ---------------------------------------------------------------------------
    // Amplitudes
    // ---------------------------------------------------------------------------

    // ---------------------------------------------------------------------------
    // Global quantities

    // Couplings
    double g_b1 = 0.24;
    double g_N = sqrt(2.) * sqrt(4. * PI * 13.81); // Nucleon coupling same for all

    // Delta p pi coupling
    double g_delta = 18.5;
    double LamPi = .9;  // 900 MeV cutoff for formfactor
    
    // ---------------------------------------------------------------------------
    // Fixed-spin pion amplitudes

    // Kinematics for b1 Δ++ final state
    reaction_kinematics kDelta (M_B1, M_DELTA);
    kDelta.set_meson_JP( 1, +1); 
    kDelta.set_baryon_JP(3, +1); 

    // Exclusive amplitude for delta
    pseudoscalar_exchange excDelta (&kDelta, M_PION, "b_{1}^{-} #Delta^{++}");
    excDelta.set_params({g_b1, g_delta});
    excDelta.set_formfactor(1, LamPi);

    // Kinematics b1 proton
    reaction_kinematics kN (M_B1, M_PROTON);
    kN.set_meson_JP( 1, +1);
    kN.set_baryon_JP(1, +1); 

    // Exclusive amplitude for nucleon
    pseudoscalar_exchange excN (&kN, M_PION, "b1 production");
    excN.set_params({g_b1, g_N});
    excN.set_formfactor(1, LamPi);

    // We now can pass this to an inclusive amplitude
    triple_regge incB1 (&excN);
    incB1.set_high_energy_approximation(false);
    
    // // ---------------------------------------------------------------------------
    // // Plotting options
    // // --------------------------------------------------------------------------

    int N = 100;

    double xmin = 2.2;   
    double xmax = 4.;

    double ymin = 0.;
    double ymax = 10.;

    std::string filename = "integrated.pdf";
    std::string ylabel   = "#sigma  [#mub]";
    std::string xlabel   = "#it{W}_{#gamma#it{p}}  [GeV]";

    // ---------------------------------------------------------------------------  
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    jpacGraph1D * plotter = new jpacGraph1D();

    // Stable delta
    auto F = [&](double w)
    {
        return excDelta.integrated_xsection(w*w) * 1.E-3; // in mub!
    };
    
    // Sill line-shape
    // Mass and width from 2106.03749
    auto Sill = [&](double w)
    {
        double width = 90.4E-3;
        double mass  = 1236.2E-3;
        double mN = 0.938, mpi = 0.134;
        double sth = (mN + mpi)*(mN+mpi);

        double gamma = width*mass / sqrt(mass*mass - sth);

        return 2.*w/M_PI * sqrt(w*w - sth) *gamma / (pow(w*w-mass*mass, 2.) + pow(sqrt(w*w - sth)*gamma, 2.));
    };

    auto H = [&](double w)
    {
        auto dH = [&](double m)
        {
            
            kDelta.set_recoil_mass(m);
            return Sill(m)*excDelta.integrated_xsection(w*w) * 1.E-3; // in mub!
        };

        ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS15);
        ROOT::Math::Functor1D wH(dH);
        ig.SetFunction(wH);

        return ig.Integral(M_PION+M_PROTON, 3.);
    };

    // Inclusive
    auto G = [&](double w)
    {
        return incB1.integrated_xsection(w*w); // in mub!
    };  

    plotter->AddEntry(N, H, {xmin, xmax},   "#it{b}_{1}^{#minus} (#Delta^{#plus#plus} #rightarrow #pi^{#plus} #it{p}) from BW", 1);
    kDelta.set_recoil_mass(M_DELTA);
    plotter->AddDashedEntry(N, F, {xmin, xmax});

    incB1.set_sigma_total(JPAC_pipp_onlyDelta);
    plotter->AddEntry(N, G, {xmin, xmax},   "#it{b}_{1}^{#minus} (#Delta^{#plus#plus} #rightarrow #pi^{#plus} #it{p}) from SAID", 1);

    incB1.set_sigma_total(JPAC_pipp_withResonances);
    plotter->AddEntry(N, G, {xmin, xmax},   "Inclusive #it{b}_{1}^{#minus}", 1);

    // ---------------------------------------------------------------------------
    // Finally make the plot pretty

    // Axes options
    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);

    // LEgend options
    plotter->SetLegend(0.5, 0.65);
    plotter->SetLegendOffset(0.3, 0.13);

    // Output to file
    plotter->Plot(filename);

    delete plotter;
};