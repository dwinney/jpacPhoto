// ---------------------------------------------------------------------------
// Prediction for Z_c and Z_b photoproduction based on pi exchange
// at low energies. 
//
// Reproduces left plot in FIG 3 of [1] 
// 
// USAGE:
// make Z_low && ./Z_low
//
// OUTPUT:
// Z_FS.pdf
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// REFERENCES:
// [1] arXiv:2008.01001 [hep-ph]
// ---------------------------------------------------------------------------


#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "amplitudes/pseudoscalar_exchange.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <cstring>

using namespace jpacPhoto;

int main( int argc, char** argv )
{

  // ---------------------------------------------------------------------------
  // Preliminaries
  // ---------------------------------------------------------------------------

  double g_NN = sqrt(2.) * sqrt(4. * M_PI * 13.81); // Nucleon coupling same for all
  double LamPi = .9;  // 900 MeV cutoff for formfactor

  // Zc(3900)
  double mZc = 3.8884; // GeV
  reaction_kinematics * kZc = new reaction_kinematics(mZc);
  kZc->set_meson_JP(AXIAL_VECTOR);

  double gc_Psi = 1.91; // psi coupling before VMD scaling
  double gc_Gamma = E * F_JPSI * gc_Psi / M_JPSI;
  std::vector<double> Zc_couplings = {gc_Gamma, g_NN};

  // Zb(10610)
  double mZb = 10.6072;
  reaction_kinematics * kZb = new reaction_kinematics(mZb);
  kZb->set_meson_JP(AXIAL_VECTOR);

  double gb_Ups1 = 0.49, gb_Ups2 = 3.30, gb_Ups3 = 9.22;
  double gb_Gamma = E * (F_UPSILON1S * gb_Ups1 / M_UPSILON1S 
                       + F_UPSILON2S * gb_Ups2 / M_UPSILON2S
                       + F_UPSILON3S * gb_Ups3 / M_UPSILON3S);  
  std::vector<double> Zb_couplings = {gb_Gamma, g_NN};

  
  // Zb(10650)
  double mZbp = 10.6522;
  reaction_kinematics * kZbp = new reaction_kinematics(mZbp);
  kZbp->set_meson_JP(AXIAL_VECTOR);

  double gbp_Ups1 = 0.21, gbp_Ups2 = 1.47, gbp_Ups3 = 4.8;
  double gbp_Gamma = E * (F_UPSILON1S * gbp_Ups1 / M_UPSILON1S 
                       +  F_UPSILON2S * gbp_Ups2 / M_UPSILON2S
                       +  F_UPSILON3S * gbp_Ups3 / M_UPSILON3S);  
  std::vector<double> Zbp_couplings = {gbp_Gamma, g_NN};
  
  // ---------------------------------------------------------------------------
  // Fixed-spin amplitudes
  // ---------------------------------------------------------------------------

  pseudoscalar_exchange Zc_fixedspin(kZc, M_PION, "#it{Z_{c}} (3900)^{+}");
  Zc_fixedspin.set_params(Zc_couplings);
  Zc_fixedspin.set_formfactor(1, LamPi);

  pseudoscalar_exchange Zb_fixedspin(kZb, M_PION,  "#it{Z_{b}} (10610)^{+}");
  Zb_fixedspin.set_params(Zb_couplings);
  Zb_fixedspin.set_formfactor(1, LamPi);

  pseudoscalar_exchange Zbp_fixedspin(kZbp, M_PION, "#it{Z'_{b}} (10650)^{+}");
  Zbp_fixedspin.set_params(Zbp_couplings);
  Zbp_fixedspin.set_formfactor(1, LamPi);

  // ---------------------------------------------------------------------------
  // Plotting options
  // ---------------------------------------------------------------------------
  
  // which amps to plot
  std::vector<amplitude*> amps;
  amps.push_back(&Zc_fixedspin);
  amps.push_back(&Zb_fixedspin);
  amps.push_back(&Zbp_fixedspin);

  int   N = 200;
  bool PRINT_TO_COMMANDLINE     = true;

  // ---------------------------------------------------------------------------
  double  xmin = 4.;
  double  xmax = 20.;

  double  ymin = 2.E-2;
  double  ymax = 100.;
  // ---------------------------------------------------------------------------

  std::string ylabel    = "#it{#sigma(#gamma p #rightarrow Z n)}  [nb]";
  std::string xlabel    = "#it{W_{#gammap}}  [GeV]";
  std::string filename  = "Z_FS.pdf";

  // ---------------------------------------------------------------------------
  // You shouldnt need to change anything below this line
  // ---------------------------------------------------------------------------

  // Plotter object
  jpacGraph1D* plotter = new jpacGraph1D();

  // ---------------------------------------------------------------------------
  // Print the desired observable for each amplitude
  for (int n = 0; n < amps.size(); n++)
  {
    auto F = [&](double x)
    {
      return amps[n]->integrated_xsection(x*x);
    };

    plotter->AddEntry(N, F, {xmin,xmax}, amps[n]->get_id(), PRINT_TO_COMMANDLINE);
  }
  
  plotter->SetXaxis(xlabel, xmin, xmax);
  plotter->SetYaxis(ylabel, ymin, ymax);
  plotter->SetYlogscale(1);
  plotter->SetLegend(0.7, 0.65);

  // Output to file
  plotter->Plot(filename);

  return 1.;
}
