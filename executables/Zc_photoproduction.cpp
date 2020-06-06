// ---------------------------------------------------------------------------
// Photoproduction of Zc by a  charged pion exchange
// Reproduces the results from [1]
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// References:
// [1] arXiv:1503.02125 [hep-ph]
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "amplitudes/pseudoscalar_exchange.hpp"

#include "jpacGraph1D.hpp"

#include <cstring>

int main( int argc, char** argv )
{
  double theta = 0.;
  double max = 25;
  double y[2] = {0., 0.1};
  int N = 100;
  std::string filename = "pseudoscalar_exchange.pdf";
  bool integ = true;
  for (int i = 0; i < argc; i++)
  {
    if (std::strcmp(argv[i],"-f")==0) filename = argv[i+1];
    if (std::strcmp(argv[i],"-c")==0) theta = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-m")==0) max = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-n")==0) N = atoi(argv[i+1]);
    if (std::strcmp(argv[i],"-y")==0) y_range(argv[i+1], y);
    if (std::strcmp(argv[i],"-diff")==0) integ = false;
  }

  // Set up kinematics for the chi_c1
  reaction_kinematics * ptr = new reaction_kinematics(4.20, "Z_{c}^{+}(4200)");

  pseudoscalar_exchange amp(ptr, mPi, "#pi meson exchange");
  amp.set_params({1.731 * 4.20, sqrt(4.*M_PI*14.4)});

  std::vector<amplitude*> amps;
  amps.push_back(&amp);

  // ---------------------------------------------------------------------------
  // You shouldnt need to change anything below this line
  // ---------------------------------------------------------------------------
  double zs = cos(theta * deg2rad);
  jpacGraph1D* plotter = new jpacGraph1D();

  for (int n = 0; n < amps.size(); n++)
  {
    std::vector<double> W, dxs;
    for (int i = 0; i <= N; i++)
    {
      double Wi = (sqrt(ptr->sth) + EPS) + double(i) * (max - (sqrt(ptr->sth) + EPS)) / N;
      W.push_back(Wi);

      double xsi, si = Wi*Wi;
      if (integ == false)
      {
        xsi = amps[n]->differential_xsection(si, zs) * 1.E-3;
      }
      else
      {
        xsi = amps[n]->integrated_xsection(si) * 1.E-3;
      }

      dxs.push_back(xsi);

      debug(i, si, xsi);
    }

    plotter->AddEntry(W, dxs, amps[n]->identifier);
  }

  std::string ylabel;
  if (integ == false)
  {
    ylabel = "d#sigma/du  (#mub GeV^{-2})";
  }
  else
  {
    ylabel = "#sigma  (#mub)";
  }

  plotter->SetYaxis(ylabel, y[0], y[1]);
  plotter->SetXaxis("W  (GeV)", sqrt(ptr->sth), max);
  plotter->SetLegend(0.6, 0.65);

  plotter->Plot(filename);
}
