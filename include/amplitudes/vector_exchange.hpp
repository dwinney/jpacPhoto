// Axial-vector meson photoproduction proceeding through a vector meson exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _AXIAL_
#define _AXIAL_

#include <string>
#include <vector>
#include <iostream>

#include <iomanip>

#include "amplitude.hpp"
#include "gamma_technology.hpp"

// ---------------------------------------------------------------------------
// vector_exchange class describes the amplitude for a fixed-spin-1 exchange
// in the t-channel. Derived in terms of simple feynman rules at tree level
//
// Initialization required a reaction_kinematics object, the mass of the exchange,
// and an optional string to identify the amplitude with.
//
//  Evaluation requires three couplings photon coupling, gGamma, and vector/tensor
// nucleon couplings, gV and gT respectively.
//
// Set couplings with amp.set_params({gGamma, gV, gT});
// ---------------------------------------------------------------------------

class vector_exchange : public amplitude
{
public:
  // Constructor
  vector_exchange(reaction_kinematics * xkinem, double mass, std::string exchange = "")
  : amplitude(xkinem, exchange), mEx2(mass*mass)
  {};

  // Copy constructor
  vector_exchange(const vector_exchange & old)
  : amplitude(old), mEx2(old.mEx2),
    gGamma(old.gGamma), gV(old.gV), gT(old.gT)
  {};

  // Setting utility
  void set_params(std::vector<double> params)
  {
    gGamma = params[0];
    gV = params[1];
    gT = params[2];
  };

  // Assemble the helicity amplitude by contracting the lorentz indices
  std::complex<double> helicity_amplitude(std::vector<int> helicities, double s, double zs);

private:
  // Mass of the exchange
  double mEx2;

  // Couplings to the axial-vector/photon and vector/tensor couplings to nucleon
  double gGamma = 0., gV = 0., gT = 0.;

  // Four-momentum of the exhange
  std::complex<double> exchange_momenta(int mu, double s, double zs);

  // Photon - Axial Vector - Vector vertex
  std::complex<double> top_vertex(int mu, int lam_gam, int lam_vec, double s, double zs);

  // Nucleon - Nucleon - Vector vertex
  std::complex<double> bottom_vertex(int nu, int lam_targ, int lam_rec, double s, double zs);

  // Vector propogator
  std::complex<double> vector_propagator(int mu, int nu, double s, double zs);
};

#endif