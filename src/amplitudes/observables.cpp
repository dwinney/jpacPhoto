// Abstract class for an amplitude. Used so we can easily build observables
// as the incoherent sum of amplitudes in s, t, and u channels.
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/amplitude.hpp"

// ---------------------------------------------------------------------------
// Differential cross section dsigma / dt
// in NANOBARN
double amplitude::differential_xsection(double s, double zs)
{
  double sum = 0.;
  for (int i = 0; i < 24; i++)
  {
    std::complex<double> amp_i = helicity_amplitude(kinematics->helicities[i], s, zs);
    sum += std::real(amp_i * conj(amp_i));
  }

  double norm = 1.;
  norm /= 64. * M_PI * s;
  norm /= real(pow(kinematics->initial.momentum("beam", s), 2.));
  norm /= (2.56819E-6); // Convert from GeV^-2 -> nb

  return norm * sum;
};

// ---------------------------------------------------------------------------
// Inegrated total cross-section
// IN NANOBARN
double amplitude::integrated_xsection(double s)
{
  auto F = [&](double zs)
  {
    double jacobian = 2.; // 2. * k * q
    jacobian *= real(kinematics->initial.momentum("beam", s));
    jacobian *= real(kinematics->final.momentum(kinematics->vector_particle, s));

    return differential_xsection(s, zs) * jacobian;
  };

  double result;
  result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(F, -1, 1, 0, 1.E-9, NULL);

  return result;
};

// ---------------------------------------------------------------------------
// Polarizatiopn asymmetry between beam and recoil proton
double amplitude::K_LL(double s, double zs)
{
  double sigmapp = 0., sigmapm = 0.;
  for (int i = 0; i < 6; i++)
  {
    std::complex<double> squarepp, squarepm;

    // Amplitudes with lam_gam = + and lam_recoil = +
    squarepp = helicity_amplitude(kinematics->helicities[2*i+1], s, zs);
    squarepp *= conj(squarepp);
    sigmapp += real(squarepp);

    // Amplitudes with lam_gam = + and lam_recoil = -
    squarepm = helicity_amplitude(kinematics->helicities[2*i], s, zs);
    squarepm *= conj(squarepm);
    sigmapm += real(squarepm);
  }

  return (sigmapp - sigmapm) / (sigmapp + sigmapm);
}

// ---------------------------------------------------------------------------
// Polarizatiopn asymmetry between beam and target proton
double amplitude::A_LL(double s, double zs)
{
  double sigmapp = 0., sigmapm = 0.;
  for (int i = 0; i < 6; i++)
  {
    std::complex<double> squarepp, squarepm;

    // Amplitudes with lam_gam = + and lam_targ = +
    squarepp = helicity_amplitude(kinematics->helicities[i+6], s, zs);
    squarepp *= conj(squarepp);
    sigmapp += real(squarepp);

    // Amplitudes with lam_gam = + and lam_targ = -
    squarepm = helicity_amplitude(kinematics->helicities[i], s, zs);
    squarepm *= conj(squarepm);
    sigmapm += real(squarepm);
  }

  return (sigmapp - sigmapm) / (sigmapp + sigmapm);
}