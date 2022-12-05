// Header file with global phyiscal constants.
// Everything is in GeV unless explicitly stated otherwise.
//
// ---------------------------------------------------------------------------
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef CONSTANT_HPP
#define CONSTANT_HPP

#include <cmath>
#include <complex>

#include "debug.hpp"

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    // Mathematical constants 

    // We use complex numbers everywhere throughout so I'll define this shortened data type
    using complex = std::complex<double>;

    const double PI       = M_PI;
    const double DEG2RAD  = (M_PI / 180.);
    const double EPS      = 1.e-6;
    const double ALPHA    = 1. / 137.;
    const double E        = sqrt(4. * PI * ALPHA);

    const complex XR  (1., 0.);
    const complex XI  (0., 1.);
    const complex IEPS(0., EPS);


    // PDG Meson masses in GeV
    const double M_PION      = 0.13957000;
    const double M_KAON      = 0.49367700;
    const double M_ETA       = 0.54753;
    const double M_RHO       = 0.77526;
    const double M_OMEGA     = 0.78265;
    const double M_PHI       = 1.01956;
    const double M_B1        = 1.229;
    const double M_JPSI      = 3.0969160;
    const double M_PSI2S     = 3.686;
    const double M_D         = 1.86965;
    const double M_DSTAR     = 2.01026;
    const double M_UPSILON1S = 9.4603;
    const double M_UPSILON2S = 10.02336;
    const double M_UPSILON3S = 10.3552;
    const double M_CHIC1     = 3.51067;

    // Exotic Meson Masses
    const double M_X3872     = 3.87169;
    const double M_Y4260     = 4.220;
    const double M_ZC3900    = 3.8884;
    const double M_ZB10610   = 10.6072;
    const double M_ZB10650   = 10.6522;

    // Meson masses squared
    const double M2_PION     = M_PION * M_PION;
    const double M2_JPSI     = M_JPSI * M_JPSI;
    const double M2_D        = M_D * M_D;
    const double M2_DSTAR    = M_DSTAR * M_DSTAR; 
    const double M2_B1       = M_B1 * M_B1;

    // Baryon masses
    const double M_PROTON    = 0.938272;
    const double M_LAMBDA    = 1.1157;
    const double M_DELTA     = 1.23055;
    const double M_LAMBDAC   = 2.28646;

    // Baryon masses squared
    const double M2_PROTON    = M_PROTON * M_PROTON;
    const double M2_LAMBDA    = M_LAMBDA * M_LAMBDA;
    const double M2_LAMBDAC   = M_LAMBDAC * M_LAMBDAC;

    // Decay constants in GeV
    const double F_JPSI      = 0.278;
    const double F_UPSILON1S = 0.23345;
    const double F_UPSILON2S = 0.16563;
    const double F_UPSILON3S = 0.1431;

    // Photon lab energy
    inline double E_beam(double W)
    {
        return (W*W / M_PROTON - M_PROTON) / 2.;
    };

    // Center of mass energy given beam energy
    inline double W_cm(double egam)
    {
        return sqrt(M_PROTON * (2. * egam + M_PROTON));
    };

    // ---------------------------------------------------------------------------
    // Kallen Triangle function

    // Only way to get a double or int Kallen is if all inputs are double/int
    template<typename T>
    inline T Kallen(T x, T y, T z)
    {
        return x*x + y*y + z*z - 2. * (x*y + x*z + y*z);
    };

    // If any of them are complex, return complex
    inline complex Kallen(complex z, double a, double b) { return Kallen<complex>(z, XR*a, XR*b); };
    inline complex Kallen(double a, complex z, double b) { return Kallen<complex>(XR*a, z, XR*b); };
    inline complex Kallen(double a, double b, complex z) { return Kallen<complex>(XR*a, XR*b, z); };

};
// ---------------------------------------------------------------------------

#endif
