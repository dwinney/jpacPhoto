// Header file with global phyiscal constants.
// Everything is in GeV unless explicitly stated otherwise.
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _DEBUG_
#define _DEBUG_

#include <iostream>
#include <iomanip>

namespace jpacPhoto
{
    // little function for printing to screen instead of having to copy this line all the time
    template<typename T>
    inline void debug(T x)
    {
        std::cout << x << std::endl;
    };

    template<typename T, typename F>
    inline void debug(T x, F y)
    {
        std::cout << std::left << std::setw(15) << x;
        std::cout << std::left << std::setw(15) << y << std::endl;
    };

    template<typename T, typename F, typename G>
    inline void debug(T x, F y, G z)
    {
        std::cout << std::left << std::setw(15) << x;
        std::cout << std::left << std::setw(15) << y;
        std::cout << std::left << std::setw(15) << z << std::endl;
    };
};

#endif

#ifndef CONSTANT
#define CONSTANT

#include <cmath>
#include <complex>

namespace jpacPhoto
{
    // Spin-parity combinations
    constexpr std::array<int,2> SCALAR{{0, +1}};
    constexpr std::array<int,2> VECTOR{{1, -1}};
    constexpr std::array<int,2> PSEUDO_SCALAR{{0, -1}};
    constexpr std::array<int,2> AXIAL_VECTOR{{1, +1}};

    enum HELICITY_CHANNEL{ S, T, U };

    // ---------------------------------------------------------------------------
    const double PI       = M_PI;
    const double DEG2RAD  = (M_PI / 180.);
    const double EPS      = 1.e-6;
    const double ALPHA    = 1. / 137.;
    const double E        = sqrt(4. * PI * ALPHA);

    const std::complex<double> XR(1., 0.);
    const std::complex<double> XI(0., 1.);
    const std::complex<double> IEPS(0., EPS);

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
    const double M_LAMBDAC   = 2.28646;

    // Baryon masses squared
    const double M2_PROTON    = M_PROTON * M_PROTON;
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

};
// ---------------------------------------------------------------------------

#endif
