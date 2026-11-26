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

#include "jpacUtils/math.hpp"

namespace jpacPhoto
{
    using jpacUtils::complex;

    // ---------------------------------------------------------------------------
    // Define 1, 0, and NaN for jpacPhoto custom data types
    
    template<typename T>
    inline T NaN(){ return std::numeric_limits<T>::quiet_NaN(); };

    template<>
    inline complex NaN(){ return complex(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()); };

    template<typename T>
    T zero();

    template<>
    inline complex zero() { return 0; };

    template<typename T> 
    T identity();

    template<>
    inline complex identity() { return 1; };

    // ---------------------------------------------------------------------------
    // PDG Meson masses in GeV

    const double M_PION      = 0.13957000;
    const double M_KAON      = 0.49367700;
    const double M_ETA       = 0.54753;
    const double M_RHO       = 0.77526;
    const double M_OMEGA     = 0.78265;
    const double M_PHI       = 1.01956;
    const double M_A2        = 1.3186;
    const double M_F2        = 1.2754;
    const double M_B1        = 1.229;
    const double M_JPSI      = 3.0969160;
    const double M_PSI2S     = 3.686;
    const double M_D         = 1.86484;
    const double M_DSTAR     = 2.00685;
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

    // ---------------------------------------------------------------------------
    // Out of the box std::complex<double> doesnt play well with ints. Here we explicitly
    // give it the functionality we want

    // Additionally the complex type is a liitle dim in c++ so we need to define int & bool multiplication
    inline complex operator*(const int& c, const complex& rhs)
    {
        return complex(c*rhs.real(), c*rhs.imag());
    };

    inline complex operator*(const complex& lhs, const int& c)
    {
        return complex(c*lhs.real(), c*lhs.imag());
    };

    inline complex operator*(const bool& c, const complex& rhs)
    {
        return (c) ? rhs : 0.;
    };

    inline complex operator*(const complex& lhs, const bool& c)
    {
        return (c) ? lhs : 0.;
    };

    inline complex operator/(const complex&c, const int& i)
    {
        return (1./i)*c;
    };

    inline complex operator/(const int& i, const complex&c)
    {
        return (1./c)*i;
    };

    inline complex operator+(const complex&c, const int& i)
    {
        return c + complex(1,0)*i;
    };

    inline complex operator+(const int& i, const complex & c)
    {
        return complex(1,0)*i + c;
    };

    inline complex operator-(const complex&c, const int& i)
    {
        return c - complex(1,0)*i;
    };

    inline complex operator-(const int& i, const complex & c)
    {
        return complex(1,0)*i - c;
    };

    inline bool operator == (const complex &z,const int n)
    {
        return (z == static_cast<double> (n));
    }
    inline bool operator == (const int n, const complex &z)
    {
        return (z == static_cast<double> (n));
    }
    inline bool operator != (const complex &z,const int n)
    {
        return (z != static_cast<double> (n));
    }
    inline bool operator != (const int n, const complex &z)
    {
        return (z != static_cast<double> (n));
    }
};
// ---------------------------------------------------------------------------

#endif
