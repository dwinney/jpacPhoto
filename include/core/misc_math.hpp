// Misc math functions that are useful to have
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _MATH_
#define _MATH_

#include "constants.hpp"
#include <iostream>
#include <complex>
#include <algorithm>

namespace jpacPhoto
{
    template <typename T>
    inline T Kallen(T x, T y, T z)
    {
        return x*x + y*y + z*z - 2. * (x*y + x*z + y*z);
    };

    std::complex<double> cgamma(std::complex<double> z, int OPT = 0);

    inline unsigned int factorial(unsigned int n) 
    {
        if (n == 0)
        return 1;
        return n * factorial(n - 1);
    };

    // ---------------------------------------------------------------------------
    // Wigner d-func coefficient of leading power
    double wigner_leading_coeff(int j, int lam1, int lam2);

    // Wigner d-function for half-integer spin
    double wigner_d_half(int j, int lam1, int lam2, double theta);

    // Wigner d-function for integer spin
    double wigner_d_int(int j, int lam1, int lam2, double theta);

    // Wigner d-function for integer spin in terms of the cosine of theta not theta
    std::complex<double> wigner_d_int_cos(int j, int lam1, int lam2, double cos);

    // Legendre function in terms of cosine theta
    double legendre(int l, double z);
};

#endif
