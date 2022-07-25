// Header file for all things gamma matrix related
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "gamma_matrices.hpp"

// ---------------------------------------------------------------------------
// Rank two gamma tensor
std::complex<double> jpacPhoto::sigma(int mu, int nu, int i, int j)
{
    std::complex<double> result = 0.;
    for (int k = 0; k < 4; k++)
    {
        result += GAMMA[mu][i][k] * GAMMA[nu][k][j];
        result -= GAMMA[nu][i][k] * GAMMA[mu][k][j];
    }

    result *= XR / 2.;

    return result;
};

// ---------------------------------------------------------------------------
// Four dimensional Levi-Civita symbol
double jpacPhoto::levi_civita(int a, int b, int c, int d)
{
    int result = (d - c) * (d - b) * (d - a) * (c - b) * (c - a) * (b - a);

    if (result == 0) return 0.;

    result /= abs(d - c) * abs(d - b) * abs(d - a) * abs(c - b) * abs(c - a) * abs(b - a);

    return result;
};
