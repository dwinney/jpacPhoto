// Phenomenological expressions for the total cross-sections of different processes
// Defined as global static functions to be used inside triple_regge amplitudes
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef SIGMA_TOT
#define SIGMA_TOT

#include <cmath>

inline double sigmatot_pi(double s)
{
    double result = 13.63 * pow(s, 0.0808) + 31.79 * pow(s, -0.425);
    return result; // convert to nb
};

#endif