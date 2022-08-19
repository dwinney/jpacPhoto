// This discontinuity relied on PW projections of the disconinuity saved to file 
// as a interpolation_2D grid
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "interpolated_discontinuity.hpp"

// ---------------------------------------------------------------------------
// Evaluate the discontinutiy by summing helicity amplitudes with corresponding d-functions

double jpacPhoto::interpolated_discontinuity::eval(double s)
{
    return 1.;
};
