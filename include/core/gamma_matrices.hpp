// Header file for all things gamma matrix related
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _GAMMA_
#define _GAMMA_

#include "constants.hpp"

#include <complex>
#include <string>
#include <vector>
#include <iostream>

namespace jpacPhoto
{
    // Mostly minus metric
    const double METRIC[4] = {1., -1., -1., -1.};

	// Gamma matrix vector in Dirac basis
	const std::complex<double> GAMMA[4][4][4] =
	{
	  //gamma0
		{ { 1., 0., 0., 0. },
	    { 0., 1., 0., 0. },
	    { 0., 0., -1., 0. },
	    { 0., 0., 0., -1. } },
	  //gamma1
		{ { 0., 0., 0., 1. },
	    { 0., 0., 1., 0. },
	    { 0., -1., 0., 0. },
	    { -1., 0., 0., 0. } },
	  //gamma2
		{ { 0., 0., 0., -XI },
	    { 0., 0., XI, 0. },
	    { 0.,  XI, 0., 0. },
	    { -XI, 0., 0., 0. } },
	  //gamma3
		{ { 0., 0., 1., 0. },
	    { 0., 0., 0., -1. },
	    { -1., 0., 0., 0. },
	    { 0., 1., 0., 0. } }
	};

	// Gamma_5
	const std::complex<double> GAMMA_5[4][4] =
	{
        { 0., 0., 1., 0. },
        { 0., 0., 0., 1. },
        { 1., 0., 0., 0. },
        { 0., 1., 0., 0. }
	};

	// ---------------------------------------------------------------------------
	// Rank two gamma tensor
	std::complex<double> sigma(int mu, int nu, int i, int j);

	// ---------------------------------------------------------------------------
	// Four dimensional Levi-Civita symbol
    // We use the Peskin convention eps^{0123} = 1
	double levi_civita(int mu, int alpha, int beta, int gamma);
	
};

#endif
