// Custom implementation of Lorentz tensors.
// These are assembled by the outer/tensor products of vectors
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ------------------------------------------------------------------------------

#include "lorentz_tensor.hpp"

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    // Non-member (and non-template) methods

    // Rank-0 tensors can be reconverted into scalars
    complex flatten(lorentz_tensor<0> T)
    {
        return T._N;
    };
    
    // While rank-1 can be converted into lorentz_tensors
    lorentz_vector flatten(lorentz_tensor<1> T)
    {
        return T._N*T._indices[0];
    };

    // ---------------------------------------------------------------------------
    // Contractions between tensors 

    // Contract two rank-0 tensors of the same size down to a scalar
    complex contract(lorentz_tensor<0> left, lorentz_tensor<0> right)
    {
        return flatten(left)*flatten(right);
    };
    // Contract two rank-1 tensors of the same size down to a scalar
    complex contract(lorentz_tensor<1> left, lorentz_tensor<1> right)
    {
        return contract(flatten(left), flatten(right));
    };

    // Calculate the contraction relevant for spin-1 exchanges
    complex contract(lorentz_tensor<1> bra, lorentz_tensor<2> T, lorentz_tensor<1> ket)
    {
        complex sum = 0.;
        for (auto mu : LORENTZ_INDICES)
        {
            for (auto nu : LORENTZ_INDICES)
            {
                sum += bra[{mu}] * metric(mu) * T[{mu, nu}] * metric(nu) * ket[{nu}];
            };
        };
        return sum;
    };
};