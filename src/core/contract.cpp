// The contract() function assembles different tensors into other structures
// At present these different interacitons must all be specified individually
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#include "contract.hpp"

namespace jpacPhoto
{
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