// The contract() function assembles different tensors into other structures
// At present these different interacitons must all be specified individually
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef CONTRACT_HPP
#define CONTRACT_HPP

#include "lorentz_vector.hpp"
#include "lorentz_tensor.hpp"
#include "dirac_spinor.hpp"

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    // Contractions between tensors

    // Contractions of multiple tensors into scalars
    complex contract(lorentz_tensor<0> left, lorentz_tensor<0> right);
    complex contract(lorentz_tensor<1> left, lorentz_tensor<1> right);
    complex contract(lorentz_tensor<1> bra,  lorentz_tensor<2> T, lorentz_tensor<1> ket);

};

#endif