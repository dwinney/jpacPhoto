// Custom implementation of Lorentz tensors.
// These are assembled by the outer/tensor products of vectors and metrics
// the two "fundamenetal" Lorentz tensors. The third structure (Levi-Civita) is defined
// as a function instead
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
    // Scalar tensors "capture" dirac data types they get multipled by a non-scalar dirac type
    // We dont use template here because compiler gets confused between complex, double, int, etc.
    // We want to always default to complex even if we write p = -2*q;


    // lorentz_tensor<complex,1> operator*(complex c, lorentz_tensor<complex,1> rhs)
    // {
    //     std::array<complex,4> new_entries = {1, 1, 1, 1};
    //     for (auto mu : LORENTZ_INDICES)
    //     {  
    //         complex old_entry = rhs(mu);
    //         new_entries[+mu]  = c * old_entry; 
    //     };

    //     return lorentz_vector<complex>(new_entries);
    // };

    // lorentz_tensor<dirac_spinor,1> operator*(dirac_spinor c, lorentz_tensor<complex,1> rhs)
    // {
    //     dirac_spinor id = identity<dirac_spinor>();
    //     std::array<dirac_spinor,4> new_entries = {id, id, id, id};
    //     for (auto mu : LORENTZ_INDICES)
    //     {  
    //         complex old_entry = rhs(mu);
    //         new_entries[+mu]  = c * old_entry; 
    //     };

    //     return lorentz_vector<dirac_spinor>(new_entries);
    // };

    // lorentz_tensor<dirac_matrix,1> operator*(dirac_matrix c, lorentz_tensor<complex,1> rhs)
    // {
    //     dirac_matrix id = identity<dirac_matrix>();
    //     std::array<dirac_matrix,4> new_entries = {id, id, id, id};
    //     for (auto mu : LORENTZ_INDICES)
    //     {  
    //         complex old_entry = rhs(mu);
    //         new_entries[+mu]  = c * old_entry; 
    //     };

    //     return lorentz_vector<dirac_matrix>(new_entries);
    // };

    // lorentz_tensor<complex,1> operator*(lorentz_tensor<complex,1> lhs, complex c)
    // {
    //     return c * lhs;
    // };

    // lorentz_tensor<dirac_matrix,1> operator*(lorentz_tensor<complex,1> lhs, dirac_matrix c)
    // {
    //     return c * lhs;
    // };

    // lorentz_tensor<dirac_spinor,1> operator*(lorentz_tensor<complex,1> lhs, dirac_spinor c)
    // {
    //     return c * lhs;
    // };
};