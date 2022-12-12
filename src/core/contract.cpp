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
#include "lorentz_tensor.hpp"

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    // Define contract between the lorentz-scalar types first

    complex contract(dirac_spinor left, dirac_spinor right)
    {
        complex sum = 0.;
        for (auto i : DIRAC_INDICES)
        {
            sum += left(i) * right(i);
        };
        return sum;
    };

    // Again, spinors need to be defined seperately becasue they dont return their own type
    complex contract(lorentz_tensor<dirac_spinor,1> left, lorentz_tensor<dirac_spinor,1> right)
    {
        complex sum = 0.;
        for (auto mu : LORENTZ_INDICES)
        {
            sum += metric(mu) * contract(left(mu), right(mu));
        };
        return sum;
    };

    // ---------------------------------------------------------------------------
    // Function to return a vector of all permutations of N indices

    std::vector<std::vector<lorentz_index>> permutations(unsigned N)
    {
        if (N == 1) return { {lorentz_index::t}, {lorentz_index::x}, 
                                                 {lorentz_index::y}, 
                                                 {lorentz_index::z}};

        // Get the previous set of permutations
        std::vector<std::vector<lorentz_index>> previous = permutations(N - 1);
        std::vector<std::vector<lorentz_index>> next;
        // and add the next instance
        for (auto single_prev : previous )
        {
            for (auto mu : LORENTZ_INDICES)
            {
                std::vector<lorentz_index> single_next(single_prev.begin(), single_prev.end());
                single_next.push_back(mu);
                next.push_back(single_next);
            };
        };
        return next;
    };
    
    // Single function to produce the metric along the diagonal
    int metric(lorentz_index mu)
    {
        return (mu == lorentz_index::t) ? 1 : -1;
    }; 

    int metric(std::vector<lorentz_index> permutations)
    {
        int prod = 1;
        for (auto mu : permutations)
        {
            if (mu != lorentz_index::t) prod *= -1;
        };
        return prod;
    };
};