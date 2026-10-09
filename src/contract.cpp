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

    // ---------------------------------------------------------------------------
    // Function to return a vector of all permutations of N indices

    std::vector<std::vector<lorentz_index>> permutations(unsigned N)
    {
        auto t = LORENTZ_INDICES[0];
        auto x = LORENTZ_INDICES[1];
        auto y = LORENTZ_INDICES[2];
        auto z = LORENTZ_INDICES[3];

        // First two cases are hard-coded for easy access
        if (N == 1) return { {t}, {x}, {y}, {z}};

        if (N == 2) return { {t, t}, {t, x}, {t, y}, {t, z}, 
                             {x, t}, {x, x}, {x, y}, {x, z}, 
                             {y, t}, {y, x}, {y, y}, {y, z}, 
                             {z, t}, {z, x}, {z, y}, {z, z} };

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

    // Return the value of the levi-civita symbol for some combination of lorentz_indices
    int levi_civita(lorentz_index mu, lorentz_index nu, lorentz_index alpha, lorentz_index beta)
    {
        // Convert indices to their ints
        int a = +mu, b = +nu, c = +alpha, d = +beta;
        int result = (d - c) * (d - b) * (d - a) * (c - b) * (c - a) * (b - a);
        return (result == 0) ? result : result / std::abs(result);
    };

    // ---------------------------------------------------------------------------
    // Covariant gamma matrix structures

    // Non-trivial dirac_matrix tensors
    // Any arbitrary tensor up to rank-2 can be built from these

    inline lorentz_tensor<dirac_matrix,1> gamma_vector()
    { 
        return lorentz_vector<dirac_matrix>({gamma_0(), gamma_1(), gamma_2(), gamma_3()});
    };
    inline lorentz_tensor<dirac_matrix,2> sigma_tensor()
    { 
        return tensor_product(gamma_vector(), gamma_vector()) - identity<dirac_matrix>()*metric_tensor(); 
    };

    inline dirac_matrix slash(lorentz_tensor<complex,1> q)
    {
        auto gamma = gamma_vector();
        dirac_matrix sum = zero<dirac_matrix>();
        for (auto mu : LORENTZ_INDICES)
        {
            sum += metric(mu) * gamma(mu) * q(mu);
        }
        return sum;
    };

    // ---------------------------------------------------------------------------
    // Special case of lorentz_tensor which mixes spinors and matrices. 
    // Here the saved subtensors, normalizations, and return values are different types
    
    template<int Rank>
    class bilinear_tensor : public tensor_object<complex>
    {
        public: 

        // Default tensor with nothing initialized
        bilinear_tensor(){};

        // Copy constructor
        bilinear_tensor(bilinear_tensor<Rank> const & old)
        : _matrix(old._matrix),
          _lhs(old._lhs), _rhs(old._rhs)
        {};

        // Implicit constructor, stores pointers to constituent tensors of smaller rank
        bilinear_tensor(dirac_spinor L, lorentz_tensor<dirac_matrix, Rank> Ts, dirac_spinor R)
        :   _lhs(L), _rhs(R),
            _matrix(Ts)
        {};

        inline complex operator()(std::vector<lorentz_index> indices)
        {
            if (indices.size() != Rank) return error("lorentz_tensor - Incorrect number of indices passed!", NaN<complex>());          

            // begin producting all the matrices to get one 
            dirac_matrix M = _matrix(indices);
            return contract(_lhs, M*_rhs); 
        };

        // Get the rank of the tensor (number of open indicies)
        inline unsigned int rank() const{ return Rank; };

        protected:

        // These are always treated as if they were tensor products 
        lorentz_tensor<dirac_matrix, Rank> _matrix;

        // Normalizations for the spinors on either side
        dirac_spinor _lhs = zero<dirac_spinor>();
        dirac_spinor _rhs = zero<dirac_spinor>();
    };

    // ---------------------------------------------------------------------------
    // Bilinear methods combine two dirac_spiors with a lorentz_tensor 

    template<int R>
    inline lorentz_tensor<complex,R> bilinear(dirac_spinor ubar, lorentz_tensor<dirac_matrix,R> Gamma, dirac_spinor u)
    {
        std::shared_ptr<tensor_object<complex>> m_ptr = std::make_shared<bilinear_tensor<R>>(ubar, Gamma, u);
        
        return lorentz_tensor<complex,R>({m_ptr}, false);
    };
};