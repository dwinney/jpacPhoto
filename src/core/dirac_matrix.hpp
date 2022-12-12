// Dirac matrixes, are 4x4 (rank-2) objects in spinor space and act on spinors
// Effectively they are operators, e.g. gamma matrices 
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef DIRAC_MATRIX_HPP
#define DIRAC_MATRIX_HPP

#include <vector> 
#include <array>

#include "constants.hpp"

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    // Dirac index handling

    // Dirac space indices which get summed over
    // Naming scheme here is p,m for +- energy solutions and u,d for up/down projection
    // The names dont actually matter much however.
    // These just need to be differentiated from Lorentz indicies for transparency 
    enum dirac_index  {pu = 0, pd = 1, mu = 2, md = 3};
    const dirac_index DIRAC_INDICES[] = {pu, pd, mu, md};

    // Convert a dirac_index to its corresponding int with a plus sign
    // e.g. int x = +dirac_index(i);
    inline constexpr unsigned operator+(dirac_index x)
    {
        return static_cast<unsigned>(x);
    };

    // -----------------------------------------------------------------------
    // These are rank-2 objects in spinor space
    // Note these are not yet combined with lorentz structure!!

    class dirac_matrix 
    {
        public: 

        // Access a single element
        inline complex operator()(dirac_index i, dirac_index j)
        {
            return _N*_entries[+i][+j];
        };

        // Assignement operations
        dirac_matrix & operator=(dirac_matrix const & G);
        dirac_matrix & operator*=(complex c);
        dirac_matrix & operator/=(complex c);
        dirac_matrix & operator+=(dirac_matrix G);
        dirac_matrix & operator*=(dirac_matrix G);

        // Negation
        dirac_matrix operator-();

        // Hermition adjoint
        dirac_matrix adjoint();
        
        private:

        // Private constructor
        // rather use generating functions like gamma_0()
        // These always start with norm = 1
        dirac_matrix(std::array<std::array<complex,4>,4> entries)
        : _entries(entries), _N(1)
        {};

        dirac_matrix(complex const & c)
        : _entries(_identity), _N(c)
        {};

        // These are the fundamental dirac_matrices
        friend dirac_matrix identity_matrix();
        friend dirac_matrix gamma_0();
        friend dirac_matrix gamma_1();
        friend dirac_matrix gamma_2();
        friend dirac_matrix gamma_3();
        friend dirac_matrix gamma_5();

        // Any arbitrary matrix will be linear comibinations of these
        friend dirac_matrix operator+(dirac_matrix, dirac_matrix);
        friend dirac_matrix operator*(complex, dirac_matrix);
        friend dirac_matrix operator*(dirac_matrix, dirac_matrix);
        friend dirac_matrix operator+(dirac_matrix, complex);
        friend dirac_matrix operator+(complex, dirac_matrix);
        friend dirac_matrix operator-(dirac_matrix, complex);

        // These always have fixed size (4x4)
        std::array<std::array<complex,4>,4> _entries;

        // Store normalization seperately to avoid unnecessary calculation
        complex _N = 1;

        static constexpr std::array<std::array<complex,4>,4>  _identity = {{ { 1,  0,  0,  0},
                                                                             { 0,  1,  0,  0},
                                                                             { 0,  0,  1,  0},
                                                                             { 0,  0,  0,  1}  }};
    };

    // ---------------------------------------------------------------------------
    // "Constructor" functions

    // 4x4 identity_matrix matrix
    dirac_matrix identity_matrix();
    
    // Gamma matrices
    dirac_matrix gamma_0();
    dirac_matrix gamma_1();
    dirac_matrix gamma_2();
    dirac_matrix gamma_3();
    dirac_matrix gamma_5();

    // Dirac_matrix filled with NaN's for error throwing
    template<>
    inline dirac_matrix NaN()
    {
        return NaN<complex>() * identity_matrix();
    };

    // ---------------------------------------------------------------------------
    // Interactions between matrices

    // Access conjugate from outside
    inline dirac_matrix conj(dirac_matrix x){ return x.adjoint(); };
    
    // Arbitrary dirac-space operators can be built from the above as linear combinations
    dirac_matrix operator+(dirac_matrix lhs, dirac_matrix rhs);
    dirac_matrix operator-(dirac_matrix lhs, dirac_matrix rhs);

    // Add a constant automatically multiplies by the identity_matrix matrix
    inline dirac_matrix operator+(dirac_matrix lhs, complex rhs){ return lhs + rhs*identity_matrix(); };
    inline dirac_matrix operator+(complex lhs, dirac_matrix rhs){ return lhs*identity_matrix() + rhs; };
    inline dirac_matrix operator-(dirac_matrix lhs, complex rhs){ return lhs - rhs*identity_matrix(); };

    // Multiply by constant 
    dirac_matrix operator*(complex c, dirac_matrix p);
    inline dirac_matrix operator*(dirac_matrix p, complex c){ return c *p; };
    inline dirac_matrix operator/(dirac_matrix p, complex c){ return (1./c) * p; };

    // Multiply two matrices together
    dirac_matrix operator*(dirac_matrix lhs, dirac_matrix rhs);
};

#endif