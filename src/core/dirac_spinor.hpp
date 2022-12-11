// The dirac_spinor class more broadly defines a 4-component, vector-like object
//  which is assumed to have the appropriate transformations in spinor space.
//
// The usual spin-1/2 wave functions e.g. are specific instances of this more 
// general object. 
// 
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef DIRAC_SPINOR_HPP
#define DIRAC_SPINOR_HPP

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
    // First we define dirac matrices, these are rank-2 objects in spinor space
    // But are not built from spinors, rather they're operators on spinors
    // e.g. gamma matrices
    // Note these are not yet combined with lorentz structure!!

    class dirac_matrix 
    {
        public: 

        // Access a single element
        complex operator[](std::array<dirac_index,2> ij)
        {
            return _N*_entries[+ij[0]][+ij[1]];
        };

        // Assignement operations
        dirac_matrix & operator=(dirac_matrix const & G);
        dirac_matrix & operator*=(complex c);
        dirac_matrix & operator/=(complex c);

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

        // These are the fundamental dirac_matrices
        friend dirac_matrix identity();
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
    };

    // 4x4 identity matrix
    dirac_matrix identity();
    
    // Gamma matrices
    dirac_matrix gamma_0();
    dirac_matrix gamma_1();
    dirac_matrix gamma_2();
    dirac_matrix gamma_3();
    dirac_matrix gamma_5();

    // Arbitrary dirac-space operators can be built from the above as linear combinations
    dirac_matrix operator+(dirac_matrix lhs, dirac_matrix rhs);
    dirac_matrix operator-(dirac_matrix lhs, dirac_matrix rhs);

    // Add a constant automatically multiplies by the identity matrix
    inline dirac_matrix operator+(dirac_matrix lhs, complex rhs){ return lhs + rhs*identity(); };
    inline dirac_matrix operator+(complex lhs, dirac_matrix rhs){ return lhs*identity() + rhs; };
    inline dirac_matrix operator-(dirac_matrix lhs, complex rhs){ return lhs - rhs*identity(); };

    // Multiply by constant 
    dirac_matrix operator*(complex c, dirac_matrix p);
    inline dirac_matrix operator*(dirac_matrix p, complex c){ return c *p; };
    inline dirac_matrix operator/(dirac_matrix p, complex c){ return (1./c) * p; };

    // Multiply two matrices together
    dirac_matrix operator*(dirac_matrix lhs, dirac_matrix rhs);

    // -----------------------------------------------------------------------
    // The actual spinor object

    class dirac_spinor
    {
        public:

        // Initialize empty
        dirac_spinor(){};

        // Initialize with given constant components
        dirac_spinor(std::array<complex,4> entries)
        : _entries(entries)
        {};

        // Access an element
        complex operator[](dirac_index a);
        inline complex operator[](int a){return operator[](static_cast<dirac_index>(a));};

        // Operations with constants
        dirac_spinor & operator*=(complex c);
        dirac_spinor & operator/=(complex c);
        
        // Negate
        dirac_spinor operator-();

        // Adjoint
        dirac_spinor adjoint();

        private: 

        std::array<complex,4> _entries;
    };

    // Arbitrary dirac-space operators can be built from the above as linear combinations
    dirac_spinor operator+(dirac_spinor lhs, dirac_spinor rhs);
    dirac_spinor operator-(dirac_spinor lhs, dirac_spinor rhs);

    // Multiply on the right by a dirac_matrix
    dirac_spinor operator*(dirac_spinor ubar, dirac_matrix M);

    // Multiply on the left by a dirac_matrix
    dirac_spinor operator*(dirac_matrix M, dirac_spinor u);
};

#endif