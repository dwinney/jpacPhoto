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

#include "dirac_spinor.hpp"

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    // Intrinsic operators of dirac_matrices

    // Assignment and scaled re-assignment 

    dirac_matrix & dirac_matrix::operator=(dirac_matrix const & G)
    {
        _N       = G._N;
        _entries = G._entries;
        return *this;
    };

    dirac_matrix & dirac_matrix::operator*=(complex c)
    {
        _N *= c;
        return *this;  
    };

    dirac_matrix & dirac_matrix::operator/=(complex c)
    {
        _N /= c;
        return *this;  
    };

    // Negation
    dirac_matrix dirac_matrix::operator-()
    {
        dirac_matrix neg = *this;
        neg *= -1;
        return neg;
    }; 
    
    // Return the hermitian adjoint of operator;
    dirac_matrix dirac_matrix::adjoint()
    {
        std::array<std::array<complex,4>,4> entries;

        for (auto i : DIRAC_INDICES)
        {
            for (auto j : DIRAC_INDICES)
            {
                entries[+i][+j] = std::conj(_entries[+j][+i]);
            };
        };

        return dirac_matrix(entries);
    };

    // ---------------------------------------------------------------------------
    // Dirac-space matrices

    dirac_matrix identity()
    {
        // Note need to use triple brace initialization here annoyingly
        return dirac_matrix({{ { 1,  0,  0,  0},
                               { 0,  1,  0,  0},
                               { 0,  0,  1,  0},
                               { 0,  0,  0,  1}  }});
    };

    dirac_matrix gamma_0()
    {
        return dirac_matrix({{ { 1,  0,   0,   0},
                               { 0,  1,   0,   0},
                               { 0,  0,  -1,   0},
                               { 0,  0,   0,  -1}  }});
    };

    dirac_matrix gamma_1()
    {
        return dirac_matrix({{ { 0,  0,   0,   1},
                               { 0,  0,   1,   0},
                               { 0, -1,   0,   0},
                               {-1,  0,   0,   0}  }});
    };

    dirac_matrix gamma_2()
    {
        return dirac_matrix({{ { 0,  0,   0,   I},
                               { 0,  0,   I,   0},
                               { 0, -I,   0,   0},
                               {-I,  0,   0,   0}  }});
    };

    dirac_matrix gamma_3()
    {
        return dirac_matrix({{ { 0,  0,   1,   0},
                               { 0,  0,   0,  -1},
                               {-1,  0,   0,   0},
                               { 0,  1,   0,   0}  }});
    };

    dirac_matrix gamma_5()
    {
        return dirac_matrix({{ { 0,  0,   1,   0},
                               { 0,  0,   0,   1},
                               { 1,  0,   0,   0},
                               { 0,  1,   0,   0}  }});
    };

    // ---------------------------------------------------------------------------
    // Linear combination operators to build arbitrarty dirac operator

    // Add two matrices together
    dirac_matrix operator+(dirac_matrix lhs, dirac_matrix rhs)
    {
        std::array<std::array<complex,4>,4> entries;

        for (auto i : DIRAC_INDICES)
        {
            for (auto j : DIRAC_INDICES)
            {
                entries[i][j] = lhs[{i,j}] + rhs[{i,j}];
            };
        };

        return dirac_matrix(entries);
    };

    // Subtract two matrices
    dirac_matrix operator-(dirac_matrix lhs, dirac_matrix rhs)
    {
        return lhs + (-rhs);
    };

    // Multiply be a constant
    dirac_matrix operator*(complex c, dirac_matrix T)
    {
        dirac_matrix cT = T;
        cT._N *= c;
        return cT;
    };

    // Multiply two matrices together
    dirac_matrix operator*(dirac_matrix lhs, dirac_matrix rhs)
    {
        std::array<std::array<complex,4>,4> entries;

        for (auto i : DIRAC_INDICES)
        {
            for (auto j : DIRAC_INDICES)
            {
                complex sum = 0.;
                for (auto k : DIRAC_INDICES)
                {
                    sum += lhs[{i,k}]*rhs[{k,j}];
                }
                entries[+i][+j] = sum;
            };
        };

        return dirac_matrix(entries);
    };

    // ---------------------------------------------------------------------------
    // Intrinsic properties of dirac_spinors

    // Access elements
    complex dirac_spinor::operator[](dirac_index a)
    {
        return _entries[+a];
    };

    // Scaled re-assignement
    dirac_spinor & dirac_spinor::operator*=(complex c)
    {
        for (auto i : DIRAC_INDICES)
        {
            _entries[+i] *= c;
        }
        return *this;
    }

    dirac_spinor & dirac_spinor::operator/=(complex c)
    {
        for (auto i : DIRAC_INDICES)
        {
            _entries[+i] /= c;
        }
        return *this;
    };

    // Negate a spinor
    dirac_spinor dirac_spinor::operator-()
    {
        dirac_spinor neg = *this;
        for (auto i : DIRAC_INDICES)
        {
            neg._entries[+i] *= -1;
        }
        return neg;
    };

    // Return the adjoint

    // ---------------------------------------------------------------------------
    // Non-member operations of dirac_spinors

    //Add two spinors together
    dirac_spinor operator+(dirac_spinor lhs, dirac_spinor rhs)
    {
        std::array<complex,4> entries;
        for (auto i : DIRAC_INDICES)
        {
            entries[+i] = lhs[i] + rhs[i];
        };
        return dirac_spinor(entries);
    };

    // Subtract two spinors 
    dirac_spinor operator-(dirac_spinor lhs, dirac_spinor rhs)
    {
        return lhs + (-rhs);
    };

    // Multiply on the right by a dirac_matrix
    dirac_spinor operator*(dirac_spinor ubar, dirac_matrix M)
    {
        std::array<complex,4> entries;
        for (auto i : DIRAC_INDICES)
        {
            complex sum = 0.;
            for (auto j : DIRAC_INDICES)
            {
                sum += ubar[j] * M[{j,i}]
            };
            entries[+i] = sum;
        };

        return dirac_spinor(entries);
    };

    // Multiply on the left by a dirac_matrix
    dirac_spinor operator*(dirac_matrix M, dirac_spinor u)
    {
        std::array<complex,4> entries;
        for (auto i : DIRAC_INDICES)
        {
            complex sum = 0.;
            for (auto j : DIRAC_INDICES)
            {
                sum += M[{i,j}] * u[j];
            };
            entries[+i] = sum;
        };

        return dirac_spinor(entries);
    };


};