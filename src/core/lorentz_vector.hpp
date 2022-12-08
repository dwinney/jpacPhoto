// Custom implementation of Lorentz vectors.
// We dont actually need more properties of but we need to be able to differentiate
// them from dirac tensor objects (i.e. in spin-projection space).
// We can also use these to construct more complicated objects like rank-2 tensors
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ------------------------------------------------------------------------------

#ifndef LORENTZ_VECTOR_HPP
#define LORENTZ_VECTOR_HPP

#include <vector>

#include "constants.hpp"

namespace jpacPhoto
{
    template<int Rank>
    class lorentz_tensor;

    // This makes looping over lorentz indices more transparent
    enum class lorentz_index: int {t = 0, x = 1, y = 2, z = 3};
    const std::array<lorentz_index,4> LORENTZ_INDICES = {lorentz_index::t, lorentz_index::x, lorentz_index::y, lorentz_index::z};

    inline const double metric(lorentz_index mu)
    {
        return (mu == lorentz_index::t) ? 1. : -1.;
    };

    inline const double metric(lorentz_index mu, lorentz_index nu)
    {
        if (mu != nu) return 0;
        return (mu == lorentz_index::t) ? 1 : -1;
    };

    // Convert lorentz_index and dirac_index to ints
    inline constexpr unsigned operator+(lorentz_index x)
    {
        return static_cast<unsigned>(x);
    };

    // This allows lorentz_indices to be cout'ed
    inline std::ostream &operator<<(std::ostream &os, lorentz_index mu)
    { 
        return os << +mu;
    };    

    // ---------------------------------------------------------------------------
    // Basic object, with one lorentz index

    class lorentz_vector
    {
        public: 

        // Constructors
        lorentz_vector(){};

        lorentz_vector(std::array<complex,4> v)
        : _components(v)
        {};

        // Copy constructor
        lorentz_vector(const lorentz_vector& p)
        : _components(p._components)
        {};  

        // Copy constructor to down-cast from tensor
        template<int R=1>
        lorentz_vector(lorentz_tensor<R> const & T)
        : _components(T._indices[0]._components)
        {};  

        complex operator[](lorentz_index mu);
        inline complex operator[](int mu)
        { return operator[](static_cast<lorentz_index>(mu)); };

        lorentz_vector & operator*=(complex c);
        lorentz_vector & operator/=(complex c);
        lorentz_vector & operator+=(lorentz_vector p);
        lorentz_vector & operator-=(lorentz_vector p);
        
        lorentz_vector & operator-();
        lorentz_vector & operator=(lorentz_vector const & p);
        lorentz_vector conj();
        
        private:

        // Zero initialization
        std::array<complex,4> _components = {0., 0., 0., 0.};
    };
    
    // -----------------------------------------------------------------------
    // Interactions of Lorentz vectors and constants

    // Add and subtract
    lorentz_vector operator+(lorentz_vector lhs, lorentz_vector rhs);
    lorentz_vector operator-(lorentz_vector lhs, lorentz_vector rhs);

    // Multiply/divide by constant
    lorentz_vector operator*(complex c, lorentz_vector p);
    inline lorentz_vector operator*(lorentz_vector p, complex c){ return c *p; };
    inline lorentz_vector operator/(lorentz_vector p, complex c){ return (1./c) * p; };

    // Contractions between two vectors
    complex contract(lorentz_vector x, lorentz_vector y);
    
    // Contraction with itself
    double  square(lorentz_vector x);

    // -----------------------------------------------------------------------
    // Methods involving the LeviCivita symbol

    // Four dimensional symbol with explicit indices
    int levi_civita(lorentz_index mu, lorentz_index alpha, lorentz_index beta, lorentz_index gamma);

    // Full contractions with four vectors
    complex levi_civita(lorentz_vector a, lorentz_vector b, lorentz_vector c, lorentz_vector d);

    // Contract three indices leaving the first index open
    lorentz_vector levi_civita(lorentz_vector a, lorentz_vector b, lorentz_vector c);
};

#endif 