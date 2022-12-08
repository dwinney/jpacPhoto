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

#include "lorentz_vector.hpp"

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    // Intrinsic methods and properties of lorentz_vectors
    
    // Access the mu-th element
    complex lorentz_vector::operator[](lorentz_index mu)
    {
        return _components[+mu]; 
    };

    lorentz_vector & lorentz_vector::operator+=(lorentz_vector p)
    {
        for (auto mu : LORENTZ_INDICES)
        {
            _components[+mu] += p[mu];
        };
        return *this;
    };

    lorentz_vector & lorentz_vector::operator-=(lorentz_vector p)
    {
        *this += (-p);
        return *this;
    };

    lorentz_vector & lorentz_vector::operator*=(complex c)
    {
        for (int i = 0; i < 4; i++)
        {
            _components[i] *= c;
        };
        return *this;
    };

    lorentz_vector & lorentz_vector::operator/=(complex c)
    {
        *this *= (1./c);
        return *this;
    };

    // Negate a 4-vector
    lorentz_vector & lorentz_vector::operator-()
    {
        *this *= -1.;
        return *this;
    };

    // Define new vectors from old ones
    lorentz_vector & lorentz_vector::operator=(lorentz_vector const & p)
    {
        _components = p._components;
        return *this;
    };

    // Return complex conjugate of a vector
    lorentz_vector lorentz_vector::conj()
    {
        lorentz_vector result;
        for (auto mu : LORENTZ_INDICES)
        {
            result._components[+mu] = std::conj(_components[+mu]);
        };
        return result;
    };

    // ---------------------------------------------------------------------------
    // Non-class methods related to lorentz_vector

    // Add two four-vectors together
    lorentz_vector operator+(lorentz_vector lhs, lorentz_vector rhs)
    {
        std::array<complex,4> result;
        for (auto mu : LORENTZ_INDICES)
        {
            result[+mu] = lhs[mu] + rhs[mu];
        };
        return lorentz_vector(result);
    };

    lorentz_vector operator-(lorentz_vector lhs, lorentz_vector rhs)
    {
        return lhs + (-rhs);
    };

    // Multiply a vector by a scalar
    lorentz_vector operator*(complex c, lorentz_vector p)
    {
        std::array<complex,4> cp;
        for (auto mu : LORENTZ_INDICES)
        {
            cp[+mu] = c * p[mu];
        };
        return lorentz_vector(cp);
    };

    // Contract two vectors together
    complex contract(lorentz_vector x, lorentz_vector y)
    {
        complex sum = 0.;
        for (auto mu : LORENTZ_INDICES)
        {
            sum += x[mu] * metric(mu) * y[mu];
        };
        return sum;
    };  

    // Contract a 4-vector with its conjugate
    double square(lorentz_vector x)
    {
        double sum = 0.;
        for (auto mu : LORENTZ_INDICES)
        {
            sum += std::norm(x[mu]) * metric(mu);
        };
        return sum;
    };  

    // ---------------------------------------------------------------------------
    // Methods involving the LeviCivita symbol

    // Four dimensional symbol with explicit indices
    int levi_civita(lorentz_index mu, lorentz_index alpha, lorentz_index beta, lorentz_index gamma)
    {
        int a = +mu;
        int b = +alpha;
        int c = +beta;
        int d = +gamma;

        int result = (d - c) * (d - b) * (d - a) * (c - b) * (c - a) * (b - a);
        if (result == 0) return 0.;
        result /= abs(d - c) * abs(d - b) * abs(d - a) * abs(c - b) * abs(c - a) * abs(b - a);
        return - result;
    };

    // Full contractions with four vectors
    complex levi_civita(lorentz_vector a, lorentz_vector b, lorentz_vector c, lorentz_vector d)
    {
        // Sum over all combinations of indices
        complex result = 0.;
        for (auto mu : LORENTZ_INDICES)
        {
            for (auto nu : LORENTZ_INDICES)
            {
                for (auto alpha : LORENTZ_INDICES)
                {
                    for (auto beta : LORENTZ_INDICES)
                    {
                        int eps = levi_civita(mu, nu, alpha, beta);
                        if (eps == 0) continue;

                        result += eps * a[mu] * b[nu] * c[alpha] * d[beta];
                    };
                };
            };
        };

        return result;
    };

    // Contract three indices leaving the first index open
    lorentz_vector levi_civita(lorentz_vector a, lorentz_vector b, lorentz_vector c)
    {
        complex t = levi_civita(lorentz_vector({1.,0.,0.,0.}), a , b, c);
        complex x = levi_civita(lorentz_vector({0.,1.,0.,0.}), a , b, c);
        complex y = levi_civita(lorentz_vector({0.,0.,1.,0.}), a , b, c);
        complex z = levi_civita(lorentz_vector({0.,0.,0.,1.}), a , b, c);
        return lorentz_vector({t, x, y, z});
    };
};