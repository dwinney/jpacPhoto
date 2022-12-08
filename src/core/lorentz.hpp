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

#ifndef LORENTZ_HPP
#define LORENTZ_HPP

#include <vector>
#include <initializer_list>

#include "constants.hpp"

namespace jpacPhoto
{
    // This makes looping over lorentz indices more transparent
    enum class lorentz_index: int {t = 0, x = 1, y = 2, z = 3};
    const std::array<lorentz_index,4> LORENTZ_INDICES = {lorentz_index::t, lorentz_index::x, lorentz_index::y, lorentz_index::z};

    inline const double metric(lorentz_index mu)
    {
        return (mu == lorentz_index::t) ? 1. : -1.;
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

    class lorentz_vector;

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

    // // -----------------------------------------------------------------------
    // // Now, we can think of objects with more than one index.
    // // In most useful contexts, these come from outer products of vectors 
    // // This class describes an arbitrary rank tensor created iteratively 
    // class lorentz_tensor 
    // {
    //     public: 

    //     // Accessor function for elements 
    //     complex operator[](std::vector<lorentz_index> vals);

    //     // Get the rank of the tensor (number of open indicies)
    //     inline int rank(){ return _rank; };

    //     // Re-assignment
    //     lorentz_tensor & operator=(lorentz_tensor const & T);

    //     // Multiply and divide by a constant
    //     lorentz_tensor & operator*=(complex c);
    //     lorentz_tensor & operator/=(complex c);
    //     lorentz_tensor & operator+=(lorentz_tensor const & T);

    //     private:

    //     // Check if an outside tensor is compatable to be added with this one
    //     bool is_compatible(lorentz_tensor T)
    //     {
    //         return (T.rank() == rank());
    //     };

    //     lorentz_tensor(std::vector<lorentz_vector> list)
    //     : _rank(list.size())
    //     {
    //         for (auto index : list)
    //         {
    //             _indices.push_back(index);
    //         };

    //         if (_indices.size() != _rank) fatal("lorentz_tensor", "Failed to initialize with correct size!");
    //     };

    //     // Make the outer product class a friend which can access the _indices
    //     // This makes it clear that although lorentz_tensor is  generic class
    //     // Currently its implementation only fully allows those tensors decomposable as
    //     // tensor products of spin-1 subspaces
    //     friend lorentz_tensor outer(std::vector<lorentz_vector>);
    //     friend lorentz_tensor outer(std::vector<lorentz_tensor>);

    //     int _rank = 0;
    //     std::vector<lorentz_vector> _indices;
    // };

    // // Outer product of any number of lorentz_vectors produces a tensor
    // inline lorentz_tensor outer(std::vector<lorentz_vector> xs)
    // {   
    //     return lorentz_tensor(xs); 
    // };

    // // Multiple tensors may also be combined in order they are listed from left to right
    // inline lorentz_tensor outer(std::vector<lorentz_tensor> list)
    // {
    //     std::vector<lorentz_vector> aggregate = {};
    //     for (auto t : list)
    //     {
    //         aggregate.insert( aggregate.end(), std::make_move_iterator(t._indices.begin()), std::make_move_iterator(t._indices.end()));
    //     };
    //     return lorentz_tensor(aggregate);
    // };
};

#endif 