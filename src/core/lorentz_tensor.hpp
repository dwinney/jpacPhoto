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

#ifndef LORENTZ_TENSOR_HPP
#define LORENTZ_TENSOR_HPP

#include <cstddef>
#include <vector>
#include <memory>

#include "constants.hpp"
#include "dirac_spinor.hpp"
#include "dirac_matrix.hpp"

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    // Lorentz index handling

    // Enclosing lorentz indices in an enum makes summing over them more transparent
    // Also has the benefit of not needing to check type or scope of a given index
    enum class lorentz_index: int {t = 0, x = 1, y = 2, z = 3};
    const std::array<lorentz_index,4> LORENTZ_INDICES = {lorentz_index::t, lorentz_index::x, lorentz_index::y, lorentz_index::z};
    
    // e.g. int x = +lorentz_index(mu);
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
    // This function is only to make the "capture" constructors below to work

    template<class In, class Out>
    inline auto capture(In x)
    {
        return identity<Out>() * x;
    };

    // The above is the behavior we want in most cases except the
    // following two exceptions which are errors
    template<>
    inline auto capture<dirac_matrix,dirac_spinor>(dirac_matrix x)
    {
        return NaN<dirac_spinor>();
    };
    template<>
    inline auto capture<dirac_spinor,dirac_matrix>(dirac_spinor x)
    {
        return NaN<dirac_matrix>();
    };


    // -----------------------------------------------------------------------
    // In most useful contexts, these come from outer products of vectors 
    // This class describes an arbitrary rank tensor created iteratively 

    template<class Type, int Rank>
    class lorentz_tensor;

    template<class Type = complex>
    class tensor_object 
    {
        public: 
        virtual const int rank() = 0;
        virtual Type operator()(std::vector<lorentz_index> indices) = 0;

        virtual std::shared_ptr<tensor_object<dirac_matrix>> matrixify(){ return nullptr; };
        virtual std::shared_ptr<tensor_object<dirac_spinor>> spinorify(){ return nullptr; };
    };

    template<class Type, int Rank>
    class lorentz_tensor : public tensor_object<Type>
    {
        // The shape is always a Rank x 4 tensor. 
        // The Type depends on the possibility of dirac structures
        // Scalar (complex/double) quantities is any object with only Lorentz structure
        // dirac_spinor and dirac_matrix tensors allow constructor of high half-integer spin wavefunctions
        using tensor_entries = std::vector<std::array<Type,4>>; 
        using tensor_ptr     = std::shared_ptr<tensor_object<Type>>;
        
        public: 

        // Default tensor with nothing initialized
        lorentz_tensor<Type,Rank>(){};

        // Copy constructor
        lorentz_tensor<Type,Rank>(lorentz_tensor<Type,Rank> const & old)
        : _entries(old._entries), _metric(old._metric),
          _lhsN(old._lhsN), _rhsN(old._rhsN), _conj(old._conj), 
          _subtensors(old._subtensors),
          _is_sum(old._is_sum)
        {};
        
        // Constructor for tensors by explicit element specification
        // This only wors for size 1 and 2, i.e. lorentz vectors and metric
        // Tensors always initialized with normalization of 1
        // Higher rank tensors should be created by outer-products
        lorentz_tensor<Type,Rank>(tensor_entries entries)
        : _entries(entries), _is_sum(false), _metric(Rank == 2)
        {};

        // Implicit constructor, stores pointers to constituent tensors of smaller rank
        lorentz_tensor<Type,Rank>(std::vector<tensor_ptr> Ts, bool sum)
        : _subtensors(Ts), _is_sum(sum)
        {};

        // Get the rank of the tensor (number of open indicies)
        const inline int rank(){ return Rank; };

        // Re-assignment, just makes a copy basically
        inline lorentz_tensor<Type,Rank> & operator=(lorentz_tensor<Type,Rank> const & T)
        {
            _entries        = T._entries;
            _lhsN           = T._lhsN;
            _rhsN           = T._rhsN;
            _conj           = T._conj;
            _subtensors     = T._subtensors;
            _is_sum         = T._is_sum;
            _metric         = T._metric;
            return *this;
        };

        // -----------------------------------------------------------------------
        // Evaluate elements of the tensor

        inline Type operator()(std::vector<lorentz_index> indices)
        {
            if (indices.size() != Rank) return error("lorentz_tensor", "Incorrect number of indices passed!", NaN<Type>());
            if (_is_sum)
            {
                Type sum = zero<Type>();
                for (auto T : _subtensors)
                {
                    sum += T->operator()(indices);
                };

                Type result = _lhsN*sum*_rhsN;
                return (_conj) ? conj(result) : result;
            };
            
            if (Rank == 1)
            {
                Type result = _lhsN*_entries[0][+indices[0]]*_rhsN;
                return (_conj) ? conj(result) : result;
            } 
            if (_metric)
            {
                Type result = _lhsN*_entries[+indices[0]][+indices[1]]*_rhsN;
                return (_conj) ? conj(result) : result;
            } 

            auto running_index = indices.begin();
            int subrank0 = _subtensors[0]->rank();
            std::vector<lorentz_index> ind0 = std::vector<lorentz_index>(running_index, running_index + subrank0);

            Type prod = _subtensors[0]->operator()(ind0);
            running_index += subrank0;

            for (int i = 1; i < _subtensors.size(); i++)
            {
                int subrank = _subtensors[i]->rank();
                prod *= _subtensors[i]->operator()( std::vector<lorentz_index>(running_index, running_index + subrank));
                running_index += subrank;
            };
            Type result = _lhsN*prod*_rhsN;
            return (_conj) ? conj(result) : result; 
        };
        inline Type operator()(lorentz_index mu){ return operator()({{mu}}); };
        inline Type operator()(lorentz_index mu, lorentz_index nu){ return operator()({{mu,nu}});};

        // -----------------------------------------------------------------------
        // Unary operations 

        // Multiply and divide by a constant
        inline lorentz_tensor<Type,Rank> & operator*=(Type c)
        {
            // Apply the constant to the first index simply
            _rhsN *= c*identity<Type>();
            return *this;
        };
        inline lorentz_tensor<Type,Rank> & operator/=(complex c)
        {
            _rhsN *= (1/c)*identity<Type>();
            return *this;
        };

        // Negate a tensor
        inline lorentz_tensor<Type,Rank> operator-()
        {
            lorentz_tensor<Type,Rank> neg = *this;
            neg._lhsN *= -1*identity<Type>();
            return neg;
        };

        // Return a copy of itself with the conj flag flipped
        inline lorentz_tensor<Type,Rank> conjugate()
        {
            lorentz_tensor<Type,Rank> copy = *this;
            copy._conj = !_conj; 
            return copy;
        };

        private:

        // Multiply tensors by non-tensor objects
        template<class T, int R>
        friend lorentz_tensor<T,R> operator*(T c, lorentz_tensor<T,R> rhs);
        template<class T, int R>
        friend lorentz_tensor<T,R> operator*(lorentz_tensor<T,R> lhs, T c);

        // "Capture" constructor require more access than other contructors
        template<int R>
        friend lorentz_tensor<dirac_matrix,R> operator*(dirac_matrix c, lorentz_tensor<complex,R> rhs);
        template<int R>
        friend lorentz_tensor<dirac_spinor,R> operator*(dirac_spinor c, lorentz_tensor<complex,R> rhs);

        // Tensor dirac_spinor * matrices
        template<int R>
        friend lorentz_tensor<dirac_spinor,R> operator*(dirac_matrix M, lorentz_tensor<dirac_spinor,R> rhs);
        template<int R>
        friend lorentz_tensor<dirac_spinor,R> operator*(lorentz_tensor<dirac_spinor,R> lhs, dirac_matrix M);

        // Whether or not to take the conjugate of the element
        bool _conj = false;
        
        // Whether this tensor started life as a metric tensor
        bool _metric = false;

        // If not, we the tensor serves as a store of constituent other tensors 
        // which are summed element wise instead
        bool _is_sum = false;

        // If this tensor represents a single rank 1 or 2, we store its explicit entries
        tensor_entries _entries;

        // If its a more complicated object (rank >= 3) we store pointers to its component tensors
        std::vector<tensor_ptr> _subtensors;

        // If the tensor is a metric, it may still be multiplied by a constant.
        // This however may be generalized depending on the substructure inside the tensor
        Type _lhsN = identity<Type>(), _rhsN = identity<Type>();

        // Both type and rank has to be the same to add
        inline void add_tensor(lorentz_tensor<Type,Rank> T)
        {
            if (!_is_sum)
            {
                warning("add_tensor()", "Cannot add tensor to pre-initialized one. Initilize a new tensor as the sum!");
                return;
            };
            _subtensors.push_back(std::make_shared<lorentz_tensor<Type,Rank>>(T)); 
        };

        // Return a dirac_spinor version of complex tensor
        inline std::shared_ptr<tensor_object<dirac_spinor>> spinorify()
        {
            dirac_spinor id = identity<dirac_spinor>();

            // If its one of the special cases (Rank 1 or a metric) go elementwise
            if (Rank == 1 || _metric)
            {
                std::vector<std::array<dirac_spinor,4>> new_entries;
                for (auto index : _entries)
                {
                    std::array<dirac_spinor,4> row;
                    for (auto mu : LORENTZ_INDICES)
                    {
                        row[+mu] = capture<Type,dirac_spinor>(_lhsN*index[+mu]*_rhsN);
                    };
                    new_entries.push_back(row);
                }
                return std::make_shared<lorentz_tensor<dirac_spinor,Rank>>(new_entries);
            };

            // Else go iteratively through sum terms first
            if (_is_sum)
            {
                std::vector<std::shared_ptr<tensor_object<dirac_spinor>>> new_to_sum;
                for (auto term : _subtensors)
                {
                    new_to_sum.push_back(term->spinorify());
                }
                // True in constructor means the subtesnors are summed not multipled
                return std::make_shared<lorentz_tensor<dirac_spinor,Rank>>(new_to_sum, true);
            };

            // If we reach this point it means we have a tensor product on our hands,
            // similar to above we simply go through each subtensor and spinorify it in turn
            std::vector<std::shared_ptr<tensor_object<dirac_spinor>>> new_subtensors;
            for (auto term : _subtensors)
            {
                new_subtensors.push_back(term->spinorify());
            }

            // Unlike above, the false means multipled
            return std::make_shared<lorentz_tensor<dirac_spinor,Rank>>(new_subtensors, false);
        };

        // Return a dirac_matrix
        // This creates a problem if mixing dirac_matrices and spinors 
        inline std::shared_ptr<tensor_object<dirac_matrix>> matrixify()
        {    
            // If its one of the special cases (Rank 1 or a metric) go elementwise
            if (Rank == 1 || _metric)
            {
                std::vector<std::array<dirac_matrix,4>> new_entries;
                for (auto index : _entries)
                {
                    std::array<dirac_matrix,4> row;
                    for (auto mu : LORENTZ_INDICES)
                    {
                        row[+mu] = capture<Type, dirac_matrix>(_lhsN*index[+mu]*_rhsN);
                    };
                    new_entries.push_back(row);
                }
                return std::make_shared<lorentz_tensor<dirac_matrix,Rank>>(new_entries);
            };

            // Else go iteratively through sum terms first
            if (_is_sum)
            {
                std::vector<std::shared_ptr<tensor_object<dirac_matrix>>> new_to_sum;
                for (auto term : _subtensors)
                {
                    new_to_sum.push_back(term->matrixify());
                }
                // True in constructor means the subtesnors are summed not multipled
                return std::make_shared<lorentz_tensor<dirac_matrix,Rank>>(new_to_sum, true);
            };

            // If we reach this point it means we have a tensor product on our hands,
            // similar to above we simply go through each subtensor and matrixify it in turn
            std::vector<std::shared_ptr<tensor_object<dirac_matrix>>> new_subtensors;
            for (auto term : _subtensors)
            {
                new_subtensors.push_back(term->matrixify());
            }

            // Unlike above, the false means multipled
            return std::make_shared<lorentz_tensor<dirac_matrix,Rank>>(new_subtensors, false);
        };
    };

    // ---------------------------------------------------------------------------
    // NON-MEMBER FUNCTIONS

    // ---------------------------------------------------------------------------
    // "Constructor" functions 

    // N vectors -> rank-N tensor
    template<class T = complex, int R>
    inline lorentz_tensor<T,R> tensor_product(std::array<lorentz_tensor<T,1>,R> vs)
    {
        std::vector< std::shared_ptr< tensor_object<T> > > ptrs;
        for (auto vec : vs)
        {
            ptrs.push_back(std::make_shared<lorentz_tensor<T,1>>(vec));
        };

        return lorentz_tensor<T,R>(ptrs, false);
    };

    // Generic tensor products for to create abritrary sized tensors
    template<class Type = complex, int R1, int R2> 
    inline lorentz_tensor<Type,R1+R2> tensor_product(lorentz_tensor<Type,R1> lhs, lorentz_tensor<Type,R2> rhs)
    {
        std::vector< std::shared_ptr< tensor_object<Type> > > ptrs;
        ptrs.push_back( std::make_shared<lorentz_tensor<Type,R1>>(lhs) );
        ptrs.push_back( std::make_shared<lorentz_tensor<Type,R2>>(rhs) );

        return lorentz_tensor<Type, R1+R2>(ptrs, false);
    };

    // "Constructor" of an arbitrary rank-1 tensor (i.e. vector)
    template<class Type = complex>
    inline lorentz_tensor<Type,1> lorentz_vector(std::array<Type,4> entries)
    {
        return lorentz_tensor<Type,1>({entries});
    };

    // The metric tensor is unique rank-2 tensor
    //  as it cannot be decomposed into a tensor product of rank-1 vectors
    template<class Type = complex>
    inline lorentz_tensor<Type,2> metric_tensor()
    {
        Type z = zero<Type>(), id = identity<Type>();
        std::vector<std::array<Type,4>> gmunu = {{id,  z,  z,  z},
                                                 {z, -id,  z,  z},
                                                 {z,  z, -id,  z},
                                                 {z,  z,  z, -id}};
        return lorentz_tensor<Type,2>(gmunu);
    };

    // ---------------------------------------------------------------------------
    // Operations between tensors

    // Add two vectors together.
    // For arbitrary size, this is done by saving tensors and adding iteratively through terms in the sum
    // when an element is asked for
    template<class T, int R>
    inline lorentz_tensor<T,R> operator+(lorentz_tensor<T,R> lhs, lorentz_tensor<T,R> rhs)
    {
        auto l = std::make_shared<lorentz_tensor<T,R>>(lhs);
        auto r = std::make_shared<lorentz_tensor<T,R>>(rhs);
        lorentz_tensor<T,R> sum({l,r}, true);
        return sum;
    };

    // Only when adding two vectors (rank-1) do we save elements explicitly
    template<class T>
    inline lorentz_tensor<T,1> operator+(lorentz_tensor<T,1> lhs, lorentz_tensor<T,1> rhs)
    {
        std::array<T,4> sum;
        for (auto mu : LORENTZ_INDICES)
        {
            sum[+mu] = lhs(mu) + rhs(mu);
        };
        return lorentz_vector(sum);
    };

    // Subtract two vectors together
    template<class T, int R>
    inline lorentz_tensor<T,R> operator-(lorentz_tensor<T,R> lhs, lorentz_tensor<T,R> rhs)
    {
        return lhs + (-rhs);
    };

    // ---------------------------------------------------------------------------
    // Scalar tensors "capture" dirac data types they get multipled by a non-scalar dirac type
    // We dont use template here because compiler gets confused between complex, double, int, etc.
    // We want to always default to complex even if we write p = -2*q;

    // complex * dirac_matrix -> dirac_matrix
    template<int R>
    inline lorentz_tensor<dirac_matrix,R> operator*(dirac_matrix c, lorentz_tensor<complex,R> rhs)
    {
        auto matrified = rhs.matrixify();
        return c * static_cast<lorentz_tensor<dirac_matrix,R>&>(*matrified);
    };
    template<int R>
    inline lorentz_tensor<dirac_matrix,R> operator*(lorentz_tensor<complex,R> lhs, dirac_matrix c)
    {
        return c * lhs;
    };

    // complex * dirac_spinor -> dirac_spinor
    template<int R>
    inline lorentz_tensor<dirac_spinor,R> operator*(dirac_spinor c, lorentz_tensor<complex,R> rhs)
    {
        auto matrified = rhs.spinorify();
        return c * static_cast<lorentz_tensor<dirac_spinor,R>&>(*matrified);
    };
    template<int R>
    inline lorentz_tensor<dirac_spinor,R> operator*(lorentz_tensor<complex,R> lhs, dirac_spinor c)
    {
        return c * lhs;
    };
    
    // --------------------------------------------------------------------------
    // Interactions with constant scalar types 

    // Multiplication by the same type
    // Order matters here only becasue of possible matrix multipication
    template<class T = complex, int R>
    inline lorentz_tensor<T,R> operator*(T c, lorentz_tensor<T,R> rhs)
    {
        rhs._lhsN *= c;
        return rhs;
    };
    template<class T = complex, int R>
    inline lorentz_tensor<T,R> operator*(lorentz_tensor<T,R> lhs, T c)
    {
        lhs._rhsN *= c;
        return lhs;
    };

    // Multiplication by an integer type, order no longer matters

    // int
    template<class T = complex, int R>
    inline lorentz_tensor<T,R> operator*(int c, lorentz_tensor<T,R> rhs)
    {
        return (c*identity<T>())*rhs;
    };
    template<class T = complex, int R>
    inline lorentz_tensor<T,R> operator*(lorentz_tensor<T,R> lhs, int c){ return c*lhs; };
    
    // we need to define complex * dirac_spinor and dirac_matrix seperately becasue else the template above gets confused
    template<int R>
    inline lorentz_tensor<dirac_spinor,R> operator*(complex c, lorentz_tensor<dirac_spinor,R> rhs){ return (c*identity<dirac_spinor>())*rhs; };
    template<int R>
    inline lorentz_tensor<dirac_spinor,R> operator*(lorentz_tensor<dirac_spinor,R> lhs, complex c){ return (c*identity<dirac_spinor>())*lhs; }; 

    template<int R>
    inline lorentz_tensor<dirac_matrix,R> operator*(complex c, lorentz_tensor<dirac_matrix,R> rhs){ return (c*identity<dirac_matrix>())*rhs; };
    template<int R>
    inline lorentz_tensor<dirac_matrix,R> operator*(lorentz_tensor<dirac_matrix,R> lhs, complex c){ return (c*identity<dirac_matrix>())*lhs; }; 

    // Division only makes sense for dirac scalar types
    template<class LType, class RType, int R>
    inline lorentz_tensor<LType,R> operator/(lorentz_tensor<LType,R> lhs, RType c){ return (1./complex(c)) * lhs; };

    
    // --------------------------------------------------------------------------
    // Interactions with constant non-scalar types 

    // Multiplying matrices and spinors with tensor structure
    template<int R>
    inline lorentz_tensor<dirac_spinor,R> operator*(dirac_matrix M, lorentz_tensor<dirac_spinor,R> rhs)
    {
        rhs._lhsN = M*rhs._lhsN;
        return rhs;
    };
    template<int R>
    inline lorentz_tensor<dirac_spinor,R> operator*(lorentz_tensor<dirac_spinor,R> lhs, dirac_matrix M)
    {
        lhs._rhsN = lhs._rhsN*M;
        return lhs;
    };
};

#endif