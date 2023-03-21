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

    // -----------------------------------------------------------------------
    // In most useful contexts, these come from outer products of vectors 
    // This class describes an arbitrary rank tensor created iteratively 

    template<class Type>
    class tensor_object 
    {
        public: 
        virtual int rank() = 0;
        virtual Type operator()(std::vector<lorentz_index> indices) = 0;
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
        
        // Constructor for tensors by explicit element specification
        // This only wors for size 1 and 2, i.e. lorentz vectors and metric
        // Tensors always initialized with normalization of 1
        // Higher rank tensors should be created by outer-products
        lorentz_tensor<Type,Rank>(tensor_entries entries)
        : _entries(entries), _is_sum(false)
        {};


        // Copy constructor
        lorentz_tensor<Type,Rank>(lorentz_tensor<Type,Rank> const & old)
        : _entries(old._entries), _subtensors(old._subtensors),
          _lhsN(old._lhsN), _rhsN(old._rhsN), _conj(old._conj)
        {};

        // Conversions between types

        operator lorentz_tensor<dirac_matrix,Rank>()
        {
            return identity<dirac_matrix>() * *this;
        };

        operator lorentz_tensor<dirac_spinor,Rank>()
        {
            return identity<dirac_spinor>() * *this;
        };

        // Vector accessor (rank 1)
        inline Type operator()(lorentz_index mu)
        {
            if (Rank != 1) return error("lorentz_tensor[]", "Incorrect number of indices passed!", NaN<Type>());
            
            Type result = _lhsN * _entries[0][+mu] * _rhsN;
            return (_conj) ? conj(result) : result; 
        };

        // Matrix accessor (rank 2)
        inline Type operator()(lorentz_index mu, lorentz_index nu)
        {
            if (Rank != 2) return error("lorentz_tensor[]", "Incorrect number of indices passed!", NaN<Type>());
            
            Type result = _lhsN * _entries[+mu][+nu] * _rhsN;
            return (_conj) ? conj(result) : result; 
        };

        // Evaluate subtensors iteratively
        inline Type operator()(std::vector<lorentz_index> indices)
        {
            if (indices.size() != Rank) return error("lorentz_tensor[]", "Incorrect number of indices passed!", NaN<Type>());
            if (Rank == 1) return operator()(indices[0]);
            if (Rank == 2) return operator()(indices[0], indices[1]);

            if (_is_sum)
            {
                Type sum = zero<Type>();
                for (auto T : _to_sum)
                {
                    sum += T(indices);
                }
                return _lhsN*sum*_rhsN;
            };

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
     
            return _lhsN*prod*_rhsN;
        };

        // Get the rank of the tensor (number of open indicies)
        inline int rank(){ return _rank; };

        // Re-assignment
        inline lorentz_tensor<Type,Rank> & operator=(lorentz_tensor<Type,Rank> const & T)
        {
            _entries        = T._entires;
            _lhsN           = T._lhsN;
            _rhsN           = T._rhsN;
            _conj           = T._conj;
            return *this;
        };

        // Multiply and divide by a constant
        inline lorentz_tensor<Type,Rank> & operator*=(Type c)
        {
            // Apply the constant to the first index simply
            _rhsN *= c;
            return *this;
        };

        inline lorentz_tensor<Type,Rank> & operator/=(complex c)
        {
            _lhsN *= (1./c);
            return *this;
        };

        // Negate a tensor
        inline lorentz_tensor<Type,Rank> operator-()
        {
            lorentz_tensor<Type,Rank> neg = *this;
            neg *= -1;
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

        // Implicit constructor, stores pointers to constituent tensors of smaller rank
        lorentz_tensor<Type,Rank>(std::vector<std::shared_ptr<tensor_object<Type>>> Ts)
        : _subtensors(Ts), _is_sum(false)
        {};

        // Sum constructor for large matrices
        lorentz_tensor<Type,Rank>(std::vector<lorentz_tensor<Type,Rank>> Ts)
        : _to_sum(Ts), _is_sum(true)
        {};

        // Rest of these classes exist outside the class but form the core ways to interact
        // with tensors

        // Sum vectors together
        template<class T, int R>
        friend lorentz_tensor<T,R> operator+(lorentz_tensor<T,R> lhs, lorentz_tensor<T,R> rhs);

        template<class T, int R>
        friend lorentz_tensor<T,R> operator-(lorentz_tensor<T,R> lhs, lorentz_tensor<T,R> rhs);

        // Tensor products 

        // Two vectors -> matrix
        template<class T>
        friend lorentz_tensor<T,2> tensor_product(lorentz_tensor<T,1>, lorentz_tensor<T,1>);

        // R1-tensor + R2-tensor -> (R1+R2)-tensor
        template<class T, int R1, int R2>
        friend lorentz_tensor<T,R1+R2> tensor_product(lorentz_tensor<T,R1>, lorentz_tensor<T,R2>);

        // N vectors -> rank-N tensor
        template<class T, int R>
        friend lorentz_tensor<T,R> tensor_product(std::array<lorentz_tensor<T,1>,R>);

        // "Capture" constructors
        
        template<class LType, int R>
        friend lorentz_tensor<LType,R> operator*(LType c, lorentz_tensor<complex,R> rhs);

        // IF we end up with a rank-0 object, return the scalar piece
        template<class T>
        friend T flatten(lorentz_tensor<T,0>);

        const static int _rank = Rank;
        bool _conj = false;
        
        // If this tensor represents a single rank 1 or 2, we store its explicit entries
        tensor_entries _entries;

        // If its a more complicated object (rank >= 3) we store pointers to its component tensors
        std::vector<tensor_ptr> _subtensors;

        // If the tensor is a metric, it may still be multiplied by a constant.
        // This however may be generalized depending on the substructure inside the tensor
        Type _lhsN = identity<Type>(), _rhsN = identity<Type>();

        // If not, we the tensor serves as a store of constituent other tensors 
        // which are summed element wise instead
        bool _is_sum = false;
        std::vector<lorentz_tensor<Type,Rank>> _to_sum;

        // Both type and rank has to be the same to add
        void add_tensor(lorentz_tensor<Type,Rank> T)
        {
            if (!_is_sum)
            {
                warning("add_tensor()", "Cannot add tensor to pre-initialized one. Initilize a new tensor as the sum!");
                return;
            };
            _to_sum.push_back(T); 
        };
    };

    // ---------------------------------------------------------------------------
    // NON-MEMBER FUNCTIONS

    // ---------------------------------------------------------------------------
    // "Constructor" functions 

    // Outer product of any number two lorentz_vectors produces a rank-2
    template<class Type>
    inline lorentz_tensor<Type,2> tensor_product(lorentz_tensor<Type,1> p, lorentz_tensor<Type,1> q)
    {   
        Type z = zero<Type>(); 
        std::vector<std::array<Type,4>> pq;
        for (auto mu : LORENTZ_INDICES)
        {  
            std::array<Type,4> row = {z, z, z, z};
            for (auto nu : LORENTZ_INDICES)
            {
                row[+nu] = p(mu) * q(nu);
            };
            pq.push_back(row);
        };
        return lorentz_tensor<Type,2>(pq); 
    };

    // N vectors -> rank-N tensor
    template<class T, int R>
    inline lorentz_tensor<T,R> tensor_product(std::array<lorentz_tensor<T,1>,R> vs)
    {
        std::vector< std::shared_ptr< tensor_object<T> > > ptrs;
        for (auto vec : vs)
        {
            ptrs.push_back( std::make_shared<lorentz_tensor<T,1>>(vec) );
        };

        return lorentz_tensor<T,R>(ptrs);
    };

    // Generic tensor products for to create abritrary sized tensors
    template<class Type, int R1, int R2> 
    inline lorentz_tensor<Type,R1+R2> tensor_product(lorentz_tensor<Type,R1> lhs, lorentz_tensor<Type,R2> rhs)
    {
        std::vector< std::shared_ptr< tensor_object<Type> > > ptrs;
        ptrs.push_back( std::make_shared<lorentz_tensor<Type,R1>>(lhs) );
        ptrs.push_back( std::make_shared<lorentz_tensor<Type,R2>>(rhs) );

        return lorentz_tensor<Type, R1+R2>(ptrs);
    };

    // "Constructor" of an arbitrary rank-1 tensor (i.e. vector)
    template<class Type>
    inline lorentz_tensor<Type,1> lorentz_vector(std::array<Type,4> entries)
    {
        return lorentz_tensor<Type,1>({entries});
    };

    // The metric tensor is unique rank-2 tensor
    //  as it cannot be decomposed into a tensor product of rank-1 vectors
    inline lorentz_tensor<complex,2> metric_tensor()
    {
        std::vector<std::array<complex,4>> gmunu = {{1,  0,  0,  0},
                                                    {0, -1,  0,  0},
                                                    {0,  0, -1,  0},
                                                    {0,  0,  0, -1}};
        return lorentz_tensor<complex,2>(gmunu);
    };

    // ---------------------------------------------------------------------------
    // Scalar tensors "capture" dirac data types they get multipled by a non-scalar dirac type

    template<class LType>
    inline lorentz_tensor<LType,1> operator*(LType c, lorentz_tensor<complex,1> rhs)
    {
        LType id = identity<LType>();
        std::array<LType,4> new_entries = {id, id, id, id};
        for (auto mu : LORENTZ_INDICES)
        {  
            complex old_entry = rhs(mu);
            new_entries[+mu]  = c * old_entry; 
        };

        return lorentz_vector<LType>(new_entries);
    };

    template<class LType>
    inline lorentz_tensor<LType,2> operator*(LType c, lorentz_tensor<complex,2> rhs)
    {
        LType id = identity<LType>();
        std::vector<std::array<LType,4>> new_entries = {};
        for (auto mu : LORENTZ_INDICES)
        {  
            std::array<LType,4> new_column = {id, id, id, id};
            for (auto nu : LORENTZ_INDICES)
            {
                complex old_entry = rhs(mu, nu);
                new_column[+nu]   = c * old_entry; 
            };
            new_entries.push_back((new_column));
        };        
        return lorentz_tensor<LType,2>(new_entries);
    };

    template<class RType>
    inline lorentz_tensor<RType,1> operator*(lorentz_tensor<complex,1> lhs, RType c)
    {
        RType id = identity<RType>();
        std::array<RType,4> new_entries = {id, id, id, id};
        for (auto mu : LORENTZ_INDICES)
        {  
            complex old_entry = lhs(mu);
            new_entries[+mu]  = old_entry * c; 
        };

        return lorentz_vector<RType>(new_entries);
    };

    template<class RType>
    inline lorentz_tensor<RType,2> operator*(lorentz_tensor<complex,2> lhs, RType c)
    {
        RType id = identity<RType>();
        std::vector<std::array<RType,4>> new_entries = {};
        for (auto mu : LORENTZ_INDICES)
        {  
            std::array<RType,4> new_column = {id, id, id, id};
            for (auto nu : LORENTZ_INDICES)
            {
                complex old_entry = lhs(mu, nu);
                new_column[+nu]   = old_entry * c; 
            };
            new_entries.push_back((new_column));
        };        
        return lorentz_tensor<RType,2>(new_entries);
    };

    // ---------------------------------------------------------------------------
    // Interactions with constants 
    // These should be handled seperately because or else int * lorentz_tensor<complex> will 
    // return a lorentz_tensor<int> which we dont want. 

    template<class LType, int R>
    inline lorentz_tensor<LType,R> operator/(lorentz_tensor<LType,R> lhs, complex c)
    {
        return (1./c) * lhs;
    };

    inline lorentz_tensor<complex,1> operator*(double c, lorentz_tensor<complex,1> rhs)
    {
        return complex(c) * rhs;
    };

    inline lorentz_tensor<complex,2> operator*(double c, lorentz_tensor<complex,2> rhs)
    {
        return complex(c) * rhs;
    };

    inline lorentz_tensor<complex,1> operator*(int c, lorentz_tensor<complex,1> rhs)
    {
        return complex(c) * rhs;
    };

    inline lorentz_tensor<complex,2> operator*(int c, lorentz_tensor<complex,2> rhs)
    {
        return complex(c) * rhs;
    };

    inline lorentz_tensor<complex,1> operator*(lorentz_tensor<complex,1> lhs, double c)
    {
        return complex(c) * lhs;
    };

    inline lorentz_tensor<complex,2> operator*(lorentz_tensor<complex,2> lhs, double c)
    {
        return complex(c) * lhs;
    };

        inline lorentz_tensor<complex,1> operator*(lorentz_tensor<complex,1> lhs, int c)
    {
        return complex(c) * lhs;
    };

    inline lorentz_tensor<complex,2> operator*(lorentz_tensor<complex,2> lhs, int c)
    {
        return complex(c) * lhs;
    };

    // ---------------------------------------------------------------------------
    // Operations between tensors

    // Add two vectors together
    template<class T, int R>
    inline lorentz_tensor<T,R> operator+(lorentz_tensor<T,R> lhs, lorentz_tensor<T,R> rhs)
    {
        lorentz_tensor<T,R> sum({lhs, rhs});
        return sum;
    };

    // Add two vectors together
    template<class T, int R>
    inline lorentz_tensor<T,R> operator-(lorentz_tensor<T,R> lhs, lorentz_tensor<T,R> rhs)
    {
        return lhs + (-rhs);
    };

    // Add two vectors together
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

    // Add matrices together
    template<class T>
    inline lorentz_tensor<T,2> operator+(lorentz_tensor<T,2> lhs, lorentz_tensor<T,2> rhs)
    {
        std::vector<std::array<T,4>> sum;
        for (auto mu : LORENTZ_INDICES)
        {  
            std::array<T,4> row;
            for (auto nu : LORENTZ_INDICES)
            {
                row[+nu] = lhs(mu, nu) + rhs(mu, nu);
            };
            sum.push_back(row);
        };
        return lorentz_tensor<T,2>(sum);
    };
};

#endif