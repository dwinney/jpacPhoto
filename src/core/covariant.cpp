// Container class for covariant quantities!
// Since not all amplitudes need these they're seperated for convience
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "covariant.hpp"

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    // Methods and properties of lorentz_vecter
    
    // Access the mu-th element
    complex lorentz_vector::operator[](lorentz_index mu)
    {
        return _components[static_cast<int>(mu)]; 
    };

    // Negate a 4-vector
    lorentz_vector lorentz_vector::operator-()
    {
        lorentz_vector result;
        for (auto MU : LORENTZ_INDICES)
        {
            int mu = static_cast<int>(MU);
            result._components[mu] = - _components[mu];
        };
        return result;
    };

    // Add two four-vectors together
    lorentz_vector lorentz_vector::operator+(lorentz_vector const & p)
    {
        lorentz_vector result;
        for (auto MU : LORENTZ_INDICES)
        {
            int mu = static_cast<int>(MU);
            result._components[mu] = _components[mu] + p._components[mu];
        };
        return result;
    };

    // Define new vectors from old ones
    lorentz_vector lorentz_vector::operator=(lorentz_vector const & p)
    {
        return lorentz_vector(p._components);
    };

    // Multiply a vector by a scalar
    lorentz_vector operator*(complex c, lorentz_vector p)
    {
        std::array<complex,4> cp;
        for (auto mu : LORENTZ_INDICES)
        {
            cp[static_cast<int>(mu)] = c * p[mu];
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
    // Methods related to our caching to minimize calculations

    // Check if energy step requires recalculating quantities
    void covariant::cache::check(double s, double t)
    {
        bool s_changed      =   unsynced(s, _s);
        bool masses_changed =   unsynced(_kinematics->get_beam_mass(), _mB)
                             || unsynced(_kinematics->get_target_mass(), _mT)
                             || unsynced(_kinematics->get_meson_mass(), _mX)
                             || unsynced(_kinematics->get_recoil_mass(), _mR);

        if ( s_changed || masses_changed ) recalculate(s);

        // Recalculating theta dependence is super easy 
        // so this doesnt need a seperate method
        double theta = _kinematics->theta_s(s, t);
        if ( are_equal(theta, _theta) )
        {
            _cos = cos(theta);
            _sin = sin(theta);
            _theta = theta;
        };

        return;
    };

    // The complicated recalculation only comes with changing energy dependence
    void covariant::cache::recalculate(double s)
    {
        // save the latest masses
        _mB = _kinematics->get_beam_mass();
        _mT = _kinematics->get_target_mass();
        _mX = _kinematics->get_meson_mass();
        _mR = _kinematics->get_recoil_mass();
        _s  = s;

        // then update momenta
        _qi = _kinematics->initial_momentum(s);
        _qf = _kinematics->final_momentum(s);

        // and energies
        _EB = _kinematics->beam_energy(s);
        _ET = _kinematics->target_energy(s);
        _EX = _kinematics->meson_energy(s);
        _ER = _kinematics->recoil_energy(s);
    };

    // ---------------------------------------------------------------------------
    // FOUR MOMENTA 

    // Beam momentum
    lorentz_vector covariant::q()
    {   
        // Aligned with the positive z-hat axis
        return lorentz_vector( {_cache._EB, 0., 0., _cache._qi} );
    };

    // Beam momentum
    lorentz_vector covariant::p()
    {   
        // Aligned with the positive z-hat axis
        return lorentz_vector( {_cache._ET, 0., 0., -_cache._qi} );
    };

    // Produced meson momentum
    lorentz_vector covariant::q_prime()
    {
        // In x-z plane with axis angle _theta from the +z
        return lorentz_vector( {_cache._EX, _cache._qf*_cache._sin, 0., _cache._qf*_cache._cos} );
    };

    // Produced meson momentum
    lorentz_vector covariant::p_prime()
    {
        // In x-z plane with axis angle _theta from the +z
        return lorentz_vector( {_cache._ER, -_cache._qf*_cache._sin, 0., -_cache._qf*_cache._cos} );
    };

};
