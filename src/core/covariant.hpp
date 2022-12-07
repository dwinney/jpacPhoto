// Container class for covariant quantities!
// Since not all amplitudes need these they're seperated for convience
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef COVARIANT_HPP
#define COVARIANT_HPP

#include "constants.hpp"
#include "kinematics.hpp"
// #include "dirac_spinor.hpp"
// #include "rarita_spinor.hpp"
// #include "two_body_state.hpp"
// #include "polarization_vector.hpp"
// #include "gamma_matrices.hpp"
#include <functional>
#include <initializer_list>

namespace jpacPhoto
{   
    // We dont actually need most properties of lorentz vectors since we will always work in the s-channel CM frame
    // we just need a basic vector class to differentiate it from spinor structures and define contractions
    class lorentz_vector
    {
        public: 
        lorentz_vector(){};

        lorentz_vector(std::array<complex,4> v)
        : _components(v)
        {};

        complex operator[](lorentz_index mu);

        inline complex operator[](int mu)
        { return operator[](static_cast<lorentz_index>(mu)); };

        lorentz_vector operator-();
        lorentz_vector operator+(lorentz_vector const & p);
        lorentz_vector operator=(lorentz_vector const & p);
        
        private:
        std::array<complex,4> _components = {};
    };
    
    // Interactions of Lorentz vectors and constants
    lorentz_vector operator*(complex c, lorentz_vector p);
    inline lorentz_vector operator*(lorentz_vector p, complex c){ return c *p; };
    inline lorentz_vector operator/(lorentz_vector p, complex c){ return (1./c) * p; };

    // Contractions between vectors
    complex contract(lorentz_vector x, lorentz_vector y);
    double  square(lorentz_vector x);

    // Now the covariants class assembles all relavant kinematic quantities and outputs 
    // them in the form they can be assembled in to amplitudes with Feynman rules
    class covariant
    {
        // ---------------------------------------------------------------------------

        public: 

        covariant(kinematics xkinem)
        : _cache(xkinem)
        {};

        // ---------------------------------------------------------------------------
        // Just like the amplitudes here we can save s, t, and helicities to avoid carrying them around

        inline void update(std::array<int,4> helicities, double s, double t)
        {
            // Save the helicities to avoid having to carry them around
            _lamB = helicities[0];
            _lamT = helicities[1];
            _lamX = helicities[2];
            _lamR = helicities[3];

            // Recalculate our cache if things have changed
            _cache.check(s, t);            
        };  

        // ---------------------------------------------------------------------------
        // These are the primary methods. They use the cached information and output
        // the appropriate Lorentz and Dirac structures

        // FOUR MOMENTA 
        // q Always refers to mesons, p always baryons

        lorentz_vector q();       // Beam momentum 
        lorentz_vector q_prime(); // Meson momentum 
        lorentz_vector p();       // Target momentum 
        lorentz_vector p_prime(); // Recoil momentum 

        // // ---------------------------------------------------------------------------
        // // Aliases for stuff!

        // // -----------------------------------------------------------------
        // // Four vectors of each particle
        // inline std::complex<double> beam_momentum(int mu)
        // {
        //     return _initial_state->q(mu, _s, 0.);
        // };
        // inline std::complex<double> meson_momentum(int mu)
        // {
        //     return _final_state->q(mu, _s, _theta);
        // };
        // inline std::complex<double> target_momentum(int mu)
        // {
        //     return _initial_state->p(mu, _s, 0.);
        // }
        // inline std::complex<double> recoil_momentum(int mu)
        // {
        //     return _final_state->p(mu, _s, _theta);
        // };

        // // -----------------------------------------------------------------
        // // Polarization vectors
        // inline std::complex<double> beam_polarization(int mu)
        // {
        //     return _beam_polarization->component(mu, _lam_gam, _s, 0.);
        // };
        // inline std::complex<double> meson_polarization(int mu)
        // {
        //     return _meson_polarization->conjugate_component(mu, _lam_vec, _s, _theta);
        // };

        // // Field strength tensors
        // inline std::complex<double> beam_field_tensor(int mu, int nu)
        // {
        //     return _beam_polarization->field_tensor(mu, nu, _lam_gam, _s, 0.);
        // }; 

        // // inline std::complex<double> field_tensor(int i, int j, int lambda, double s, double theta)
        // // {
        // //     std::complex<double> result;
        // //     result  = _state->q(i, s, theta) * component(j, lambda, s, theta);
        // //     result -= _state->q(j, s, theta) * component(i, lambda, s, theta);

        // //     return result;
        // // };

        // // -----------------------------------------------------------------
        // // Spinors 
        // inline std::complex<double> target_spinor(int i)
        // {
        //     return _target_spinor->component(i, _lam_tar, _s, 0.);
        // };

        // // If called with only one index assume we mean the dirac spinor
        // inline std::complex<double> recoil_spinor(int i)
        // {
        //     return _recoil_spinor->adjoint_component(i, _lam_rec, _s, _theta);
        // };
        // // Else with two indices we use the rarita spinor
        // inline std::complex<double> recoil_spinor(int i, int mu)
        // {
        //     return _recoil_rarita_spinor->adjoint_component(i, mu, _lam_rec, _s, _theta);
        // };

        // // -----------------------------------------------------------------
        // // momentum transfer 4-vectors

        // inline std::complex<double> s_momentum(int mu)
        // {
        //     return _initial_state->q(mu, _s, 0.) + _initial_state->p(mu, _s, 0.);
        // };


        // inline std::complex<double> t_momentum(int mu)
        // {
        //     return _initial_state->q(mu, _s, 0.) - _final_state->q(mu, _s, _theta);
        // };

        // inline std::complex<double> u_momentum(int mu)
        // {
        //     return _initial_state->q(mu, _s, 0) - _final_state->p(mu, _s, _theta);
        // };

        // // -----------------------------------------------------------------
        // // Slashed quantities 

        // // Take in a Lorentz vector and output the contraction with the gamma vector
        // inline std::complex<double> slash(int i, int j, std::function<std::complex<double>(int)> vec)
        // {
        //     std::complex<double> result = 0.;
        //     for (int mu = 0; mu < 4; mu++)
        //     {
        //         std::complex<double> temp;
        //         temp  = GAMMA[mu][i][j];
        //         temp *= METRIC[mu];
        //         temp *= vec(mu);

        //         result += temp;
        //     }

        //     return result;
        // };

        // inline std::complex<double> slashed_beam_polarization(int i, int j)
        // {
        //     return slash(i, j, [this](int mu){ return beam_polarization(mu); });
        // };

        // inline std::complex<double> slashed_meson_polarization(int i, int j)
        // {
        //     return slash(i, j, [this](int mu){ return meson_polarization(mu); });
        // };
        
        // inline std::complex<double> slashed_s_momentum(int i, int j)
        // {
        //     return slash(i, j, [this](int mu){ return s_momentum(mu); });
        // };

        // inline std::complex<double> slashed_u_momentum(int i, int j)
        // {
        //     return slash(i, j, [this](int mu){ return u_momentum(mu); });
        // };

        // inline std::complex<double> slashed_beam_momentum(int i, int j)
        // {
        //     return slash(i, j, [this](int mu){ return beam_momentum(mu); });
        // };

        // inline std::complex<double> slashed_meson_momentum(int i, int j)
        // {
        //     return slash(i, j, [this](int mu){ return meson_momentum(mu); });
        // };

        // ---------------------------------------------------------------------------
        private:

        // Core member structures that make up any other combination in our amplitudes

        // four_momenta structure is shared_ptr because the rest depend on it
        // std::shared_ptr<four_momenta>        _initial, _final;

        // The rest are unique_ptrs so they get deleted appropriately
        // std::unique_ptr<polarization_vector> _beam_polarization, _meson_polarization;
        // std::unique_ptr<dirac_spinor>        _target_spinor,     _recoil_spinor;
        // std::unique_ptr<rarita_spinor>       _recoil_rarita_spinor;

        // Save helicities current helicities so they dont need to be passed around
        int _lamB, _lamT, _lamX, _lamR; 

        struct cache
        {
            cache(kinematics xkinem)
            : _kinematics(xkinem)
            {};

            kinematics _kinematics;

            void check(double s, double t);
            void recalculate(double s);

            // Caching tolerance
            constexpr static double _tolerance = 1.E-5;
            inline bool unsynced(double x, double y)
            { return ( abs(x - y) < _tolerance ) ? false : true; };

            // Primary energy and angular variables
            double _s, _theta;

            // Masses 
            double _mB, _mT, _mX, _mR;

            // Energy of all three particles
            complex _EB, _ET, _EX, _ER;

            // 3-momenta of initial and final states
            complex _qi, _qf;

            // cos(theta) and sin(theta) of s-channel scattering angle
            double _cos, _sin;
        };

        cache _cache;
    };
};

#endif