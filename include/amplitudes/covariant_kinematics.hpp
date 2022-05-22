// Container class for covariant quantities!
// Since not all amplitudes need these they're seperated for convience
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef COVARIANTS
#define COVARIANTS

#include "reaction_kinematics.hpp"

namespace jpacPhoto
{
    class covariant_kinematics
    {
        // ---------------------------------------------------------------------------

        public: 

        covariant_kinematics(reaction_kinematics * xkinem)
        : _kinematics(xkinem)         
        {
            _initial_state = new two_body_state(xkinem->get_beam_mass(),  xkinem->get_target_mass());
            _final_state   = new two_body_state(xkinem->get_meson_mass(), xkinem->get_recoil_mass());

            _beam_polarization  = new polarization_vector(_initial_state);
            _meson_polarization = new polarization_vector(_final_state);

            _target_spinor = new dirac_spinor(_initial_state);
            _recoil_spinor = new dirac_spinor(_final_state);
        };

        ~covariant_kinematics()
        {
            delete _initial_state;
            delete _final_state;
            delete _beam_polarization;
            delete _meson_polarization;
            delete _target_spinor;
            delete _recoil_spinor;
        }; 

        // ---------------------------------------------------------------------------
        // Just like the amplitudes here we can save s, t, and helicities to avoid carrying them around

        inline void update(std::array<int,4> helicities, double s, double t)
        {
            _lam_gam = helicities[0];
            _lam_tar = helicities[1];
            _lam_vec = helicities[2];
            _lam_rec = helicities[3];

            _s = s; _t = t, _theta = _kinematics->theta_s(s, t);

            // Grab masses from the kinematics
            _mX = _kinematics->get_meson_mass();
            _mB = _kinematics->get_beam_mass();
            _mT = _kinematics->get_target_mass();
            _mR = _kinematics->get_recoil_mass();

            // Pass them into the appropriate structures
            _initial_state->set_meson_mass(_mB);
            _initial_state->set_baryon_mass(_mT);
            _final_state->set_meson_mass(_mX);
            _final_state->set_baryon_mass(_mR);
        };  

        // ---------------------------------------------------------------------------
        // Aliases for stuff!

        // Four vectors of each particle
        inline std::complex<double> beam_momentum(int mu)
        {
            return _initial_state->q(mu, _s, 0.);
        };
        inline std::complex<double> meson_momentum(int mu)
        {
            return _final_state->q(mu, _s, _theta);
        };
        inline std::complex<double> target_momentum(int mu)
        {
            return _initial_state->p(mu, _s, 0.);
        }
        inline std::complex<double> recoil_momentum(int mu)
        {
            return _final_state->p(mu, _s, _theta);
        };

        // Polarization vectors
        inline std::complex<double> beam_polarization(int mu)
        {
            return _beam_polarization->component(mu, _lam_gam, _s, 0.);
        };
        inline std::complex<double> meson_polarization(int mu)
        {
            return _meson_polarization->conjugate_component(mu, _lam_vec, _s, _theta);
        };

        // Field strength tensors
        inline std::complex<double> beam_field_tensor(int mu, int nu)
        {
            return _beam_polarization->field_tensor(mu, nu, _lam_gam, _s, 0.);
        }; 

        // Spinors 
        inline std::complex<double> target_spinor(int i)
        {
            return _target_spinor->component(i, _lam_tar, _s, 0.);
        };
        inline std::complex<double> recoil_spinor(int i)
        {
            return _recoil_spinor->adjoint_component(i, _lam_rec, _s, _theta);
        };


        // momentum transfer 4-vectors
        inline std::complex<double> t_momentum(int mu)
        {
            return _initial_state->q(mu, _s, 0.) - _final_state->q(mu, _s, _theta);
        };

        inline std::complex<double> u_momentum(int mu)
        {
            return _initial_state->q(mu, _s, 0) - _final_state->p(mu, _s, _theta);
        };

        // ---------------------------------------------------------------------------
        private:
        
        // Pointer to kinematics object to get masses and things
        reaction_kinematics * _kinematics;

        // Core member structures
        two_body_state * _initial_state,  * _final_state;
        polarization_vector * _beam_polarization, * _meson_polarization;
        dirac_spinor * _target_spinor, * _recoil_spinor;

        // Energies
        double _s, _t, _theta;

        // Helicities
        int _lam_gam, _lam_vec;  
        int _lam_rec, _lam_tar;

        // Masses
        double _mB, _mX, _mR, _mT; 
    };
};

#endif