// Class to contain all relevant kinematic quantities. The kinematics of the reaction
// gamma p -> X p' is entirely determined by specifying the mass of the vector particle.
//
// Additional options to include virtual photon and different baryons (e.g. gamma p -> X Lambda_c) 
// also available
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------


#ifndef _KINEMATICS_
#define _KINEMATICS_

#include "constants.hpp"
#include "misc_math.hpp"
#include "two_body_state.hpp"
#include "dirac_spinor.hpp"
#include "polarization_vector.hpp"
#include "helicities.hpp"

#include "TMath.h"

#include <array>
#include <vector>
#include <string>
#include <cmath>

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    // The reaction kinematics object is intended to have all relevant kinematic quantities
    // forthe reaction. Here you'll find the momenta and energies of all particles,
    //  spinors for the baryons and polarization vectors for the gamma and produced meson
    // ---------------------------------------------------------------------------

    class reaction_kinematics
    {
        public: 
        // Empty constructor,
        // defaults to compton scattering: gamma p -> gamma p
        reaction_kinematics()
        {
            _initial_state   = new two_body_state(0., M2_PROTON);
            _eps_gamma       = new polarization_vector(_initial_state);
            _target          = new dirac_spinor(_initial_state);

            _final_state     = new two_body_state(0., M2_PROTON);
            _eps_vec         = new polarization_vector(_final_state);
            _recoil          = new dirac_spinor(_final_state);
        };

        // Constructor with a set mX and JP
        // defaults to proton as baryon and real photon
        reaction_kinematics(double mX, std::array<int,2> jp = VECTOR)
        : _mX(mX), _mX2(mX*mX)
        {
            _initial_state   = new two_body_state(0., M2_PROTON);
            _eps_gamma       = new polarization_vector(_initial_state);
            _target          = new dirac_spinor(_initial_state);

            _final_state     = new two_body_state(mX*mX, M2_PROTON);
            _eps_vec         = new polarization_vector(_final_state);
            _recoil          = new dirac_spinor(_final_state);

            set_JP(jp[0], jp[1]);
        };

        // Constructor with a set mX and baryon mass mR
        // defaults to real photon
        reaction_kinematics(double mX, double mR)
        : _mX(mX), _mX2(mX*mX),
          _mR(mR), _mR2(mR*mR)
        {
            _initial_state   = new two_body_state(0., M2_PROTON);
            _eps_gamma       = new polarization_vector(_initial_state);
            _target          = new dirac_spinor(_initial_state);

            _final_state     = new two_body_state(mX*mX, mR*mR);
            _eps_vec         = new polarization_vector(_final_state);
            _recoil          = new dirac_spinor(_final_state);
        };

        // Constructor with a set mV and baryon mass mR
        // and massive incoming scatterer and target
        reaction_kinematics(double mX, double mR, double mT, double mB = 0.)
        : _mX(mX), _mX2(mX*mX),
          _mR(mR), _mR2(mR*mR),
          _mB(mB), _mB2(mB*mB),
          _mT(mT), _mT2(mT*mT)
        {
            if (mB > 0.) _photon = false;
            _initial_state   = new two_body_state(mB*mB, mT*mT);
            _eps_gamma       = new polarization_vector(_initial_state);
            _target          = new dirac_spinor(_initial_state);

            _final_state     = new two_body_state(mX*mX, mR*mR);
            _eps_vec         = new polarization_vector(_final_state);
            _recoil          = new dirac_spinor(_final_state);
        };

        // destructor
        ~reaction_kinematics()
        {
            delete _initial_state;
            delete _final_state;
            delete _eps_gamma;
            delete _eps_vec;
            delete _target;
            delete _recoil;
        }

        // ---------------------------------------------------------------------------
        // Masses
        bool _photon = true;
        double _mB = 0., _mB2 = 0.;       // mass and mass squared of the "beam" 
        double _mX = 0., _mX2 = 0.;       // mass and mass squared of the produced particle

        double _mT = M_PROTON, _mT2 = M2_PROTON;  // mass of the target, assumed to be proton unless overriden
        double _mR = M_PROTON, _mR2 = M2_PROTON;  // mass of the recoil baryon, assumed to be proton unless overriden

        inline double Wth(){ return (_mX + _mR); }; // square root of the threshold
        inline double sth(){ return Wth() * Wth(); }; // final state threshold

        // Change the meson mass
        inline void set_mX(double m)
        {
            _mX  = m;
            _mX2 = m*m;

            // also update the meson mass in two_body_state
            _final_state->set_mV2(m*m);
        };

        inline void set_mX2(double m2)
        {
            _mX  = sqrt(m2);
            _mX2 = m2;

            // also update the meson mass in two_body_state
            _final_state->set_mV2(m2);
        };

        // Change virtuality of the photon
        // Q2 > 0
        inline void set_Q2(double q2)
        {
            if (q2 < 0) { std::cout << "Caution! set_Q2(x) requires x > 0! \n"; }
            _mB2 = -q2;
            _initial_state->set_mV2(-q2);
        };

        // ---------------------------------------------------------------------------
        // Quantum numbers of produced meson. 
        std::array<int,2> _jp{{1,1}};
        inline void set_JP(int J, int P)
        { 
            _jp = {J, P};
            _helicities = get_helicities(J, _mB);
            _nAmps = _helicities.size();
        };
        
        inline void set_JP(std::array<int,2> jp)
        { 
            _jp = jp;
            _helicities = get_helicities(jp[0], _mB);
            _nAmps = _helicities.size();
        };

        // Helicity configurations
        // Defaults to spin-1
        // Photon [0], Incoming Proton [1], Produced meson [2], Outgoing Proton [3]
        int _nAmps = 24; 
        std::vector< std::array<int, 4> > _helicities = SPIN_ONE_HELICITIES;

        //--------------------------------------------------------------------------
        two_body_state * _initial_state,  * _final_state;
        polarization_vector * _eps_vec, * _eps_gamma;
        dirac_spinor * _target, * _recoil;

        // Get s-channel scattering angle from invariants
        inline double z_s(double s, double t)
        {
            std::complex<double> qdotqp = _initial_state->momentum(s) * _final_state->momentum(s);
            std::complex<double> E1E3   = _initial_state->energy_V(s) * _final_state->energy_V(s);

            double result = t - _mX2 - _mB2 + 2.*real(E1E3);
            result /= 2. * real(qdotqp);

            return result;
        };

        // Scattering angle in the s-channel
        // Use TMath::ACos instead of std::acos because its safer at the end points
        inline double theta_s(double s, double t)
        {
            return TMath::ACos( z_s(s, t) );
        };

        // Invariant variables
        inline double t_man(double s, double theta)
        {
            std::complex<double> qdotqp = _initial_state->momentum(s) * _final_state->momentum(s);
            std::complex<double> E1E3   = _initial_state->energy_V(s) * _final_state->energy_V(s);

            return _mX2 + _mB2 - 2. * real(E1E3) + 2. * real(qdotqp) * cos(theta);
        };

        inline double u_man(double s, double theta)
        { 
            return _mX2 + _mB2 + _mT2 + _mR2 - s - t_man(s, theta);
        };

        // Scattering angles in t and u channel frames
        inline std::complex<double> z_t(double s, double theta)
        {
            double t = t_man(s, theta);
            double u = u_man(s, theta);

            std::complex<double> result;
            result  = t * (s - u) + (_mB2 - _mX2) * (_mT2 - _mR2);
            result /=  sqrt(XR * Kallen(t, _mX2, _mB2)) * sqrt(XR * Kallen(t, _mT2, _mR2));

            return result;
        };

        inline std::complex<double> z_u(double s, double theta)
        {
            double t = t_man(s, theta);
            double u = u_man(s, theta);

            std::complex<double> result;
            result  = u * (t - s) + (_mB2 - _mR2) * (_mT2 - _mX2);
            result /=  sqrt(XR * Kallen(u, _mR2, _mB2)) * sqrt(XR * Kallen(u, _mT2, _mX2));

            return result;
        };

        // momentum transfer 4-vectors
        inline std::complex<double> t_exchange_momentum(int mu, double s, double theta)
        {
            std::complex<double> qGamma_mu, qA_mu;
            qGamma_mu   = _initial_state->q(mu, s, 0.);
            qA_mu       = _final_state->q(mu, s, theta);

            return (qGamma_mu - qA_mu);
        };

        inline std::complex<double> u_exchange_momentum(int mu, double s, double theta)
        {
            std::complex<double> qGamma_mu, qRec_mu;
            qGamma_mu   = _initial_state->q(mu, s, 0);
            qRec_mu     = _final_state->p(mu, s, theta + PI);

            return qGamma_mu - qRec_mu;
        };

        inline double parity_phase(std::array<int, 4> helicities, HELICITY_CHANNEL channel)
        {
            int s_a, s_b, s_c, s_d;
            int eta_a, eta_b, eta_c, eta_d;
            int lam, lamp;

            // a is always the photon
            s_a = 2; eta_a = 1; // spin multiplied by two because of spin 1/2 baryons

            switch (channel)
            {
                case HELICITY_CHANNEL::S :
                {
                    s_b =  1;           eta_b = 1;                         // proton
                    s_c =  2*_jp[0];    eta_c = _jp[1] * pow(-1, _jp[0]);  // produced meson
                    s_d =  1;           eta_d = 1;                         // recoil baryon

                    lam =  double(2 * helicities[0] - helicities[1]);
                    lamp = double(2 * helicities[2] - helicities[3]);

                    break;
                }
                case HELICITY_CHANNEL::T :
                {
                    s_b =  2*_jp[0];    eta_b = _jp[1] * pow(-1, _jp[0]);   // produced meson
                    s_c =  1;           eta_c = 1;                          // proton
                    s_d =  1;           eta_d = 1;                          // recoil baryon

                    lam =  double(2 * (helicities[0] - helicities[2]));
                    lamp = double(helicities[1] - helicities[3]);

                    break;
                }
                case HELICITY_CHANNEL::U :
                {
                    s_b =  1;           eta_b = 1;                          // recoil baryon
                    s_c =  1;           eta_c = 1;                          // proton
                    s_d =  2*_jp[0];    eta_d = _jp[1] * pow(-1, _jp[0]);   // produced meson

                    lam =  double(2 * helicities[0] - helicities[3]);
                    lamp = double(2 * helicities[2] - helicities[1]);

                    break;
                }

                default:
                {
                    return 0.;
                }
            };

            int eta = eta_a * eta_b * eta_c * eta_d *  pow(-1, (lam - lamp)/2) * pow(-1., (s_c + s_d - s_a - s_b)/2);

            return double(eta);
        };
    };
};

#endif
