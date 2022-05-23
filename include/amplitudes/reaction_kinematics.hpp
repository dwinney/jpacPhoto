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


#ifndef KINEMATICS
#define KINEMATICS

#include "constants.hpp"
#include "misc_math.hpp"
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
    // forthe reaction. Here you'll find the momenta and energies of all particles, angles,
    // invariant, etc.
    //
    // These are the basis for more complicated structures such as covariants or the
    // amplitudes themselves
    // ---------------------------------------------------------------------------
    
    class reaction_kinematics
    {
        public: 
        // Empty constructor
        // defaults to compton scattering: gamma p -> gamma p
        reaction_kinematics()
        {};

        // Constructor to fully specify the final state
        reaction_kinematics(double mX, double mR = M_PROTON)
        : _mX(mX), _mX2(mX*mX), _mR(mR), _mR2(mR*mR),
          _mB(0.), _mB2(0.), _mT(M_PROTON), _mT2(M2_PROTON)
        {};

        // Constructor to fully specify both final and initial states
        reaction_kinematics(double mB, double mT, double mX, double mR)
        : _mX(mX), _mX2(mX*mX), _mR(mR), _mR2(mR*mR),
          _mB(mB), _mB2(mB*mB), _mT(mT), _mT2(mT*mT)
        { if (mB > 0.) _photon = false; };

        // ---------------------------------------------------------------------------
        // Masses 

        inline double Wth(){ return (_mX + _mR); }; // square root of the threshold
        inline double sth(){ return Wth() * Wth(); }; // final state threshold

        // Assessor functions for masses
        inline double get_meson_mass() { return _mX; };
        inline double get_recoil_mass(){ return _mR; };
        inline double get_target_mass(){ return _mT; };
        inline double get_beam_mass()  { return _mB; };
        
        // Setters for masses
        inline void set_meson_mass(double x) {_mX = x; _mX2 = x*x; };
        inline void set_target_mass(double x){_mT = x; _mT2 = x*x; };
        inline void set_recoil_mass(double x){_mR = x; _mR2 = x*x; };
        inline void set_beam_mass(double x){  _mB = x; _mB2 = x*x;};

        // Get whether current kinematics has photon beam
        inline bool is_photon(){ return _photon; };
        inline void set_Q2(double x)
        {
            if (!is_photon()) 
            {
                std::cout << "Trying to set Q2 without initializing as a photon first! \n";
                std::cout << "Initialize reaction_kinematics with massless beam to indicate photon then use set_Q2()!" << std::endl;
                exit(1);
            }
            _virtual = true; 
            set_beam_mass(x);
        }

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
        // Beam [0], Target [1], Produced Meson [2], Recoil Baryon [3]
        int _nAmps = 24; 
        std::vector< std::array<int, 4> > _helicities = SPIN_ONE_HELICITIES;

        //--------------------------------------------------------------------------
        // Other quantities

        // Moduli of the initial and final state 3-momenta 

        inline double initial_momentum(double s)
        {
            return sqrt( Kallen(s, _mB2, _mT2)) / (2. * sqrt(s));
        };

        inline double final_momentum(double s)
        {
            return sqrt( Kallen(s, _mR2, _mX2)) / (2. * sqrt(s));
        };

        // Energies of all the particles

        inline double beam_energy(double s)
        {
            return (s - _mT2 + _mB2) / sqrt(4. * s);
        };

        inline double target_energy(double s)
        {
            return (s + _mT2 - _mB2) / sqrt(4. * s);
        };

        inline double meson_energy(double s)
        {
            return (s - _mR2 + _mX2) / sqrt(4. * s);
        };

        inline double recoil_energy(double s)
        {
            return (s + _mR2 - _mX2) / sqrt(4. * s);
        };

        // Get s-channel scattering angle from invariants
        inline double z_s(double s, double t)
        {
            std::complex<double> qdotqp = initial_momentum(s) * final_momentum(s);
            std::complex<double> E1E3   = beam_energy(s) * meson_energy(s);

            double result = t - _mX2 - _mB2 + 2.*real(E1E3);
            result /= 2. * real(qdotqp);

            return result;
        };

        // Scattering angle in the s-channel
        // Use TMath::ACos instead of std::acos because its safer at the end points
        inline double theta_s(double s, double t)
        {
            double zs = z_s(s, t);
            if (std::abs(zs - 1) < 1.E-5){zs =  1.;}
            if (std::abs(zs + 1) < 1.E-5){zs = -1.;}
            return TMath::ACos( zs );
        };

        // Invariant variables
        inline double t_man(double s, double theta)
        {
            std::complex<double> qdotqp = initial_momentum(s) * final_momentum(s);
            std::complex<double> E1E3   = beam_energy(s) * meson_energy(s);

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


        // Phase relating lambda_gamma = +1 and lambda_gamma = -1 amplitudes 
        // Depends on the channel with respect to which the helicities are defined
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

                    lam =  (2 * helicities[0] - helicities[1]);
                    lamp = (2 * helicities[2] - helicities[3]);

                    break;
                }
                case HELICITY_CHANNEL::T :
                {
                    s_b =  2*_jp[0];    eta_b = _jp[1] * pow(-1, _jp[0]);   // produced meson
                    s_c =  1;           eta_c = 1;                          // proton
                    s_d =  1;           eta_d = 1;                          // recoil baryon

                    lam =  (2 * (helicities[0] - helicities[2]));
                    lamp = (helicities[1] - helicities[3]);

                    break;
                }
                case HELICITY_CHANNEL::U :
                {
                    s_b =  1;           eta_b = 1;                          // recoil baryon
                    s_c =  1;           eta_c = 1;                          // proton
                    s_d =  2*_jp[0];    eta_d = _jp[1] * pow(-1, _jp[0]);   // produced meson

                    lam =  (2 * helicities[0] - helicities[3]);
                    lamp = (2 * helicities[2] - helicities[1]);

                    break;
                }

                default: { return 0.; }
            };

            int eta = eta_a * eta_b * eta_c * eta_d *  pow(-1, (lam - lamp)/2) * pow(-1., (s_c + s_d - s_a - s_b)/2);

            return double(eta);
        };

        private:

        // Masses are private to prevent them from being changed mid calculation 
        // Instead should be manipulated with public accessor and settor above
        bool _photon = true, _virtual = false;    // whether we have a photon and if its virtual
        double _mB = 0.,       _mB2 = 0.;         // mass and mass squared of the "beam" 
        double _mX = 0.,       _mX2 = 0.;         // mass and mass squared of the produced particle
        double _mT = M_PROTON, _mT2 = M2_PROTON;  // mass of the target, assumed to be proton unless overriden
        double _mR = M_PROTON, _mR2 = M2_PROTON;  // mass of the recoil baryon, assumed to be proton unless overriden
    };
};

#endif
