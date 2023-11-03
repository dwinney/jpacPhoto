// Class to contain all relevant kinematic quantities. 
// (e.g. angles, momenta, helicity info)
//
// Unlike single_meson::kinematics this class represents a 2->3 process
// with two pseudo-scalars and a proton in the final state
// Generalizations to more complicated final states can be implemented lated
//
// Since we dont need to specify spins, everything is determined by each external mass:
// _mB = beam particle mass (default 0)
// _mT = target recoil mass (default M_PROTON)
// _m1 = "meson 1"
// _m2 = "meson 2"
// _mR = recoil baryon mass (default M_PROTON)
//
// --------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
//               University of Bonn (HISKP)
// Email:        winney@hiskp.uni-bonn.de
// --------------------------------------------------------------------------------

#ifndef KINEMATICS2_HPP
#define KINEMATICS2_HPP

#include "constants.hpp"
#include "helicities.hpp"

#include "TMath.h"

#include <array>
#include <vector>
#include <string>
#include <cmath>
#include <memory>

namespace jpacPhoto
{
    namespace two_meson
    {
        // ---------------------------------------------------------------------------
        // The reaction kinematics object is intended to have all relevant kinematic quantities
        // forthe reaction. Here you'll find the momenta and energies of all particles, angles,
        // invariant, etc.
        //
        // These are the basis for more complicated structures such as covariants or the
        // amplitudes themselves
        // ---------------------------------------------------------------------------
        
        // Forward declaration so we can rename ptr to kinematics as just kinematics
        // WE do this because we basically never want to work with a raw instance, but pass around a pointer
        class raw_kinematics;
        using kinematics = std::shared_ptr<jpacPhoto::two_meson::raw_kinematics>;

        // Constructors are made private to prevent raw_kinematics objects from being initialized
        // All user interactions should be through the pointer class kinematics
        class masses
        {
            private:
            
            friend kinematics new_kinematics(double, double);
            friend kinematics new_kinematics(double, double, double);
            friend kinematics new_kinematics(double, double, double, double, double);

            explicit masses(std::array<double,5> m)
            : _mB(m[0]), _mT(m[1]),
              _m1(m[2]), _m2(m[3]), _mR(m[4])
            {};

            public:

            double _mB = 0.;        // mass and mass squared of the "beam" 
            double _mT = M_PROTON;  // mass of the target, assumed to be proton unless overriden
            double _mR = M_PROTON;  // mass of the recoil baryon, assumed to be proton unless overriden

            // Meson masses
            double _m1 = M_PION, _m2 = M_PION;
        };


        // Specifying three masses determines the final state
        // Set proton gamma as initial state.
        // Third argument (recoil mass) defaults to proton
        inline kinematics new_kinematics(double m1, double m2, double m3 = M_PROTON)
        {
            return std::make_shared<raw_kinematics>(masses{{0., M_PROTON, m1, m2, m3}});
        };
        
        // With 5 arguments we can specify all masses
        inline kinematics new_kinematics(double m1, double m2, double m3, double m4, double m5)
        {
            return std::make_shared<raw_kinematics>(masses{{m1, m2, m3, m4, m5}});
        };

        class raw_kinematics
        {
            public: 

            // -----------------------------------------------------------------------
            // Only constructor requires the masses struct which is private
            // This is to prevent raw_kinematics objects from being created on the stack
            // We only want kinematics pointer objects to exist.

            raw_kinematics(masses m)
            : _masses(m),  _photon(!(m._mB>0.))
            {};

            // Delete copy and move operators so this is only passable by pointer
            raw_kinematics(const raw_kinematics&)                     = delete;
            const raw_kinematics &operator = (const raw_kinematics &) = delete;

            // The real "constructor" should be this factory method here.
            // This outputs shared_pointer objects which then are passed around

            // -----------------------------------------------------------------------
            // Accessors and setters for masses 

            inline double Wth(){ return (_masses._m1 + _masses._m2 + _masses._mR); }; // square root of the threshold
            inline double sth(){ return Wth() * Wth(); }; // final state threshold

            // Assessor functions for masses
            inline std::array<double,2> get_meson_masses() { return {_masses._m1, _masses._m2}; };
            inline double get_recoil_mass(){ return _masses._mR; };
            inline double get_target_mass(){ return _masses._mT; };
            inline double get_beam_mass()  { return _masses._mB; };
            
            // Setters for masses
            inline void set_meson_masses(std::array<double, 2> x) { _masses._m1 = x[0]; _masses._m2 = x[1]; };
            inline void set_target_mass(double x){ _masses._mT = x; };
            inline void set_recoil_mass(double x){ _masses._mR = x; };
            inline void set_beam_mass(double x)  { _masses._mB = x; };

            // Get whether current kinematics has photon beam
            inline bool is_photon() { return _photon;  };

            // ---------------------------------------------------------------------------
            // Accessing the helicity combinations 
            inline int N_amps(){ return _nAmps; };
            inline std::vector<std::array<int,3>> helicities(){ return _helicities; };
            inline std::array<int,3> helicities(int i){ return _helicities[i]; };
            int helicity_index(std::array<int,3> hel);

            //--------------------------------------------------------------------------
            // Other quantities

            // -----------------------------------------------------------------------
            private:

            // Container class holding all the masses of the system
            // Additionally serves as a pass-key class in the constructor
            masses _masses;

            // Masses are private to prevent them from being changed mid calculation 
            // Instead should be manipulated with public accessor and settor above
            bool _photon = true;    // whether we have a photon beam

            // Helicity configurations
            // Assume the two mesons are scalars 
            // Can be extended later if needed
            int _nAmps = 8; 
            std::vector<std::array<int,3>> _helicities =
            {
            //  { γ,  p,  p'}
                { 1, -1, -1}, // 0
                { 1, -1,  1}, // 1
                { 1,  1, -1}, // 2
                { 1,  1,  1}, // 3
                {-1, -1, -1}, // 4 
                {-1, -1,  1}, // 5
                {-1,  1, -1}, // 6
                {-1,  1,  1}  // 7
            };
        };
    };
};

#endif
