// Abstract class for an amplitude. Used so we can easily build observables
// as the incoherent sum of amplitudes in s, t, and u channels.
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _AMPLITUDE_
#define _AMPLITUDE_

// ---------------------------------------------------------------------------
// Abstract class to define helicity amplitudes. This will allow multiple different
// classes (for s, t, and u- channels but also multiple contibutions in each channel)
// to be added together and evaluated in observables.
//
// Any generic amplitude needs a reaction_kinematics object
// and a way to evaluate the helicity amplitude for given set of helicities,
// CoM energy and scattering angle.
// ---------------------------------------------------------------------------

#include "reaction_kinematics.hpp"

#include "Math/GSLIntegrator.h"
#include "Math/IntegrationTypes.h"
#include "Math/Functor.h"

#include <string>
#include <algorithm>

namespace jpacPhoto
{
    class amplitude
    {
        public:
        // Constructor with nParams for backward compatibility (now depricated)
        amplitude(reaction_kinematics * xkinem, std::string id = "", int n = 0)
        : _kinematics(xkinem), _identifier(id)
        {};

        // Kinematics object for thresholds and etc.
        reaction_kinematics * _kinematics;

        // saved energies and angle 
        double _s, _t, _theta;

        // Some saveable string by which to identify the amplitude
        std::string _identifier;

        // How the calculate the helicity amplitude
        // Must be given a specific implementation in a user derived class
        virtual std::complex<double> helicity_amplitude(std::array<int, 4> helicities, double s, double t) = 0;

        // Sum amplitudes get special treatment (for example in the check_cache() method)
        bool _isSum = false;
        
        // ---------------------------------------------------------------------------
        int _debug = 0;
        inline void set_debug(int d)
        {
            _debug = d;
        };

        // ---------------------------------------------------------------------------
        // Observables
        // Evaluatable in terms of s and t or an event object (see reaction_kinematics.hpp)

        // Modulus of the amplitude summed over all helicity combinations
        double probability_distribution(double s, double t);

        // Differential and total cross-section
        double differential_xsection(double s, double t);

        // integrated crossection
        double integrated_xsection(double s);

        // Spin asymmetries
        double A_LL(double s, double t); // Beam and target
        double K_LL(double s, double t); // Beam and recoil

        // Spin density matrix elements
        std::complex<double> SDME(int alpha, int lam, int lamp, double s, double t);

        // Beam Asymmetries
        double beam_asymmetry_y(double s, double t);    // Along the y direction
        double beam_asymmetry_4pi(double s, double t);  // integrated over decay angles

        // Parity asymmetry
        double parity_asymmetry(double s, double t);

        // ---------------------------------------------------------------------------
        // If helicity amplitudes have already been generated for a value of mV, s, t 
        // store them
        double _cached_mX2 = 0., _cached_s = 0., _cached_t = 0.;
        std::vector<std::complex<double>> _cached_helicity_amplitude;

        void check_cache(double s, double t);

        virtual int parity_phase(std::array<int,4> helicities)
        {
            return 0;
        };

        // ---------------------------------------------------------------------------
        // nParams error message
        int _nParams = 0;
        inline void set_nParams(int N){ _nParams = N; };
        inline void check_nParams(std::vector<double> params)
        {
            if (params.size() != _nParams)
            {
                std::cout << "\nWarning! Invalid number of parameters (" << params.size() << ") passed to " << _identifier << ".\n";
            }
        };

        // ---------------------------------------------------------------------------
        // Each amplitude must supply a function which returns a vector of allowed 2-tuples {J, P}
        virtual std::vector<std::array<int,2>> allowedJP() = 0;
        virtual std::vector<std::array<int,2>> allowedJP_Regge()
        { 
            return {};
        };
        
        // Allowed JP error message
        inline void check_JP(std::array<int,2> JP, bool REGGE = false)
        {
           std::vector<std::array<int,2>> allowed_JP;
           (REGGE) ? (allowed_JP = allowedJP_Regge()) : (allowed_JP = allowedJP());
           
            if (std::find(allowed_JP.begin(), allowed_JP.end(), JP) == allowed_JP.end())
            {
                std::cout << "Error! Amplitude for spin: " << JP[0] << " and parity " << JP[1] << " for " << _identifier << " and REGGE = " << REGGE << " is unavailable.\n";
                exit(0);
            }      
        };
    };
};

#endif
