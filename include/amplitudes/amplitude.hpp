// Abstract class for an amplitude. Used so we can easily build observables
// as the incoherent sum of amplitudes in s, t, and u channels.
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef AMPLITUDE
#define AMPLITUDE

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
#include "covariant_kinematics.hpp"

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
        amplitude(reaction_kinematics * xkinem, std::string amp_name, std::string id = "", int n = 0)
        : _kinematics(xkinem), _identifier(id), _classname(amp_name)
        {
            _covariants = new covariant_kinematics(xkinem);
        };

        // Destructor needs to delete the covariants object
        ~amplitude()
        {
            delete _covariants;
        }

        // Kinematics object for thresholds and etc.
        reaction_kinematics * _kinematics;

        // Every amplitude gets a covariant_kinematics instance even if they dont use it
        covariant_kinematics * _covariants;

        // Allow each amplitude to carry a string id
        // This is user set and can help differentiate amplitudes
        // e.g. in plot legends
        inline void set_id(std::string id) { _identifier = id; };
        inline std::string get_id() { return _identifier; };

        // How the calculate the helicity amplitude
        // Must be given a specific implementation in a user derived class
        virtual std::complex<double> helicity_amplitude(std::array<int, 4> helicities, double s, double t) = 0;

        // Constant string which is used to differenciate derived classes 
        inline std::string amplitude_name()
        {
            return this->_classname;
        };
        
        // ---------------------------------------------------------------------------
        // Integer flag each amplitude carries for debugging options 
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

        virtual int parity_phase(std::array<int,4> helicities)
        {
            return 0;
        };

        // ---------------------------------------------------------------------------
        // nParams error message
        inline void set_nParams(int N){ _nParams = N; };

        // ---------------------------------------------------------------------------
        // Each amplitude must supply a function which returns a vector of allowed 2-tuples {J, P}
        virtual std::vector<std::array<int,2>> allowedJP() = 0;
        virtual std::vector<std::array<int,2>> allowedJP_Regge()
        { 
            return {};
        };
        
        // Small function just to know whether a given instance is a sum 
        inline bool is_sum(){ return _isSum; };

        // Update the cache
        void update_cache(double s, double t);
        std::complex<double> get_cached_helicity_amplitude(int i)
        {
            return _cached_helicity_amplitude[i];
        };

        // ---------------------------------------------------------------------------
        protected:
        
        inline void check_nParams(std::vector<double> params)
        {
            if (params.size() != _nParams)
            {
                std::cout << "\nWarning! Invalid number of parameters (" << params.size() << ") passed to " << _identifier << ".\n";
            }
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

        // Private string identifier to return which class the called object instance belongs to
        // used to differenciate different derived classes 
        // i.e. different physics amplitudes
        std::string _classname = "";

        // Debugging flag for each amplitude
        int _debug = 0;

        // Some saveable string by which to identify the amplitude
        std::string _identifier;

        // Bool only to differentiate the amplitude_sum class from everything else. 
        // This class gets special treatment in the helicity amplitude caching
        bool _isSum = false;

        // Number of parameters this amplitude takes in 
        int _nParams = 0;

        // ---------------------------------------------------------------------------
        // Each amplitude may save these values to avoid having to pass them as arguments
        // at each amplitude evaluation step

        inline void update(std::array<int,4> helicities, double s, double t)
        {
            _lam_gam = helicities[0];
            _lam_tar = helicities[1];
            _lam_vec = helicities[2];
            _lam_rec = helicities[3];

            _s = s; _t = t, _theta = _kinematics->theta_s(s, t);

            _mX = _kinematics->get_meson_mass();
            _mB = _kinematics->get_beam_mass();
            _mT = _kinematics->get_target_mass();
            _mR = _kinematics->get_recoil_mass();

            _covariants->update(helicities, s, t);
        };  

        // Energies
        double _s, _t, _theta;

        // Helicities
        int _lam_gam, _lam_vec;  
        int _lam_rec, _lam_tar;

        // Masses
        double _mB, _mX, _mR, _mT; 

        // ---------------------------------------------------------------------------
        // If helicity amplitudes have already been generated for a value of mV, s, t 
        // store them
        double _cache_tolerance = 1.E-4;
        double _cached_mX = 0., _cached_s = 0., _cached_t = 0.;
        std::vector<std::complex<double>> _cached_helicity_amplitude;
    };
};

#endif
