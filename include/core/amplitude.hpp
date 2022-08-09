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
        inline void set_debug(int d){ _debug = d; };

        // Force amplitude to evaluate using s-channel covariants if possible
        inline void force_covariant(bool x){ _useCovariant = x; };

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

        // this calculates the SDME in the amplitudes 'natural' frame 
        // i.e. doesnt require any additional rotations
        std::complex<double> SDME(int alpha, int lam, int lamp, double s, double t);

        // Rotate the SDMEs to another frame by and angle theta
        std::complex<double> rotated_SDME(int alpha, int lam, int lamp, double s, double t, double theta);
       
        // These check what the natural frame and include rotations if required to ensure
        // SDME is defined in the Helicity or Gottfried-Jackson frames
        std::complex<double> SDME_H (int alpha, int lam, int lamp, double s, double t);
        std::complex<double> SDME_GJ(int alpha, int lam, int lamp, double s, double t);

        // Beam Asymmetries
        double beam_asymmetry_y(double s, double t);    // Along the y direction
        double beam_asymmetry_4pi(double s, double t);  // integrated over decay angles

        // Parity asymmetry
        double parity_asymmetry(double s, double t);

        // ---------------------------------------------------------------------------
        // nParams error message
        inline void set_nParams(int N){ _nParams = N; };

        // ---------------------------------------------------------------------------
        // Each amplitude must supply a function which returns a vector of allowed 2-tuples {J, P}
        virtual std::vector<std::array<int,2>> allowed_meson_JP()  = 0;
        virtual std::vector<std::array<int,2>> allowed_baryon_JP() = 0;
        
        // Small function just to know whether a given instance is a sum 
        inline bool is_sum(){ return _isSum; };

        // Update the cache
        void update_cache(double s, double t);
        std::complex<double> get_cached_helicity_amplitude(int i){ return _cached_helicity_amplitude[i]; };

        // Each amplitude needs to define which frame the helicities are defined in
        virtual helicity_channel helicity_CM_frame() = 0;

        // With this phase only half of the helicity amplitudes need to be calculated and cached
        int parity_phase(std::array<int,4> helicities)
        {
            return _kinematics->parity_phase(helicities, this->helicity_CM_frame());
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
        inline void check_JP(reaction_kinematics * kinem)
        {
            // Get all the allowed JP's from the amplitude
            std::vector<std::array<int,2>> allowed_meson_JP, allowed_baryon_JP;
            allowed_meson_JP  = this->allowed_meson_JP();
            allowed_baryon_JP = this->allowed_baryon_JP();

            // Grab the requested JPs from the reaction_kinematics
            std::array<int,2> requested_meson_JP, requested_baryon_JP;
            requested_meson_JP  = kinem->get_meson_JP();
            requested_baryon_JP = kinem->get_baryon_JP();

            // Check if they are allowed
            bool meson_passes, baryon_passes;
            meson_passes  = !(std::find(allowed_meson_JP.begin(),  allowed_meson_JP.end(),  requested_meson_JP)  == allowed_meson_JP.end());
            baryon_passes = !(std::find(allowed_baryon_JP.begin(), allowed_baryon_JP.end(), requested_baryon_JP) == allowed_baryon_JP.end());

            if (!meson_passes || !baryon_passes)
            {
                std::cout << "Fatal error! Amplitude " << _classname << " not available for meson (J, P) = (" << requested_meson_JP[0] << ", " << requested_meson_JP[1] << ")";
                std::cout << " and baryon (J, P) = (" << requested_baryon_JP[0] << "/2, " << requested_baryon_JP[1] << ")! Quitting..." << std::endl;
                exit(1);
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

        // Flag to force amplitude to use covariant evaluation
        bool _useCovariant = false;

        // Flag to let amplitude know if it's a reggeized amplitude or not
        bool _reggeized = false;

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
        double _cached_mX = 0., _cached_mR = 0.;
        double _cached_s = 0., _cached_t = 0.;
        std::vector<std::complex<double>> _cached_helicity_amplitude;
    };
};

#endif
