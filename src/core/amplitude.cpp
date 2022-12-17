// The core object of interest is the amplitude class which encorporates
// the dynamical physics model for a process
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#include "amplitude.hpp"
#include "constants.hpp"
#include "kinematics.hpp"
#include <cstddef>
#include <string>

namespace jpacPhoto
{
    // ------------------------------------------------------------------------------
    // External methods for amplitude sums

    bool are_compatible(amplitude a, amplitude b)
    {
        // Two amplitudes are compatible if they point to the same
        // kinematics object (i.e. are synced together)
        bool same_kinem= (a->_kinematics == b->_kinematics);

        // But also they must have helicites defined in the same frame
        // TODO: add helicity crossing relations so this is not a necessity
        bool same_frame   = (a->native_helicity_frame() == b->native_helicity_frame());

        
        std::string error_msg = "Attempted to add incompatible amplitudes: " + a->id() + " and " + b->id() + "!";
        if (!same_frame) 
        {
            return error(error_msg + " (Helicites defined in different frames!)", false);
        };

        if (!same_kinem)
        {
            return error(error_msg + " (Contain different kinematics objects!)", false);
        };

        return true;
    };
    
    // Sum two amplitudes into a new amplitude 
    amplitude operator+(amplitude a, amplitude b)
    {
        if (!are_compatible(a,b)) return nullptr;
        
        // New id will be the "sum" of the individual id's
        std::string id = a->id() + " + " + b->id(); 

        // vector of subamps
        std::vector<amplitude> amps = {a, b};
        
        // Initialize and return a new amplitude
        return std::make_shared<raw_amplitude>(amplitude_key(), amps, id);
    };

    // Take an existing sum and add a new amplitude
    void operator +=(amplitude a, amplitude b)
    {
        if ( !a->is_sum() )
        {
            warning("amplitude::operator+= ", 
                    "Attempting to sum two non-sums! \nPlease initialize a sum of amplitudes with = first (auto c = a + b;) then use the += operator (c += d;)!");
            return;
        };  

        a->add(b);
    };

    // ------------------------------------------------------------------------------
    // Internal methods for amplitude sums

    // Check if a given amplitude is compatible to be summed with *this
    bool raw_amplitude::is_compatible(amplitude new_amp)
    {
        // Two amplitudes are compatible if they point to the same
        // kinematics object (i.e. are synced together)
        bool same_kinem = (new_amp->_kinematics == _kinematics);

        // But also they must have helicites defined in the same frame
        // TODO: add helicity crossing relations so this is not a necessity
        bool same_frame   = (new_amp->native_helicity_frame() == native_helicity_frame());

        
        std::string error_msg = "Attempted to add incompatible amplitudes: " + id() + " and " + new_amp->id() + "!";
        if (!same_frame) 
        {
            return error(error_msg + " (Helicites defined in different frames!)", false);
        };

        if (!same_kinem)
        {
            return error(error_msg + " (Contain different kinematics objects!)", false);
        };
        return true;
    };

    // Add a new amplitude to the list
    void raw_amplitude::add(amplitude new_amp)
    {
        // If our list is empty, capture the first entry we encounter
        if (_subamplitudes.size() == 0)
        {
            _kinematics = new_amp->_kinematics;
            _subamplitudes.push_back(new_amp);
            return;
        };

        // Each subsequent entry gets compared to the first
        if (is_compatible(new_amp))
        {
            // If new_amp is itself a sum of amps,
            // reach in and pull out the components 
            if (new_amp->is_sum())
            {
                add(new_amp->_subamplitudes);
                return;
            };

            _subamplitudes.push_back(new_amp);
            _N_pars += new_amp->N_pars();
            return;
        };
    };

    // Add a whole vector, this just loops over each entry
    void raw_amplitude::add(std::vector<amplitude> new_amps)
    {
        for (auto amp : new_amps)
        {
            add(amp);
        };
    };

    // ------------------------------------------------------------------------------
    // Updators for the different caches we have

    // Save current kinematic quantities from _kinematics
    void raw_amplitude::store(std::array<int,4> helicities, double s, double t)
    {
        _lamB = helicities[0];
        _lamT = helicities[1];
        _lamX = helicities[2];
        _lamR = helicities[3];

        _s = s; _t = t; _theta = _kinematics->theta_s(s,t);

        _mB = _kinematics->get_beam_mass();
        _mT = _kinematics->get_target_mass();
        _mX = _kinematics->get_meson_mass();
        _mR = _kinematics->get_recoil_mass();

        return;
    };

    // Check for any changes within the set tolerance 
    // And recalculate helicity amplitudes if needed
    void raw_amplitude::update_cache(double s, double t)
    {
        bool s_changed  = !are_equal(_cached_s, s, _cache_tolerance);
        bool t_changed  = !are_equal(_cached_t, t, _cache_tolerance);
        bool mX_changed = !are_equal(_cached_mX, _kinematics->get_meson_mass(),  _cache_tolerance);
        bool mR_changed = !are_equal(_cached_mR, _kinematics->get_recoil_mass(), _cache_tolerance);

        bool need_update = s_changed || t_changed || mX_changed || mR_changed;

        if (need_update)
        {
            _cached_helicity_amplitudes.clear();

            // Total number of amplitudes
            int n = _kinematics->N_amps();

            // However we ill only calculate the first half
            // with the rest related by the parity_phase
            for (int i = 0; i < n; i++)
            {
                complex hel_amp;
                if (i < n/2) hel_amp = helicity_amplitude(_kinematics->helicities(i), s, t);
                else         hel_amp = parity_phase(i) * _cached_helicity_amplitudes[n - 1 - i];
                _cached_helicity_amplitudes.push_back(hel_amp); 
            }

            if (_cached_helicity_amplitudes.size() != n)
            {
                warning(name()+"::update_cache", "Cached size doesn't match expected number of helicity amplitude!");
            };
        };

        // Update the cache info as well
        _cached_mX = _kinematics->get_meson_mass();
        _cached_mR = _kinematics->get_recoil_mass();
        _cached_s  = s; _cached_t = t;
    };

    // ------------------------------------------------------------------------------
    // Phase properties of helicity amplitudes

    double raw_amplitude::parity_phase(std::array<int,4> helicities)
    {
        return _kinematics->parity_phase(helicities, this->native_helicity_frame());
    };

    double raw_amplitude::parity_phase(int i)
    {
        return _kinematics->parity_phase(i, this->native_helicity_frame());
    };

    // ------------------------------------------------------------------------------
    // Parameter handling

    // Ability to change number of expected parameters from inside a class
    void raw_amplitude::set_N_pars(int N)
    {
        _N_pars = N;
    };

    // Simple check that a given vector is of the expected size
    bool raw_amplitude::correct_size(std::vector<double> pars)
    {
        if (pars.size() != _N_pars)
        {
            return error(name()+"::set_parameters", "Number of parameters passed not the expected size!", false);
        };
        return true;
    };

    // ------------------------------------------------------------------------------
    // Compatability checks

    // At construction, check the kinematics with set quantum numbers is compatable
    void raw_amplitude::check_QNs(kinematics kinem)
    {
        // Get all the allowed JP's from the amplitude
        std::vector<std::array<int,2>> allowed_meson_JP, allowed_baryon_JP;
        allowed_meson_JP  = this->allowed_meson_JP();
        allowed_baryon_JP = this->allowed_baryon_JP();

        // Grab the requested JPs from the raw_kinematics
        std::array<int,2> requested_meson_JP, requested_baryon_JP;
        requested_meson_JP  = kinem->get_meson_JP();
        requested_baryon_JP = kinem->get_baryon_JP();

        // Check if they are allowed
        bool meson_fails, baryon_fails;
        meson_fails  = std::find(allowed_meson_JP.begin(),  allowed_meson_JP.end(),  requested_meson_JP)  == allowed_meson_JP.end();
        baryon_fails = std::find(allowed_baryon_JP.begin(), allowed_baryon_JP.end(), requested_baryon_JP) == allowed_baryon_JP.end();

        if (meson_fails)
        {
            warning(name()+"::check_QNs", "Requested meson quantum numbers (J=" + std::to_string(requested_meson_JP[0]) + ", P=" + std::to_string(requested_meson_JP[1])+") not available!");
        };
    
        if (baryon_fails)
        {
            warning(name()+"::check_QNs", "Requested baryon quantum numbers (J=" + std::to_string(requested_baryon_JP[0]) + "/2, P=" + std::to_string(requested_baryon_JP[1])+") not available!");
        };    
    };

    // ------------------------------------------------------------------------------
    // OBSERVABLES
    
    // Square of the spin averaged amplitude squared
    double raw_amplitude::probability_distribution(double s, double t)
    {
        // Check we have the right amplitudes cached
        update_cache(s, t);

        double sum = 0.;
        for (auto amp : _cached_helicity_amplitudes)
        {
            sum += std::real(amp * conj(amp));
        }

        return sum;
    };

    // Differential cross section dsigma / dt
    // in NANOBARN
    double raw_amplitude::differential_xsection(double s, double t)
    {
        if (s < _kinematics->sth()) return 0.;

        double sum = probability_distribution(s, t);

        double norm = 1.;
        norm /= 64. * PI * s;
        norm /= std::real(pow(_kinematics->initial_momentum(s), 2.));
        norm /= (2.56819E-6); // Convert from GeV^-2 -> nb

        // Average over initial helicities
        if (_kinematics->is_photon()) norm /= 4.;
        else                          norm /= 6.;

        return norm * sum;
    };

    // Integrated total cross-section
    // IN NANOBARN
    double raw_amplitude::integrated_xsection(double s)
    {
        if (s < _kinematics->sth()) return 0.;

        auto F = [&](double t)
        {
            return differential_xsection(s, t);
        };

        ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS61);
        ROOT::Math::Functor1D wF(F);
        ig.SetFunction(wF);

        double t_min = _kinematics->t_min(s);
        double t_max = _kinematics->t_max(s);

        return ig.Integral(t_max, t_min);
    };


};