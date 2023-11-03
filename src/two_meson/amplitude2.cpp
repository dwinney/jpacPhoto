// Analogue to the jpacPhoto::amplitude class in /core except extended to 2 meson
// final states. This is a generic class and requires explicit implementations
// given by the user
//
// --------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
//               University of Bonn (HISKP)
// Email:        winney@hiskp.uni-bonn.de
// --------------------------------------------------------------------------------

#include "amplitude2.hpp"

namespace jpacPhoto
{
    namespace two_meson
    {
        // ------------------------------------------------------------------------------
        // External methods for amplitude sums

        bool are_compatible(amplitude a, amplitude b)
        {
            // Two amplitudes are compatible if they point to the same
            // kinematics object (i.e. are synced together)
            bool same_kinem= (a->_kinematics == b->_kinematics);
            
            std::string error_msg = "Attempted to add incompatible amplitudes: " + a->id() + " and " + b->id() + "!";
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

            // Extract the sub amplitudes from each vector being summed
            // This allows the saved pointers to all ways be "base" level amplitudes in stead of other sums
            std::vector<amplitude> from_a = (a->is_sum()) ? a->_subamplitudes : std::vector<amplitude>{{a}};
            std::vector<amplitude> from_b = (b->is_sum()) ? b->_subamplitudes : std::vector<amplitude>{{b}};
            from_a.insert(from_a.end(), from_b.begin(), from_b.end());

            // Initialize and return a new amplitude
            return std::make_shared<raw_amplitude>(amplitude_key(), from_a, id);
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
            
            std::string error_msg = "Attempted to add incompatible amplitudes: " + id() + " and " + new_amp->id() + "!";
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
                _N_pars     = new_amp->N_pars();
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
        void raw_amplitude::store(std::array<int,3> helicities, double s, double s12, double t)
        {
            _lamB = helicities[0];
            _lamT = helicities[1];
            _lamR = helicities[2];

            _s = s; _s12 = s12; _t = t;

            _mB = _kinematics->get_beam_mass();
            _mT = _kinematics->get_target_mass();
            _mR = _kinematics->get_recoil_mass();

            std::array<double,2> m = _kinematics->get_meson_masses();
            _m1 = m[0]; _m2 = m[1];
            return;
        };

        // Check for any changes within the set tolerance 
        // And recalculate helicity amplitudes if needed
        void raw_amplitude::update_cache(double s, double s12, double t)
        {
            bool s_changed    = !are_equal(_cached_s, s, _cache_tolerance);
            bool s12_changed  = !are_equal(_cached_s12, s12, _cache_tolerance);
            bool t_changed    = !are_equal(_cached_t, t, _cache_tolerance);
            bool mR_changed   = !are_equal(_cached_mR, _kinematics->get_recoil_mass(), _cache_tolerance);

            bool need_update = _parameters_changed || s_changed || s12_changed || t_changed || mR_changed;

            if (need_update)
            {
                _cached_helicity_amplitudes.clear();

                // Total number of amplitudes
                int n = (helicity_independent()) ? 1 : _kinematics->N_amps();

                if (n == 1)
                {
                    // Arbitrary pick the first helicity set to evaluate, skip parity phase calculations
                    _cached_helicity_amplitudes.push_back(helicity_amplitude(_kinematics->helicities(0), s, s12, t));
                }
                // else
                // {
                //     // Otherwise we only calculate the first half amplitudes
                //     // with the rest related by the parity_phase
                //     for (int i = 0; i < n; i++)
                //     {
                //         complex hel_amp;
                //         if (i < n/2) hel_amp = helicity_amplitude(_kinematics->helicities(i), s, s12, t);
                //         else         hel_amp = parity_phase(i) * _cached_helicity_amplitudes[n - 1 - i];
                //         _cached_helicity_amplitudes.push_back(hel_amp); 
                //     }
                // };

                if (_cached_helicity_amplitudes.size() != n)
                {
                    warning(id()+"::update_cache", "Cached size doesn't match expected number of helicity amplitude!");
                };
            };

            // Update the cache info as well
            _cached_s  = s; _cached_s12 = s12; _cached_t = t;
            _cached_mR = _kinematics->get_recoil_mass();
            _parameters_changed = false;
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
                return error(id()+"::set_parameters", "Number of parameters passed not the expected size!", false);
            };
            return true;
        };
    };
};