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
        bool same_kinem   = (a->_kinematics == b->_kinematics);

        // But also they must have helicites defined in the same frame
        // TODO: add helicity crossing relations so this is not a necessity
        bool same_frame   = (a->native_helicity_frame() == b->native_helicity_frame());

        std::string error_msg = "Attempted to add incompatible amplitudes: " + a->id() + " and " + b->id() + "!";
        if (!same_frame) return error(error_msg + " (Helicites defined in different frames!)", false);
        if (!same_kinem) return error(error_msg + " (Contain different kinematics objects!)", false);

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
        return std::make_shared<raw_amplitude>(key(), from_a, id);
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
        bool same_kinem   = (new_amp->_kinematics == _kinematics);

        // But also they must have helicites defined in the same frame
        // TODO: add helicity crossing relations so this is not a necessity
        bool same_frame   = (new_amp->native_helicity_frame() == native_helicity_frame());

        std::string error_msg = "Attempted to add incompatible amplitudes: " + id() + " and " + new_amp->id() + "!";
        if (!same_frame) return error(error_msg + " (Helicites defined in different frames!)", false);
        if (!same_kinem) return error(error_msg + " (Contain different kinematics objects!)", false);
        return true;
    };

    // Add a new amplitude to the list
    void raw_amplitude::add(amplitude new_amp)
    {
        // If our list is empty, capture the first entry we encounter
        if (_subamplitudes.size() == 0)
        {
            _kinematics = new_amp->_kinematics; _N_pars     = new_amp->N_pars();
            _subamplitudes.push_back(new_amp);
            return;
        };

        // Each subsequent entry gets compared to the first
        if (is_compatible(new_amp))
        {
            // If new_amp is itself a sum of amps,
            // reach in and pull out the components 
            if (new_amp->is_sum()){ add(new_amp->_subamplitudes); return; };

            _subamplitudes.push_back(new_amp); _N_pars += new_amp->N_pars();
            return;
        };
    };

    // Add a whole vector, this just loops over each entry
    void raw_amplitude::add(std::vector<amplitude> new_amps)
    {
        for (auto amp : new_amps) add(amp);
    };

    // ------------------------------------------------------------------------------
    // Updators for the different caches we have

    // Save current kinematic quantities from _kinematics
    void raw_amplitude::store(std::array<int,4> helicities, double s, double t)
    {
        _lamB = helicities[0]; _lamT = helicities[1];
        _lamX = helicities[2]; _lamR = helicities[3];
    
        _s = s; _t = t; _theta = _kinematics->theta_s(s,t);

        _mB = _kinematics->get_beam_mass();  _mT = _kinematics->get_target_mass();
        _mX = _kinematics->get_meson_mass(); _mR = _kinematics->get_recoil_mass();
        return;
    };

    // Check for any changes within the set tolerance 
    // And recalculate helicity amplitudes if needed
    void raw_amplitude::update_cache(double s, double t)
    {
        bool st_changed   = !are_equal(_cached_s, s, _cache_tolerance) || !are_equal(_cached_t, t, _cache_tolerance);
        bool mass_changed = !are_equal(_cached_mX, _kinematics->get_meson_mass(),  _cache_tolerance) && !are_equal(_cached_mR, _kinematics->get_recoil_mass(), _cache_tolerance);

        bool need_update = _parameters_changed || st_changed || mass_changed;

        if (need_update)
        {
            _cached_helicity_amplitudes.clear();

            // Total number of amplitudes
            int n = (native_helicity_frame() == HELICITY_INDEPENDENT) ? 1 : _kinematics->N_amps();

            if (n == 1)
            {
                // Arbitrary pick the first helicity set to evvaluate, skip parity phase calculations
                _cached_helicity_amplitudes.push_back(helicity_amplitude(_kinematics->helicities(0), s, t));
            }
            else
            {
                // Otherwise we only calculate the first half amplitudes
                // with the rest related by the parity_phase
                for (int i = 0; i < n; i++)
                {
                    complex hel_amp;
                    if (i < n/2) hel_amp = helicity_amplitude(_kinematics->helicities(i), s, t);
                    else         hel_amp = parity_phase(i) * _cached_helicity_amplitudes[n - 1 - i];
                    _cached_helicity_amplitudes.push_back(hel_amp); 
                }
            };

            if (_cached_helicity_amplitudes.size() != n) warning(id()+"::update_cache", "Cached size doesn't match expected number of helicity amplitude!");
        };

        // Update the cache info as well
        _cached_mX = _kinematics->get_meson_mass(); _cached_mR = _kinematics->get_recoil_mass();
        _cached_s  = s; _cached_t = t;
        _parameters_changed = false;
    };

    // ------------------------------------------------------------------------------
    // Phase properties of helicity amplitudes

    double raw_amplitude::parity_phase(std::array<int,4> helicities){ return _kinematics->parity_phase(helicities, this->native_helicity_frame()); };
    double raw_amplitude::parity_phase(int i){ return _kinematics->parity_phase(i, this->native_helicity_frame()); };

    // ------------------------------------------------------------------------------
    // Parameter handling

    // Ability to change number of expected parameters from inside a class
    void raw_amplitude::set_N_pars(int N){ _N_pars = N; };

    // Simple check that a given vector is of the expected size
    bool raw_amplitude::correct_size(std::vector<double> pars)
    {
        if (pars.size() != _N_pars) return error(id()+"::set_parameters", "Number of parameters passed not the expected size!", false);
        return true;
    };

    // ------------------------------------------------------------------------------
    // Compatability checks

    // At construction, check the kinematics with set quantum numbers is compatable
    void raw_amplitude::check_QNs(kinematics kinem)
    {
        // Get all the allowed JP's from the amplitude
        std::vector<quantum_numbers> allowed_mesons, allowed_baryons;
        allowed_mesons  = this->allowed_mesons(); allowed_baryons = this->allowed_baryons();

        // Grab the requested JPs from the raw_kinematics
        quantum_numbers requested_meson, requested_baryon;
        requested_meson  = kinem->get_meson(); requested_baryon = kinem->get_baryon();

        // Check if they are allowed
        bool meson_fails, baryon_fails;
        meson_fails  = std::find(allowed_mesons.begin(),  allowed_mesons.end(),  requested_meson)  == allowed_mesons.end();
        baryon_fails = std::find(allowed_baryons.begin(), allowed_baryons.end(), requested_baryon) == allowed_baryons.end();
    
        auto requested_meson_JP = kinem->get_meson_JP(); auto requested_baryon_JP = kinem->get_meson_JP();

        if (meson_fails)  warning(id()+"::check_QNs", "Requested meson quantum numbers (J=" + std::to_string(requested_meson_JP[0]) + ", P=" + std::to_string(requested_meson_JP[1])+") not available!");
        if (baryon_fails) warning(id()+"::check_QNs", "Requested baryon quantum numbers (J=" + std::to_string(requested_baryon_JP[0]) + "/2, P=" + std::to_string(requested_baryon_JP[1])+") not available!");
    };

    // ------------------------------------------------------------------------------
    // OBSERVABLES
    
    // -------------------------------------------------------
    // Cross sections 

    // Square of the spin averaged amplitude squared
    double raw_amplitude::probability_distribution(double s, double t)
    {
        // Check we have the right amplitudes cached
        double sum = 0.; auto cache = get_cache(s, t);
        for (auto amp : cache) sum += std::norm(amp);
        return sum;
    };

    // Differential cross section dsigma / dt
    // in NANOBARN
    double raw_amplitude::differential_xsection(double s, double t)
    {
        if (s < _kinematics->sth()) return 0.;

        double sum = probability_distribution(s, t);
        double norm = 64. * PI * s * pow(_kinematics->initial_momentum(s), 2.) * (2.56819E-6); // Convert from GeV^-2 -> nb

        // Average over initial helicities
        if (native_helicity_frame() !=  HELICITY_INDEPENDENT) norm *= 4*(_kinematics->is_photon()) + 6*(!_kinematics->is_photon());

        return sum / norm;
    };

    // Integrated total cross-section
    // IN NANOBARN
    double raw_amplitude::integrated_xsection(double s)
    {
        if (s < _kinematics->sth()) return 0.;

        ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS61);
        ROOT::Math::Functor1D wF([&](double t){ return differential_xsection(s, t); });
        ig.SetFunction(wF);

        double t_min = _kinematics->t_min(s); double t_max = _kinematics->t_max(s);
        return ig.Integral(t_max, t_min);
    };


    // Differential cross section dsigma_perp/para / dt
    // in NANOBARN
    double raw_amplitude::polarized_differential_xsection(int perp_or_para, double s, double t)
    {
        if (s < _kinematics->sth()) return 0.;
        if (abs(perp_or_para) != 1) return std::nan("");
        
        // Sum first half of amplitudes which are lam_gamma = +1
        auto cache = get_cache(s, t); int n = cache.size()/2;

        double sum = 0.;
        for (int i = 0; i < n; i++) sum += std::norm(cache[i] + perp_or_para * cache[n + i]);
        double norm = 64. * PI * s * pow(_kinematics->initial_momentum(s), 2.) * (2.56819E-6); // Convert from GeV^-2 -> nb

        // Average over initial helicities
        if (native_helicity_frame() !=  HELICITY_INDEPENDENT) { norm *= 4*(_kinematics->is_photon()) + 6*(!_kinematics->is_photon()); }

        return sum / norm;
    };

    // -------------------------------------------------------
    // Asymmetries 

    // Polarization asymmetry between beam and recoil baryon
    double raw_amplitude::K_LL(double s, double t)
    {
        double sum = 0; 
        auto cache = get_cache(s, t); int n = cache.size();
        for (int i = 0; i < n; i++)
        {
            auto hel = _kinematics->helicities(i);
            int eta  = (1 - hel[0]*hel[3] / abs(hel[0]*hel[3])) / 2;
            sum += pow(-1, eta) * std::norm(cache[i]);
        }
        return sum / probability_distribution(s, t);
    };

    // Polarization asymmetry between beam and target proton
    double raw_amplitude::A_LL(double s, double t)
    {
        double sum = 0; 
        auto cache = get_cache(s, t); int n = cache.size();
        for (int i = 0; i < n; i++)
        {
            auto hel = _kinematics->helicities(i);
            int eta  = (1 - hel[0]*hel[1] / abs(hel[0]*hel[1])) / 2;
            sum += pow(-1, eta) * std::norm(cache[i]);
        }
        return sum / probability_distribution(s, t);
    };

    // Integrated beam asymmetry sigma_4pi
    double raw_amplitude::beam_asymmetry_4pi(double s, double t)
    {
        double sum = 0;
        auto cache = get_cache(s, t); int n = cache.size()/2;
        for (int i = 0; i < n; i++) sum += 2.*std::real(cache[i]*conj(cache[i + n]));
        return sum / probability_distribution(s, t);
    };

    // -------------------------------------------------------
    // baryon SDME's

    // Baryon SDME
    complex raw_amplitude::bSDME(unsigned int alpha, int lam, int lamp, double s, double t)
    {
        if (alpha > 2) return error("amplitude::bSDME", "Invalid SDME (alpha = " + std::to_string(alpha) + ") requested!", std::nan(""));

        int J = _kinematics->get_baryon_JP()[0];
        if (std::abs(lam)  > J) return error("amplitude::bSDME", "Invalid SDME (lam = "  + std::to_string(lam)  + ") requested!", std::nan(""));
        if (std::abs(lamp) > J) return error("amplitude::bSDME", "Invalid SDME (lam' = " + std::to_string(lamp) + ") requested!", std::nan(""));

        // Check lam > lamp and lam > 0
        bool CONJ = false; int phase = 1;

        if (std::abs(lam) < abs(lamp))
        {
            int temp = lam; lam = lamp; lamp = temp; // Swap them
            CONJ = true; // Conjugate at the end
        };

        if (lam < 0)
        {
            lam *= -1; lamp *= -1; // Negate both
            phase = pow(-1, (lam - lamp)/ 2 + (alpha == 2)); // Add phase at end
        };

        auto cache = get_cache(s, t); int n = cache.size()/2;
        complex sum = 0;
        for (int i = 0; i < n; i++)
        {
            auto hi = _kinematics->helicities(i); 
            if (hi[3] != lam) continue;

            for (int j = 0; j < n; j++)
            {
                // First n amplitudes have lg = +1 second half have lg = -1
                auto hj = _kinematics->helicities(j); 
                if (hj[3] != lamp || hi[1] != hj[1]) continue;

                switch (alpha)
                {
                    case 0: sum += cache[i] * conj(cache[j])   + cache[i+n]*conj(cache[j+n]); break;
                    case 1: sum += cache[i] * conj(cache[j+n]) + cache[i+n]*conj(cache[j]);   break;
                    case 2: sum += cache[i] * conj(cache[j+n]) - cache[i+n]*conj(cache[j]);   break;
                };
            };
        };
        
        if (alpha == 2) sum *= I;
        if (CONJ) sum = std::conj(sum);
        return phase * sum / probability_distribution(s, t);        
    };  

    complex raw_amplitude::rotated_bSDME(unsigned int alpha, int lam, int lamp, double s, double t, double theta)
    {
        // First get the spin of the particle we're rotating
        int J = _kinematics->get_baryon_JP()[0];

        complex result = 0.;
        for (int m = -J; m <= J; m += 2)
        {
            for (int mp = -J; mp <= J; mp += 2)
            {
                result += wigner_d_half(J, lam, m, -theta) * bSDME(alpha, m, mp, s, t) * wigner_d_half(J, mp, lamp, theta);
            };
        };

        return result;
    };

    // Calculate the Helicity frame SDME
    complex raw_amplitude::bSDME_H(unsigned int alpha, int lam, int lamp, double s, double t)
    {
        switch (this->native_helicity_frame())
        {
            case S_CHANNEL: return bSDME(alpha, lam, lamp, s, t);
            case T_CHANNEL: return rotated_bSDME(alpha, lam, lamp, s, t, -_kinematics->bH_to_GJ_angle(s, t));
            case U_CHANNEL: return error("bSDME_H", "Rotations from u-channel CM frame to Helicty frame not yet implemented... Returning 0.", std::nan(""));
        };

        return std::nan("");
    };

    // Calculate the Gottfried-Jackson SDME
    complex raw_amplitude::bSDME_GJ(unsigned int alpha, int lam, int lamp, double s, double t)
    {
        switch (this->native_helicity_frame())
        {
            case S_CHANNEL: return rotated_bSDME(alpha, lam, lamp, s, t, _kinematics->bH_to_GJ_angle(s, t));
            case T_CHANNEL: return bSDME(alpha, lam, lamp, s, t);
            case U_CHANNEL: return error("bSDME_GJ", "Rotations from u-channel CM frame to Gottfried-Jackson frame not yet implemented... Returning 0.", std::nan(""));
            default: return std::nan("");
        };

        return std::nan("");
    };

    // -------------------------------------------------------
    // Meson SDME's

    // "raw" SDME
    complex raw_amplitude::mSDME(unsigned int alpha, int lam, int lamp, double s, double t)
    {
        if (alpha > 2) return error("amplitude::mSDME", "Invalid SDME (alpha = " + std::to_string(alpha) + ") requested!", std::nan(""));

        int J = _kinematics->get_meson_JP()[0];
        if (std::abs(lam)  > J) return error("amplitude::mSDME", "Invalid SDME (lam = "  + std::to_string(lam)  + ") requested!", std::nan(""));
        if (std::abs(lamp) > J) return error("amplitude::mSDME", "Invalid SDME (lam' = " + std::to_string(lamp) + ") requested!", std::nan(""));

        // Check lam > lamp and lam > 0
        bool CONJ = false; int phase = 1;

        if (std::abs(lam) < abs(lamp))
        {
            int temp = lam; lam = lamp; lamp = temp; // Swap them
            CONJ = true; // Conjugate at the end
        };

        if (lam < 0)
        {
            lam *= -1; lamp *= -1; // Negate both
            phase = pow(-1, (lam - lamp) + (alpha == 2)); // Add phase at end
        };

        auto cache = get_cache(s, t); int n = cache.size()/2;
        complex sum = 0;
        for (int i = 0; i < n; i++)
        {
            auto hi = _kinematics->helicities(i); 
            if (hi[2] != lam) continue;

            for (int j = 0; j < n; j++)
            {
                // First n amplitudes have lg = +1 second half have lg = -1
                auto hj = _kinematics->helicities(j); 
                if (hj[2] != lamp || hi[1] != hj[1]) continue;

                switch (alpha)
                {
                    case 0: sum += cache[i] * conj(cache[j])   + cache[i+n]*conj(cache[j+n]); break;
                    case 1: sum += cache[i] * conj(cache[j+n]) + cache[i+n]*conj(cache[j]);   break;
                    case 2: sum += cache[i] * conj(cache[j+n]) - cache[i+n]*conj(cache[j]);   break;
                };
            };
        };
        
        if (alpha == 2) sum *= I;
        if (CONJ) sum = std::conj(sum);
        return phase * sum / probability_distribution(s, t);        
    };  

    complex raw_amplitude::rotated_mSDME(unsigned int alpha, int lam, int lamp, double s, double t, double theta)
    {
        // First get the spin of the particle we're rotating
        int J = _kinematics->get_meson_JP()[0];

        // Sum over all the rest
        complex result = 0.;
        for (int m = -J; m <= J; m++)
        {
            for (int mp = -J; mp <= J; mp++)
            {
                result += wigner_d_int(J, lam, m, -theta) * mSDME(alpha, m, mp, s, t) * wigner_d_int(J, mp, lamp, theta);
            };
        };
        return result;
    };

    // Calculate the Helicity frame SDME
    complex raw_amplitude::mSDME_H(unsigned int alpha, int lam, int lamp, double s, double t)
    {
        switch (this->native_helicity_frame())
        {
            case S_CHANNEL: return mSDME(alpha, lam, lamp, s, t);
            case T_CHANNEL: return rotated_mSDME(alpha, lam, lamp, s, t, -_kinematics->mH_to_GJ_angle(s, t));
            case U_CHANNEL: return error("mSDME_H", "Rotations from u-channel CM frame to Helicty frame not yet implemented... Returning 0.", std::nan(""));
            default: return std::nan("");
        };
        return std::nan("");
    };

    // Calculate the Gottfried-Jackson SDME
    complex raw_amplitude::mSDME_GJ(unsigned int alpha, int lam, int lamp, double s, double t)
    {
        switch (this->native_helicity_frame())
        {
            case S_CHANNEL: return rotated_mSDME(alpha, lam, lamp, s, t, _kinematics->mH_to_GJ_angle(s, t));
            case T_CHANNEL: return mSDME(alpha, lam, lamp, s, t);
            case U_CHANNEL: return error("mSDME_GJ", "Rotations from u-channel CM frame to Gottfried-Jackson frame not yet implemented... Returning 0.", std::nan(""));
            default: return std::nan("");
        };
        return std::nan("");
    };
};