// Implementation of a partial wave which takes in an amplitude and numerically computes
// the partial wave projection integral
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2026)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Instituto de Ciencias Nuclears (ICN),
//               Universidad Nacional Autonoma de México (UNAM)
// Email:        daniel.winney@nucleares.unam.mx
// ------------------------------------------------------------------------------

#ifndef PROJECTION_HPP
#define PROJECTION_HPP

#include "partial_wave.hpp"

namespace jpacPhoto
{
    class projected_amplitude; 

    partial_wave project(uint J, amplitude to_project)
    {
        partial_wave amp_ptr = std::make_shared<projected_amplitude>(key(), J, to_project);
        return amp_ptr;
    };

    class projected_amplitude : public raw_partial_wave
    {
        projected_amplitude(key k, uint J, amplitude to_project)
        : raw_partial_wave(k, to_project->get_kinematics(), J, "projected_amplitude"),
        _amplitude(to_project)
        {
            switch (to_project->native_helicity_frame())
            {
                case helicity_frame::S_CHANNEL:
                case helicity_frame::HELICITY_INDEPENDENT: _amplitude = to_project; break;
                case helicity_frame::T_CHANNEL: _amplitude = cross_to<helicity_frame::S_CHANNEL>(to_project); break;
                default: warning("project - Unknown native_helicity_frame()!");
            };
            set_N_pars(0);
        };

        // These are always assumed to be s-channel helicities so this is fixed
        helicity_frame native_helicity_frame()
        {
            bool hel_indep = _amplitude->native_helicity_frame() == helicity_frame::HELICITY_INDEPENDENT;
            return (hel_indep) ? helicity_frame::HELICITY_INDEPENDENT : helicity_frame::S_CHANNEL;
        };

        std::vector<quantum_numbers> allowed_mesons(){  return (_amplitude == nullptr) ? std::vector<quantum_numbers>() : _amplitude->allowed_mesons(); };
        std::vector<quantum_numbers> allowed_baryons(){ return (_amplitude == nullptr) ? std::vector<quantum_numbers>() : _amplitude->allowed_baryons(); };
        
        inline void allocate_parameters(int x)
        {
            warning("allocate_parameters - Partial-waves created using project() cannot change parameters! Change them in the amplitude being projected");
            return;
        };


        // Evaluate the partial-wave projection integral numerically and return only the s-dependent piece
        complex partial_wave(std::array<int,4> helicities, double s)
        {
            // Quick check that J is large enough for given helicities
            int lam  = 2 * helicities[0] - helicities[1]; // Photon - Target
            int lamp = 2 * helicities[2] - helicities[3]; // Meson  - Recoil

            if (_amplitude->native_helicity_frame() == helicity_frame::S_CHANNEL)
            {
                if ( std::abs(lam) > _J || std::abs(lamp) > _J ) return 0; 
                // if it is calculate the PWA integral
                auto F = [&](double theta)
                {
                    std::complex<double> integrand;
                    integrand  = sin(theta);
                    integrand *= wigner_d_half(_J, lam, lamp, theta);
                    integrand *= _amplitude->helicity_amplitude(helicities, s, _kinematics->t_man(s, theta));
                    return integrand/2;
                };
        
                return boost::math::quadrature::gauss_kronrod<double, 15>::integrate(F, 0., PI, 0., 1.E-6, NULL);
        }

            // if it is calculate the PWA integral
            auto F = [&](double theta)
            {
                std::complex<double> integrand;
                integrand  = sin(theta);
                integrand *= legendre(_J, theta);
                integrand *= _amplitude->helicity_amplitude(helicities, s, _kinematics->t_man(s, theta));
                return integrand/2;
            };

            return boost::math::quadrature::gauss_kronrod<double, 15>::integrate(F, 0., PI, 0., 1.E-6, NULL);
        };
        

        protected:

        // Pointer to a full amplitude and which they project
        amplitude _amplitude = nullptr;
    };
};

#endif