// Extension of the raw_amplitude which supports individual partial waves in the s-channel
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#include "partial_wave.hpp"

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    // Add two partial waves together
    
    amplitude operator+(partial_wave a, partial_wave b)
    {
        amplitude wa = std::static_pointer_cast<raw_amplitude>(a);
        amplitude wb = std::static_pointer_cast<raw_amplitude>(b);
        
        return wa + wb;
    };

    amplitude project(int J, amplitude to_project)
    {
        amplitude amp_ptr = std::make_shared<raw_partial_wave>(key(), J, to_project);
        return amp_ptr;
    };

    // ---------------------------------------------------------------------------
    // Member methods of partial wave amplitudes

    // Contribution to the full amplitude from this partial partial wave
    // We multiply by the appropriate dfunct and normalization factors
    complex raw_partial_wave::helicity_amplitude(std::array<int,4> helicities, double s, double t)
    {
        // s-channel scattering angle
        double theta = _kinematics->theta_s(s, t);

        switch (_amplitude->native_helicity_frame())
        {
            case helicity_frame::S_CHANNEL: 
            {
                // Net helicities
                int lam  = 2 * helicities[0] - helicities[1]; // Photon - Target
                int lamp = 2 * helicities[2] - helicities[3]; // Meson  - Recoil
                
                return (_J + 1) * wigner_d_half(_J, lam, lamp, theta) * partial_wave(helicities, s);
            }
            case helicity_frame::HELICITY_INDEPENDENT:
            {
                return (2*_J+1) * legendre(_J, cos(theta)) * partial_wave(s);
            };
            default: return NaN<complex>();
        };
    };


    // Evaluate the partial-wave projection integral numerically and return only the s-dependent piece
    complex raw_partial_wave::partial_wave(std::array<int,4> helicities, double s)
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
};