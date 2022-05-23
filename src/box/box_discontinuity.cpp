// Discontinuity across unitarity cut for vector production via a one-loop box diagram
// Used as a container class in integration processes
// 
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "box/box_discontinuity.hpp"

// ---------------------------------------------------------------------------
// Evaluate the product of sub-amplitudes integrating over intermediate phase-space
double jpacPhoto::box_discontinuity::eval(double s)
{
    // Check the intermediate state matched
    if (_matchError) return 0.;

    // if below threshold return 0
    if (s < _initialAmp->_kinematics->sth()) return 0.;
    
    // make sure the external helicites get passed correctly
    int lam_gam = _external_helicities[0];
    int lam_tar = _external_helicities[1];
    int lam_vec = _external_helicities[2];
    int lam_rec = _external_helicities[3];
    
    auto dF = [&](const double * x)
    {
        double theta_gam = x[0];
        double phi_gam   = x[1];
            
        // calculate the sub-process momentum transfers
        double t_gam        = _initialAmp->_kinematics->t_man(s, theta_gam);
    
        double costheta_vec = cos(_external_theta) * cos(theta_gam) + sin(_external_theta) * sin(theta_gam) * cos(phi_gam);
        double theta_vec    = TMath::ACos(costheta_vec);
        double t_vec        = _finalAmp->_kinematics->t_man(s, theta_vec);

        // Sum over intermediate helicities 
        std::complex<double> result = 0.;
        for (int i = 0; i < 4*_jp_left[0]+2; i++)
        {
            int lam_meson  = _intermediate_helicities[i][2];
            int lam_baryon = _intermediate_helicities[i][3];

            std::complex<double> left, right;
            left  = _initialAmp->helicity_amplitude({lam_gam, lam_tar, lam_meson, lam_baryon}, s, t_gam);
            right = _finalAmp->helicity_amplitude(  {lam_vec, lam_rec, lam_meson, lam_baryon}, s, t_vec);
            result += left * right;
        };

        double jacobian = sin(theta_gam);
        return jacobian * real(result);
    };

    // Integrate over theta_gamma = [0, pi] and phi = [0, 2pi]
    double min[2] = {0., 0.};
    double max[2] = {PI, 2.*PI};
    
    ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kDEFAULT);
    ROOT::Math::Functor wF(dF, 2);
    ig.SetFunction(wF, 2);

    double result = ig.Integral(min, max);

    // Muliply by the two-body phase space
    double phase_space;
    phase_space  = std::real(_initialAmp->_kinematics->final_momentum(s));
    phase_space *= 2. / sqrt(s);

    return result * phase_space / (64. * PI*PI);
};