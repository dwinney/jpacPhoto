// Parameterization of a resonant amplitude in the s-channel
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/baryon_resonance.hpp"

// Combined amplitude as a Breit-Wigner with the residue as the prodect of hadronic and photo-couplings
std::complex<double> jpacPhoto::baryon_resonance::helicity_amplitude(std::array<int, 4> helicities, double s, double t)
{
    int lam_i = 2 * helicities[0] - helicities[1];
    int lam_f = 2 * helicities[2] - helicities[3];

    // update save values of energies and angle
    update(helicities, s, t);

    std::complex<double> residue = 1.;
    residue  = photo_coupling(lam_i);
    residue *= hadronic_coupling(lam_f);
    residue *= threshold_factor(1.5);

    residue *= wigner_d_half(_resJ, lam_i, lam_f, _theta);
    residue /= (s + XI * _mRes * _gamRes - _mRes*_mRes);

    return residue;
};

// Ad-hoc threshold factor to kill the resonance at threshold
double jpacPhoto::baryon_resonance::threshold_factor(double beta)
{
    double result = pow((_s - _kinematics->sth()) / _s, beta);
    result /= pow((_mRes*_mRes - _kinematics->sth()) / (_mRes*_mRes), beta);

    return result;
};

// Photoexcitation helicity amplitude for the process gamma p -> R
std::complex<double> jpacPhoto::baryon_resonance::photo_coupling(int lam_i)
{
    // For spin-1/2 exchange no double flip
    if (_resJ == 1 && abs(lam_i) > 1) return 0.;

    // A_1/2 or A_3/2 depending on ratio R_photo
    double a;
    (std::abs(lam_i) == 1) ? (a = _photoR) : (a = sqrt(1. - _photoR * _photoR));

    // Electromagnetic decay width given by VMD assumption
    std::complex<double> emGamma = (_xBR * _gamRes) * pow(F_JPSI / M_JPSI, 2.);
    emGamma *= pow(XR * _pibar / _pfbar, double(2 * _lmin + 1)) * _pt;

    // Photo-coupling overall size of |A_1/2|^2 + |A_3/2|^2 is restriced from VMD
    std::complex<double> A_lam = emGamma * PI * _mRes * double(_resJ + 1) / (2. * M_PROTON * _pibar * _pibar);
    A_lam = sqrt(XR * A_lam);

    std::complex<double> result = sqrt(XR * _s) * _pibar / _mRes;
    result *= sqrt(XR * 8. * M_PROTON * _mRes / _kinematics->initial_momentum(_s));
    result *= A_lam * a;

    // FACTOR OF 4 PI SOMETIMES FACTORED OUT
    result *= sqrt(4. * PI * ALPHA);

    return result;
};

// Hadronic decay helicity amplitude for the R -> J/psi p process
std::complex<double> jpacPhoto::baryon_resonance::hadronic_coupling(int lam_f)
{
    // Hadronic coupling constant g, given in terms of branching ratio xBR
    std::complex<double> g;
    g  = 8. * PI * _xBR * _gamRes;
    g *=  _mRes * _mRes * double(_resJ + 1) / 6.;
    g /= pow(_pfbar, double(2 * _lmin + 1));
    g = sqrt(XR * g);

    std::complex<double> gpsi;
    gpsi = g * pow(_kinematics->final_momentum(_s), _lmin);

    (lam_f < 0) ? (gpsi *= double(_naturality)) : (gpsi *= 1.);

    return gpsi;
};
