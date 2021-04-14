// Abstract class for an amplitude. Used so we can easily build observables
// as the incoherent sum of amplitudes in s, t, and u channels.
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/amplitude.hpp"

// ---------------------------------------------------------------------------

void jpacPhoto::amplitude::check_cache(double s, double t)
{
    // check if saved version its the one we want
    if (  (abs(_cached_s - s) < 0.00001) && 
          (abs(_cached_t - t) < 0.00001) &&
          (abs(_cached_mX2 - _kinematics->_mX2) < 0.00001) // important to make sure the value of mX2 hasnt chanced since last time
       )
    {
        return; // do nothing
    }
    else // save a new set
    {
        _cached_helicity_amplitude.clear();

        int n = _kinematics->_nAmps;

        // Save the first half for lam_gam = +1 
        for (int i = 0; i < n/2; i++)
        {
            std::complex<double> amp_gamp = helicity_amplitude(_kinematics->_helicities[i], s, t);
            _cached_helicity_amplitude.push_back(amp_gamp);
        };

        // Rest given by parity 
        for (int i = n/2; i <= n; i++)
        {
            int phase = _kinematics->_jp[1] * pow(-1, _kinematics->_jp[0]);
            _cached_helicity_amplitude.push_back(double(phase) *  _cached_helicity_amplitude[i - n/2]);
        };

        // update cache info
        _cached_mX2 = _kinematics->_mX2; _cached_s = s; _cached_t = t;
    }

    return;
};

// ---------------------------------------------------------------------------
// Square of the spin averaged amplitude squared
double jpacPhoto::amplitude::probability_distribution(double s, double t)
{
    // Check we have the right amplitudes cached
    check_cache(s, t);

    double sum = 0.;
    for (int i = 0; i < _kinematics->_nAmps; i++)
    {
        std::complex<double> amp_i = _cached_helicity_amplitude[i];
        sum += std::real(amp_i * conj(amp_i));
    }

    return sum;
};

// ---------------------------------------------------------------------------
// Differential cross section dsigma / dt
// in NANOBARN
double jpacPhoto::amplitude::differential_xsection(double s, double t)
{
    double sum = probability_distribution(s, t);

    double norm = 1.;
    norm /= 64. * PI * s;
    norm /= real(pow(_kinematics->_initial_state->momentum(s), 2.));
    norm /= (2.56819E-6); // Convert from GeV^-2 -> nb
    norm /= 4.; // Average over initial state helicites

    return norm * sum;
};

// ---------------------------------------------------------------------------
// Inegrated total cross-section
// IN NANOBARN
double jpacPhoto::amplitude::integrated_xsection(double s)
{
    auto F = [&](double t)
    {
        return differential_xsection(s, t);
    };

    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE);
    ROOT::Math::Functor1D wF(F);
    ig.SetFunction(wF);

    double t_min = _kinematics->t_man(s, 0.);
    double t_max = _kinematics->t_man(s, PI);

    return ig.Integral(t_max, t_min);
};

// ---------------------------------------------------------------------------
// Polarizatiopn asymmetry between beam and recoil proton
double jpacPhoto::amplitude::K_LL(double s, double t)
{
    // Check we have the right amplitudes cached
    check_cache(s, t);

    double sigmapp = 0., sigmapm = 0.;
    for (int i = 0; i < 6; i++)
    {
        std::complex<double> squarepp, squarepm;

        // Amplitudes with lam_gam = + and lam_recoil = +
        squarepp  = _cached_helicity_amplitude[2*i+1];
        squarepp *= conj(squarepp);
        sigmapp  += real(squarepp);

        // Amplitudes with lam_gam = + and lam_recoil = -
        squarepm  = _cached_helicity_amplitude[2*i];
        squarepm *= conj(squarepm);
        sigmapm  += real(squarepm);
    }

    return (sigmapp - sigmapm) / (sigmapp + sigmapm);
}

// ---------------------------------------------------------------------------
// Polarization asymmetry between beam and target proton
double jpacPhoto::amplitude::A_LL(double s, double t)
{
    // Check we have the right amplitudes cached
    check_cache(s, t);

    double sigmapp = 0., sigmapm = 0.;
    for (int i = 0; i < 6; i++)
    {
        std::complex<double> squarepp, squarepm;

        // Amplitudes with lam_gam = + and lam_targ = +
        squarepp  = _cached_helicity_amplitude[i+6];
        squarepp *= conj(squarepp);
        sigmapp  += real(squarepp);

        // Amplitudes with lam_gam = + and lam_targ = -
        squarepm  = _cached_helicity_amplitude[i];
        squarepm *= conj(squarepm);
        sigmapm  += real(squarepm);
    }

    return (sigmapp - sigmapm) / (sigmapp + sigmapm);
}

// ---------------------------------------------------------------------------
// Photon spin-density matrix elements
std::complex<double> jpacPhoto::amplitude::SDME(int alpha, int lam, int lamp, double s, double t)
{
    if (alpha < 0 || alpha > 2 || std::abs(lam) > 1 || std::abs(lamp) > 1)
    {
        std::cout << "\nError! Invalid parameter passed to SDME. Returning 0!\n";
        return 0.;
    };

    if (_kinematics->_mB > 1.E-5) 
    {
        std::cout << "\nError! SDME only valid for photon in the initial state. Returning 0!\n";
        return 0.;
    };

    // Phase and whether to conjugate at the end
    bool CONJ = false;
    double phase = 1.;

    // if first index smaller, switch them
    if (std::abs(lam) < std::abs(lamp))
    {
        int temp = lam;
        lam = lamp;
        lamp = temp;

        CONJ = true;
    }

    // if first index is negative, flip to positive
    if (lam < 0)
    {
        lam *= -1;
        lamp *= -1;

        phase *= pow(-1., double(lam - lamp));
    }

    // Normalization (sum over all amplitudes squared)
    double norm = probability_distribution(s, t);

    // These are the indexes of the amplitudes in reaction_kinematics that have
    // lambda_V = +1
    std::vector<int> pos_iters = {0, 1, 6, 7, 12, 13, 18, 19};
    std::vector<int> neg_iters = {12, 13, 18, 19, 0, 1, 6, 7};

    // j and k filter the right helicity combinations for 00, 1-1, 10, 11
    int j, k;
    (lam == 0) ? (k = 2) : (k = 0);

    switch (lamp)
    {
        case -1: { j = 4; break; }
        case 0:  { j = 2; break; }
        case 1:  { j = 0; break; }
        default:
        {
            std::cout << "\nSDME: Invalid parameter. J/Psi helicity projection lamp = 0 or 1! Returning zero.";
            return 0.;;
        }
    }

    // Sum over the appropriate amplitude combinations
    std::complex<double> result = 0.;
    for (int i = 0; i < 8; i++)
    {
        int index;
        (alpha == 0) ? (index = pos_iters[i]) : (index = neg_iters[i]);

        std::complex<double> amp_i, amp_j;
        amp_i = _cached_helicity_amplitude[index + k];
        amp_j = _cached_helicity_amplitude[pos_iters[i] + j];

        (alpha == 2) ? (amp_j *= XI * double(_kinematics->_helicities[pos_iters[i] + j][0])) : (amp_j *= XR);

        result += real(amp_i * conj(amp_j));
    }

    if (CONJ == true)
    {
        result = conj(result);
    }

    result /= norm;
    result *= phase;

    return result;
};

// ---------------------------------------------------------------------------
// Integrated beam asymmetry sigma_4pi
double jpacPhoto::amplitude::beam_asymmetry_4pi(double s, double t)
{
    double rho100 = real(SDME(1, 0, 0, s, t));
    double rho111 = real(SDME(1, 1, 1, s, t));
    double rho000 = real(SDME(0, 0, 0, s, t));
    double rho011 = real(SDME(0, 1, 1, s, t));

    return -(rho100 + 2. * rho111) / (rho000 + 2. * rho011);
};
// ---------------------------------------------------------------------------
// Beam asymmetry along y axis sigma_y 
double jpacPhoto::amplitude::beam_asymmetry_y(double s, double t)
{
    double rho111  = real(SDME(1, 1,  1, s, t));
    double rho11m1 = real(SDME(1, 1, -1, s, t));
    double rho011  = real(SDME(0, 1,  1, s, t));
    double rho01m1 = real(SDME(0, 1, -1, s, t));

    return (rho111 + rho11m1) / (rho011 + rho01m1);
};

// ---------------------------------------------------------------------------
// Parity asymmetry P_sigma
double jpacPhoto::amplitude::parity_asymmetry(double s, double t)
{
    double rho100  = real(SDME(1, 0,  0, s, t));
    double rho11m1 = real(SDME(1, 1, -1, s, t));

    return 2. * rho11m1 - rho100;
};
