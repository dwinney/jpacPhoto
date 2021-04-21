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
        
        // If this is a single helicity ampltiude we can use the parity relation to only calculate half of the amplitudes
        if (!_isSum)
        {

            for (int i = 0; i < n/2; i++)
            {
                std::complex<double> amp_gamp = helicity_amplitude(_kinematics->_helicities[i], s, t);
                _cached_helicity_amplitude.push_back(amp_gamp);
            };

            for (int i = 0; i < n/2; i++)
            {
                std::complex<double> amp_gamp = _cached_helicity_amplitude[n/2 - 1 - i];
                double eta = double(parity_phase(_kinematics->_helicities[i]));
                _cached_helicity_amplitude.push_back( eta * amp_gamp);
            };
        }
        else
        {
            for (int i = 0; i < n; i++)
            {
                std::complex<double> amp_gamp = helicity_amplitude(_kinematics->_helicities[i], s, t);
                _cached_helicity_amplitude.push_back(amp_gamp);
            };
        };

        if (_cached_helicity_amplitude.size() != n)
        {
            std::cout << "Error! cache size not equal to expected number of helicity ampltiudes! Returning 0!\n";
            exit(1);
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

    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS61);
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
    if (alpha < 0 || alpha > 2 || std::abs(lam) > 2 || std::abs(lamp) > 2)
    {
        std::cout << "\nError! Invalid parameter passed to SDME. Returning 0!\n";
        return 0.;
    };

    if (!_kinematics->_photon) 
    {
        std::cout << "\nError! SDME only valid for photon in the initial state. Returning 0!\n";
        return 0.;
    };

    // If spin is too small return 0 automatically
    int j = _kinematics->_jp[0];
    if (j < 2 && (abs(lam) == 2 || abs(lamp) == 2)) return 0.;
    if (j < 1 && (abs(lam) >= 1 || abs(lamp) >= 1)) return 0.;

    // Phase and whether to conjugate at the end
    bool CONJ = false;
    double phase = 1.;

    // if first index smaller, switch them
    if (std::abs(lam) < std::abs(lamp))
    {
        int temp = lam;
        lam  = lamp;
        lamp = temp;

        CONJ = true;
    }

    // if first index is negative, flip to positive
    if (lam < 0)
    {
        lam  *= -1;
        lamp *= -1;

        phase *= pow(-1., double(lam - lamp));

        if (alpha == 2){phase *= -1.;};
    }

    // Normalization (sum over all amplitudes squared)
    double norm = probability_distribution(s, t);

    // k filters first index to be  0, 1, 2
    // l filters second index to be 0, 1, 2
    // m filters sign of second index
    int k, l, m;
    std::array<std::vector<int>, 2> iters = get_iters(j);
    
    int lamlamp = 10 * lam + abs(lamp);
    switch(lamlamp)
    {
        case  0: {k =  2; l =  2; break;};
        case 10: {k =  0; l =  2; break;};
        case 11: {k =  0; l =  0; break;};
        case 20: {k = -2; l =  2; break;};
        case 21: {k = -2; l =  0; break;};
        case 22: {k = -2; l = -2; break;};
        default: exit(0);
    };
    
    (lamp < 0) ? (m = 4 * j) : (m = 0);

    // Sum over the appropriate amplitude combinations
    std::complex<double> result = 0.;
    for (int i = 0; i < 8; i++)
    {
        int index;
        (alpha == 0) ? (index = iters[0][i]) : (index = iters[1][i]);
        std::complex<double> amp, amp_star, temp;
        amp      = _cached_helicity_amplitude[index + k];
        amp_star = _cached_helicity_amplitude[iters[0][i] + l + m];

        temp = real(amp * conj(amp_star));
        if (alpha == 2)
        {
            temp *= XI * double(_kinematics->_helicities[iters[0][i] + k][0]);
        }
        
        result += temp;
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
    double rho122 = real(SDME(1, 2, 2, s, t));
    double rho000 = real(SDME(0, 0, 0, s, t));
    double rho011 = real(SDME(0, 1, 1, s, t));
    double rho022 = real(SDME(0, 2, 2, s, t));

    return -(rho100 + 2. * rho111 + 2. * rho122) / (rho000 + 2. * rho011 + 2. * rho022);
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
    double rho12m2 = real(SDME(1, 2, -2, s, t));

    return 2. * rho11m1 - 2. * rho12m2 - rho100;
};
