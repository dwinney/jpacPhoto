// Discontinuity across unitarity cut for vector production via a one-loop box diagram
// Used as a container class in integration processes
// 
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "box_discontinuity.hpp"

// ---------------------------------------------------------------------------
// Default behavior of the helicity amplitude which is to fix
std::complex<double> jpacPhoto::brute_force_discontinuity::helicity_amplitude(std::array<int,4> helicities, double s, double t)
{
    // Save external values
    _external_helicities = helicities; 
    _external_theta      = _kinematics->theta_s(s, t);

    // The one free parameter is the dispersion cutoff
    double s_cut = _params[0];

    if (s_cut < _kinematics->sth() + EPS)
    {
        std::cout<< "Warning! Dispersion integral cutoff less than threshold. Returning 0..." << std::endl;
        return 0.;
    };

    // Store the invariant energies to avoid having to pass them around 
    double sub = eval(s);
    auto F = [&](double sp)
    {
        std::complex<double> result = (eval(sp) - sub) / (sp - s - IEPS);
        return result;
    };

    std::complex<double> intpiece = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(F, _kinematics->sth() + EPS, s_cut, 0, 1.E-6, NULL);
    std::complex<double> logpiece = sub * (log(s_cut - s - IEPS) - log(_kinematics->sth() + EPS - s - IEPS));
    std::complex<double> result =  (intpiece + logpiece) / M_PI;

    return result;
};
// ---------------------------------------------------------------------------
// Evaluate the product of sub-amplitudes integrating over intermediate phase-space
double jpacPhoto::brute_force_discontinuity::eval(double s)
{
    // if below threshold return 0
    if (s < _kinematics->sth()) return 0.;
    
    // make sure the external helicites get passed correctly
    int lam_gam = _external_helicities[0];
    int lam_tar = _external_helicities[1];
    int lam_vec = _external_helicities[2];
    int lam_rec = _external_helicities[3];

    // Spins of the intermediate particles
    int mJ = _left_amp->_kinematics->get_meson_JP()[0];
    int bJ = _left_amp->_kinematics->get_baryon_JP()[0];
    int n_amps = (bJ+1)*(2*mJ+1); // Number of inermediate spin combinations

    auto dF = [&](const double * x)
    {
        // We're integrating over the angles of the left_amp
        double theta_gam = x[0];
        double phi_gam   = x[1];

        // Related to the other angle by subsequent rotations            
        double costheta_vec = cos(_external_theta) * cos(theta_gam) + sin(_external_theta) * sin(theta_gam) * cos(phi_gam);
        double theta_vec    = TMath::ACos(costheta_vec);

        // Calculate the sub-process momentum transfers
        double t_gam        =  _left_amp->_kinematics->t_man(s, theta_gam);
        double t_vec        = _right_amp->_kinematics->t_man(s, theta_vec);

        // Sum over intermediate helicities 
        double result = 0.;
        for (int i = 0; i < n_amps; i++)
        {
            int lam_meson  = _intermediate_helicities[i][2];
            int lam_baryon = _intermediate_helicities[i][3];

            std::complex<double> left, right;
            left    = _left_amp->helicity_amplitude( {lam_gam, lam_tar, lam_meson, lam_baryon}, s, t_gam);
            right   = _right_amp->helicity_amplitude({lam_vec, lam_rec, lam_meson, lam_baryon}, s, t_vec);

            debug(s, t_gam, t_vec);
            result += std::real(left * std::conj(right));
        };

        return sin(theta_gam) * result / (8. * PI);
    };

    // Integrate over theta_gamma = [0, pi] and phi = [0, 2pi]
    double min[2] = {0., 0.};
    double max[2] = {PI, 2.*PI};
    
    ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kDEFAULT);
    ROOT::Math::Functor wF(dF, 2);
    ig.SetFunction(wF, 2);

    double result = ig.Integral(min, max);

    return result * phase_space(s);
};

// ---------------------------------------------------------------------------
// Check the two amplitudes given are compatible with each other and with the desired overall process
void jpacPhoto::brute_force_discontinuity::compatibility_check(reaction_kinematics * kinem, amplitude * left, amplitude * right)
{
    // Grab all the quantum numbers for all three processes
    std::array<int,2> mjp_left  = left->_kinematics->get_meson_JP(),  bjp_left  = left->_kinematics->get_baryon_JP(); 
    std::array<int,2> mjp_right = right->_kinematics->get_meson_JP(), bjp_right = right->_kinematics->get_baryon_JP(); 
    if ( (mjp_left != mjp_right) || (bjp_left != bjp_right)) qn_error(left, right);

    // Now we should check all the masses involved
    double mX, mR, mB, mT;
    double mXL, mRL, mBL, mTL;
    double mXR, mRR, mBR, mTR;

    double mass_check_tolerance = 1.E-3;

    // Overall proces masses
    mX = kinem->get_meson_mass(); 
    mR = kinem->get_recoil_mass();
    mB = kinem->get_beam_mass();
    mT = kinem->get_target_mass();

    // Sub amplitude masses
    mXL = left->_kinematics->get_meson_mass();  mXR = right->_kinematics->get_meson_mass();
    mRL = left->_kinematics->get_recoil_mass(); mRR = right->_kinematics->get_recoil_mass();
    mBL = left->_kinematics->get_beam_mass();   mBR = right->_kinematics->get_beam_mass();
    mTL = left->_kinematics->get_target_mass(); mTR = right->_kinematics->get_target_mass();

    // Check the intermediate particles match
    if (   (abs(mXL - mXR) > mass_check_tolerance) || (abs(mRL - mRR) > mass_check_tolerance)
        || (abs(mBL - mB ) > mass_check_tolerance) || (abs(mTL - mT ) > mass_check_tolerance)
        || (abs(mBR - mX ) > mass_check_tolerance) || (abs(mTR - mR ) > mass_check_tolerance)
       )  mass_error(kinem, left, right);

    if (!_match_error) _intermediate_helicities = get_helicities(mjp_left[0], bjp_left[0]);
};

void jpacPhoto::brute_force_discontinuity::qn_error(amplitude * left, amplitude * right)
{
    // Grab all the quantum numbers for all three processes
    std::array<int,2> mjp_left  = left->_kinematics->get_meson_JP(),  bjp_left  = left->_kinematics->get_baryon_JP(); 
    std::array<int,2> mjp_right = right->_kinematics->get_meson_JP(), bjp_right = right->_kinematics->get_baryon_JP(); 

    std::cout << std::left;
    std::cout << std::setw(40) << "box_discontinuity: Error! Given sub-amplitude quantum numbers don't match!" << std::endl;
    std::cout << std::setw(20) << "Left  sub-process:" << "\t ("  << mjp_left[0]  << ", " << mjp_left[1] 
                                                          << ") & (" << bjp_left[0]  << ", " << bjp_left[1]  << ")" << std::endl;
    std::cout << std::setw(20) << "Right sub-process:" << "\t ("  << mjp_right[0] << ", " << mjp_right[1] 
                                                          << ") & (" << bjp_right[0] << ", " << bjp_right[1] << ")" << std::endl;
    std::cout << "Returning 0!" << std::endl;

    _match_error = true;
};

void jpacPhoto::brute_force_discontinuity::mass_error(reaction_kinematics * kinem, amplitude * left, amplitude * right)
{
    // Now we should check all the masses involved
    double mX, mR, mB, mT;
    double mXL, mRL, mBL, mTL;
    double mXR, mRR, mBR, mTR;

    // Overall proces masses
    mX = kinem->get_meson_mass(); 
    mR = kinem->get_recoil_mass();
    mB = kinem->get_beam_mass();
    mT = kinem->get_target_mass();

    // Sub amplitude masses
    mXL = left->_kinematics->get_meson_mass();  mXR = right->_kinematics->get_meson_mass();
    mRL = left->_kinematics->get_recoil_mass(); mRR = right->_kinematics->get_recoil_mass();
    mBL = left->_kinematics->get_beam_mass();   mBR = right->_kinematics->get_beam_mass();
    mTL = left->_kinematics->get_target_mass(); mTR = right->_kinematics->get_target_mass();

    std::cout << std::left;
    std::cout << std::setw(40) << "box_discontinuity: Error! Given amplitudes particle masses dont match!" << std::endl;
    std::cout << std::setw(20) << "Overall  process:  { " << mB  << " , " << mT  << " , " << mX  << " , " << mR  << "} \n";
    std::cout << std::setw(20) << "Left  sub-process: { " << mBL << " , " << mTL << " , " << mXL << " , " << mRL << "} \n";
    std::cout << std::setw(20) << "Right sub-process: { " << mBR << " , " << mTR << " , " << mXR << " , " << mRR << "} \n";
    std::cout << "Returning 0!" << std::endl;

    _match_error = true;
};