// // Class to contain all relevant kinematic quantities. 
// (e.g. angles, momenta, helicity info)
//
// Full kinematics of the process are entirely determiend by the masses of all particles:
// mB - beam mass 
// mT - target (baryon) mass
// mX - produced meson mass
// mR - recoil baryon mass
//
// As well as quantum numbers of the final states
// All of these are allowed to float with which combinations are allowed handled by individual amplitudes
// 
// --------------------------------------------------------------------------------
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// --------------------------------------------------------------------------------

#include "kinematics.hpp"

namespace jpacPhoto
{
    namespace one_meson
    {
        // ---------------------------------------------------------------------------
        // Quantum number handling

        particle raw_kinematics::get_meson()
        {
            auto JP = get_meson_JP();
            int  qns =  10 * JP[0] + (JP[1] == 1);
            switch(qns)
            {
                case (11): return AXIALVECTOR;
                case (10): return VECTOR;
                case ( 0): return PSEUDOSCALAR;
                case ( 1): return SCALAR;
            };  
            return PARTICLE_ERROR;
        };

        void raw_kinematics::set_meson_JP(particle x)
        {
            int J, P;
            switch (x)
            {
                case ( AXIALVECTOR): J = 1, P = +1; break;
                case (      VECTOR): J = 1, P = -1; break;
                case (      SCALAR): J = 0, P = +1; break;
                case (PSEUDOSCALAR): J = 0, P = -1; break;
                default: return;
            };

            set_meson_JP(J, P);
        };

        particle raw_kinematics::get_baryon()
        {
            auto JP = get_baryon_JP();
            int  qns =  10 * JP[0] + (JP[1] == 1);
            switch(qns)
            {
                case (11): return HALFPLUS;
                case (10): return HALFMINUS;
                case (31): return THREEPLUS;
                case (30): return THREEMINUS;
            };  
            return PARTICLE_ERROR;
        };

        void raw_kinematics::set_baryon_JP(particle x)
        {
            int J, P;
            switch (x)
            {
                case (  HALFPLUS): J = 1, P = +1; break;
                case ( HALFMINUS): J = 1, P = -1; break;
                case ( THREEPLUS): J = 3, P = +1; break;
                case (THREEMINUS): J = 3, P = -1; break;
                default: return;
            };

            set_baryon_JP(J, P);
        };


        // ---------------------------------------------------------------------------
        // Find a set of helicities given an index

        std::array<int, 4> raw_kinematics::helicities(int i)
        {
            if (i < 0 || i >= _nAmps) 
            {
                fatal("kinematics", "Can't find helicities with index " + std::to_string(i) + "!");
            }
            return _helicities[i];
        };

        // ---------------------------------------------------------------------------
        // Modulus of three-momenta (in the s-channel)

        double raw_kinematics::initial_momentum(double s)
        {
            return sqrt(Kallen(s, _masses._mB2, _masses._mT2)) / sqrt(4.*s);
        };

        double raw_kinematics::final_momentum(double s)
        {
            return sqrt(Kallen(s, _masses._mR2, _masses._mX2)) / sqrt(4.*s);
        };

        // ------------------------------------------------------------------------------
        // Energies of all particles in s-channel CM frame

        double raw_kinematics::beam_energy(double s)
        {
            return (s - _masses._mT2 + _masses._mB2) / sqrt(4.*s);
        };

        double raw_kinematics::target_energy(double s)
        {
            return (s + _masses._mT2 - _masses._mB2) / sqrt(4.*s);
        };

        double raw_kinematics::meson_energy(double s)
        {
            return (s - _masses._mR2 + _masses._mX2) / sqrt(4.*s);
        };

        double raw_kinematics::recoil_energy(double s)
        {
            return (s + _masses._mR2 - _masses._mX2) / sqrt(4.*s);
        };

        // ---------------------------------------------------------------------------
        // s-channel scattering length from invariants

        double raw_kinematics::z_s(double s, double t)
        {
            complex qdotqp = initial_momentum(s) * final_momentum(s);
            complex E1E3   = beam_energy(s) * meson_energy(s);

            double result = t - _masses._mX2 - _masses._mB2 + 2.*real(E1E3);
            result /= 2. * real(qdotqp);

            return result;
        };

        // Scattering angle in the s-channel
        // Use TMath::ACos instead of std::acos because its safer at the end points
        double raw_kinematics::theta_s(double s, double t)
        {
            double zs = z_s(s, t);
            if (std::abs(zs - 1) < 1.E-5){zs =  1.;}
            if (std::abs(zs + 1) < 1.E-5){zs = -1.;}
            return TMath::ACos( zs );
        };

        // ---------------------------------------------------------------------------
        // Other Mandelstam invariants and scattering angles

        double raw_kinematics::t_man(double s, double theta)
        {
            complex qdotqp = initial_momentum(s) * final_momentum(s);
            complex E1E3   = beam_energy(s) * meson_energy(s);

            return _masses._mX2 + _masses._mB2 - 2. * real(E1E3) + 2. * real(qdotqp) * cos(theta);
        };

        // Scattering angles in t and u channel frames
        double raw_kinematics::z_t(double s, double theta)
        {
            double t = t_man(s, theta);
            double u = u_man(s, theta);

            double result;
            result  = t * (s - u) + (_masses._mB2 - _masses._mX2) * (_masses._mT2 - _masses._mR2);
            result /=  sqrt(Kallen(t, _masses._mX2, _masses._mB2) * Kallen(t, _masses._mT2, _masses._mR2));

            return result;
        };

        double raw_kinematics::u_man(double s, double theta)
        { 
            return _masses._mX2 + _masses._mB2 + _masses._mT2 + _masses._mR2 - s - t_man(s, theta);
        };

        double raw_kinematics::z_u(double s, double theta)
        {
            double t = t_man(s, theta);
            double u = u_man(s, theta);

            double result;
            result  = u * (t - s) + (_masses._mB2 - _masses._mR2) * (_masses._mT2 - _masses._mX2);
            result /=  sqrt(Kallen(u, _masses._mR2, _masses._mB2) * Kallen(u, _masses._mT2, _masses._mX2));

            return result;
        };

        // ------------------------------------------------------------------------------
        // Momenta in the t-channel center-of-mass frame
        
        complex raw_kinematics::initial_momentum_tframe(double t)
        {
            return csqrt(Kallen(t, _masses._mB2, _masses._mX2)/(4.*t) );
        };

        complex raw_kinematics::final_momentum_tframe(double t)
        {
            return csqrt( Kallen(t, _masses._mR2, _masses._mT2)/(4.*t) );
        };

        // ------------------------------------------------------------------------------
        // Phases associated with helicity relations 

        // Intrisnic parity obeyed by helicity amplitude 
        // This depends on which scattering chan we look at and quantum numbers of all particles
        int raw_kinematics::intrinsic_parity(helicity_frame channel)
        {
            int s_a, s_b, s_c, s_d;
            int eta_a, eta_b, eta_c, eta_d;

            // a is always the photon
            s_a = 2; eta_a = -1; // spin multiplied by two because of spin 1/2 baryons

            switch (channel)
            {
                case helicity_frame::S_CHANNEL :
                {
                    s_b =  1;            eta_b = +1;         // proton
                    s_c =  2*_mjp[0];    eta_c = _mjp[1];   // produced meson
                    s_d =  _bjp[0];      eta_d = _bjp[1];   // recoil baryon
                    break;
                }
                case helicity_frame::T_CHANNEL :
                {
                    s_b =  2*_mjp[0];   eta_b = _mjp[1];    // produced meson
                    s_c =  1;           eta_c = +1;          // proton
                    s_d =  _bjp[0];     eta_d = _bjp[1];    // recoil baryon
                    break;
                }
                case helicity_frame::U_CHANNEL :
                {
                    s_b =  _bjp[0];      eta_b = _bjp[1];    // recoil baryon
                    s_c =  1;            eta_c = +1;          // proton
                    s_d =  2*_mjp[0];    eta_d = _mjp[1];    // produced meson
                    break;
                }

                default: { return 0; }
            };

            int eta = eta_a * eta_b * eta_c * eta_d * pow(-1, (s_c + s_d - s_a - s_b)/2);
            
            return eta;
        };

        // Phase relating lambda_gamma = +1 and lambda_gamma = -1 amplitudes 
        // Depends on the channel with respect to which the helicities are defined
        int raw_kinematics::parity_phase(std::array<int, 4> helicities, helicity_frame channel)
        {
            int lam, lamp;
            switch (channel)
            {
                case helicity_frame::S_CHANNEL :
                {
                    lam =  (2 * helicities[0] - helicities[1]);
                    lamp = (2 * helicities[2] - helicities[3]);
                    break;
                }
                case helicity_frame::T_CHANNEL :
                {
                    lam =  (2 * (helicities[0] - helicities[2]));
                    lamp = (helicities[1] - helicities[3]);
                    break;
                }
                case helicity_frame::U_CHANNEL :
                {
                    lam =  (2 * helicities[0] - helicities[3]);
                    lamp = (2 * helicities[2] - helicities[1]);
                    break;
                }

                default: { return 0; }
            };

            int eta = intrinsic_parity(channel) *  pow(-1, (lam - lamp)/2 );
            return eta;
        };

        // Wigner rotation angle connecting helicity and gottfried-jackson frames 
        double raw_kinematics::H_to_GJ_angle(double s, double t)
        {
            double beta = final_momentum(s) / meson_energy(s);
            double zs = z_s(s, t);

            double cosAlpha = (beta - zs) / (beta * zs - 1.);

            return TMath::ACos( cosAlpha );
        };
    };

    namespace two_meson
    {
        // ---------------------------------------------------------------------------
        // Given a set of helicities, return the index its saved under
        int raw_kinematics::helicity_index(std::array<int,3> helicities)
        {
            auto iterator = std::find(_helicities.begin(), _helicities.end(), helicities);
            if (iterator != _helicities.end()) return iterator - _helicities.begin();
            return error("find_helicity", "Cannot find helicities: " + jpacPhoto::print_helicities(helicities) + "!", -1);
        };    

        // ---------------------------------------------------------------------------
        // Energies of all particles in GJ frame

        // Beam 
        double raw_kinematics::EB(double s, double t, double s12)
        {
            return (s12 + _mB*_mB - t) / (2.*sqrt(s12));
        };

        // Target 
        double raw_kinematics::ET(double s, double t, double s12)
        {
            double u = s12 + _mB*_mB + _mR*_mR + _mT*_mT - s - t;
            return (s12 + _mT*_mT - u) / (2.*sqrt(s12));
        };

        // Recoil 
        double raw_kinematics::ER(double s, double t, double s12)
        {
            return - (s12 + _mR*_mR - s) / (2.*sqrt(s12));
        };

        // These particles are in their CM frame
        double raw_kinematics::E1(double s, double t, double s12)
        {
            return (s12 + _m1*_m1 - _m2*_m2) / (2.*sqrt(s12));
        };
        double raw_kinematics::E2(double s, double t, double s12)
        {
            return (s12 + _m2*_m2 - _m1*_m1) / (2.*sqrt(s12));
        };

        // ---------------------------------------------------------------------------
        // Angles between beam and baryons

        // Angle of beam and target
        double raw_kinematics::cosXi(double s, double t, double s12)
        {
            double E_b = EB(s, t, s12), q_b = sqrt(E_b*E_b - _mB*_mB);
            double E_t = ET(s, t, s12), q_t = sqrt(E_t*E_t - _mT*_mT);

            return (s - 2.*E_b*E_t - _mT*_mT - _mB*_mB) / (2.*q_b*q_t);
        };

        double raw_kinematics::cosEps(double s, double t, double s12)
        {
            double E_b = EB(s, t, s12), q_b = sqrt(E_b*E_b - _mB*_mB);
            double E_r = ER(s, t, s12), q_r = sqrt(E_r*E_r - _mR*_mR);

            return (s - 2.*E_b*E_r - _mR*_mR - s12 + t) / (2.*q_b*q_r);
        };

        // ---------------------------------------------------------------------------
        // Subsystem invariant masses for one meson and the recoil 

        // Particle 1 - Recoil subsystem
        double raw_kinematics::s1(double s, double t, double s12, double thetaGJ, double phiGJ)
        {
            return s + _m1*_m1 + _m2*_m2 + _mR*_mR - s2(s, t, s12, thetaGJ, phiGJ) - s12;
        };
        
        // Particle 2 - Recoil subsystem
        double raw_kinematics::s2(double s, double t, double s12, double thetaGJ, double phiGJ)
        {
            double E_r = ER(s, t, s12), q_r = sqrt(E_r*E_r - _mR*_mR);
            double E_2 = E2(s, t, s12), q_2 = sqrt(E_2*E_2 - _m2*_m2);

            double cos_eps = cosEps(s, t, s12);
            double sin_eps = -sqrt(1. - cos_eps*cos_eps);

            return _m2*_m2 + _mR*_mR + 2.*E_r*E_2 - 2.*q_r*q_2*(sin_eps*sin(thetaGJ)*cos(phiGJ) + cos_eps*cos(thetaGJ));
        };

        // ---------------------------------------------------------------------------
        // Momentum tranfers between beam and a meson

        // Beam - particle 1 transfer
        double raw_kinematics::t1(double s, double t, double s12, double thetaGJ, double phiGJ)
        {
            double E_1 = E1(s, t, s12), q_1 = sqrt(E_1*E_1 - _m1*_m1);
            double E_b = EB(s, t, s12), q_b = sqrt(E_b*E_b - _mB*_mB);
            return _mB*_mB + _m1*_m1 - 2.*E_1*E_b + 2.*q_1*q_b*cos(thetaGJ);
        };

        // Beam - particle 2 transfer
        double raw_kinematics::t2(double s, double t, double s12, double thetaGJ, double phiGJ)
        {
            return t + _m1*_m1 + _m2*_m2 + _mB*_mB - t1(s, t, s12, thetaGJ, phiGJ) - s12;
        };
    };
};