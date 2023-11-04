// Class to contain all relevant kinematic quantities. 
// (e.g. angles, momenta, helicity info)
//
// Unlike single_meson::kinematics this class represents a 2->3 process
// with two pseudo-scalars and a proton in the final state
// Generalizations to more complicated final states can be implemented lated
//
// Since we dont need to specify spins, everything is determined by each external mass:
// _mB = beam particle mass (default 0)
// _mT = target recoil mass (default M_PROTON)
// _m1 = "meson 1"
// _m2 = "meson 2"
// _mR = recoil baryon mass (default M_PROTON)
//
// --------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
//               University of Bonn (HISKP)
// Email:        winney@hiskp.uni-bonn.de
// --------------------------------------------------------------------------------

#include "kinematics2.hpp"

namespace jpacPhoto
{
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
            return (s12 + _mR*_mR - s) / (2.*sqrt(s12));
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
            double s_2 = s2(s, t, s12, thetaGJ, phiGJ);
            return s - s_2 - s12 + _m1*_m1 + _m2*_m2 + _mR*_mR;
        };
        
        // Particle 2 - Recoil subsystem
        double raw_kinematics::s2(double s, double t, double s12, double thetaGJ, double phiGJ)
        {
            double E_r = ER(s, t, s12), q_r = sqrt(E_r*E_r - _mR*_mR);
            double E_2 = E2(s, t, s12), q_2 = sqrt(E_2*E_2 - _m2*_m2);

            double cos_eps = cosEps(s, t, s12);
            double sin_eps = sqrt(1. - cos_eps*cos_eps);

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
            double t_1 = t1(s, t, s12, thetaGJ, phiGJ);
            return t - t_1 - s12 + _m1*_m1 + _m2*_m2 + _mB*_mB;
        };
    };
};