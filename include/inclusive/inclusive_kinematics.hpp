// Extension of the reaction_kinematics class to include quantities revevant
// for semi-inclusive reactions at high energies.
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef INC_KINEM
#define INC_KINEM

#include "constants.hpp"
#include "misc_math.hpp"

namespace jpacPhoto
{
    class inclusive_kinematics 
    {
        public:
        // Empty constructor
        inclusive_kinematics()
        {};

        // Construtor with produced meson mass
        inclusive_kinematics(double mX)
        : _mX(mX), _mX2(mX * mX)
        {};

        inline void set_minM2(double xM2)
        {
            _minM2 = xM2;
        };

        inline double get_mT2()
        {
            return _mT2;
        };

        double _mX, _mX2;                         // Mass of the produce (observed particle)
        double _mT = M_PROTON, _mT2 = M2_PROTON;  // Mass of the target

        double _minM2 = M2_PROTON; // Default to proton is minimum mass unobserved

        // ---------------------------------------------------------------------------
        // Center-of-mass kinematics

        // Momenta
        inline double pX(double s, double M2)
        {
            return sqrt(Kallen(s, _mX2, M2)) / (2. * sqrt(s));
        };

        inline double pX_max(double s)
        {
            return sqrt(Kallen(s, _mX2, _minM2)) / (2. * sqrt(s));
        };

        inline double pGamma(double s)
        {
            return sqrt(Kallen(s, 0., M2_PROTON)) / (2. * sqrt(s));
        };

        // Kinematic bounds for t and M2 with the other fixed
        inline double t_bounds(int pm, double s, double M2)
        {
            double result; 
            result  = _mX2 - (s - _mT2) * (s - M2 + _mX2) / (2. * s);
            result += double(pm) * 2. * pX(s, M2) * pGamma(s);

            return result;
        };

        inline double t_bounds(int pm, double s)
        {
            return (-_minM2*_mT2 + _mT2*_mX2 + (_minM2+_mT2+_mX2)*s-s*s-double(pm)*(s-_mT2)*sqrt(_minM2+pow(_mX2 - s, 2.)-2.*_minM2*(_mX2 + s))) /2./s;
        };

        inline double M2_bounds(int pm, double s, double t)
        {
            if (pm == -1)
            {
                return _minM2;
            } 
            else 
            {
                return (_mT2 + _mX2 - s - t) * (_mT2*_mX2 - s*t) / (_mT2 - s) / (_mX2 - t);
            }
        };

        inline double M2_bounds(int pm, double s)
        {
            if (pm == -1) return _minM2;
            else return pow((sqrt(s) - _mX), 2.);
        };

        // ---------------------------------------------------------------------------
        // Starting from (x, pT2) -> (t, M2)

        inline double M2_from_xpT2(double s, double x, double pT2)
        {
            double pPara = x * pX_max(s);
            double EX    = sqrt(_mX2 + pT2 + pPara*pPara);

            return s + _mX2 - 2.*sqrt(s)*EX;
        };

        inline double cosTheta_from_xpT2(double s, double x, double pT2)
        {
            double M2 = M2_from_xpT2(s, x, pT2);
            return x * pX_max(s) / pX(s, M2);
        };

        inline double t_from_xpT2(double s, double x, double pT2)
        {
            double M2       = M2_from_xpT2(s, x, pT2);
            double cosTheta = cosTheta_from_xpT2(s, x, pT2); 
            double sigma = _mT2 + _mX2 + M2;

            double result;
            result  = 2. * pGamma(s) * pX(s, M2) * cosTheta;
            result += (sigma - s) / 2. + _mT2*(_mX2 - M2)/(2.*s);
            
            return result;
        };
    };
};

#endif