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

        inline double M2(double s, double x, double pT2)
        {
            double pPara = x * pX_max(s);
            double EX    = sqrt(_mX2 + pT2 + pPara*pPara);

            return s + _mX2 - _minM2 - 2.*sqrt(s)*EX;
        };

        inline double cosTheta(double s, double x, double pT2)
        {
            return x * pX_max(s) / pX_CM(s, M2(s, x, pT2));
        };

        inline double t_man(double s, double x, double pT2)
        {
            double m2       = M2(s, x, pT2);
            double CosTheta = cosTheta(s, x, pT2); 

            double pGamma = sqrt(Kallen(s, 0., M2_PROTON)) / (2. * sqrt(s));
            double pX     = pX_CM(s, m2);

            double sigma = _mT2 + _mX2 + m2;

            double result;
            result  = 2. * pGamma * pX * CosTheta;
            result += (sigma - s) / 2. + _mT2*(_mX2 - m2)/(2.*s);
            
            return result;
        };


        inline double pX_CM(double xs, double xM2)
        {
            return sqrt(Kallen(xs, _mX2, xM2)) / (2. * sqrt(xs));
        };

        inline double pX_max(double xs)
        {
            return sqrt(Kallen(xs, _mX2, _minM2)) / (2. * sqrt(xs));
        };
    };
};

#endif