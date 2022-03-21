// Analogous to the reaction_kinematics class for exclusive reactions,
// here all the information of the kinematics of the semi-inclusive reaction are stored.
// This includes conversions between the different sets of variables e.g. (s, t, M2) and (s, x, pT2), etc.
// as well as jacobians and masses.
//
// Author:       Daniel Winney (2022)
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
        inclusive_kinematics(double mX, std::string id = "")
        : _mX(mX), _mX2(mX * mX), _id(id)
        {};

        // String label for the produced particle
        std::string _id;

        double _mX, _mX2;                         // Mass of the produce (observed particle)
        double _mT = M_PROTON, _mT2 = M2_PROTON;  // Mass of the target

        double _minM2 = M2_PROTON; // Default to proton is minimum mass unobserved

        double _s = 100.; // Total center-of-mass energies 

        // ---------------------------------------------------------------------------
        // Photon Quantites
        // These dont depend on anything other than s and are characteristic of the initial state

        // Momentum of the photon
        inline double qGamma()
        {
            return (_s - _mT2) / sqrt(4. * _s);
        };

        // ---------------------------------------------------------------------------
        // Quantities in the exclusive limit
        
        // Max momentum of the produced particle
        inline double pMax()
        {
            return sqrt(Kallen(_s, _mX2, _minM2)) / (2. * sqrt(_s));
        };

        // ---------------------------------------------------------------------------
        // POLAR VARIABLES (r, cos)
        
        // Energy of the particle X
        inline double EfromRCOS(double r, double cos)
        {
            double e2 = _mX2 + pMax()*pMax()*r*r;
            return sqrt(e2);
        };

        // Jacobian in polar coordinates
        inline double jacobianRCOS(double r, double cos)
        {
            double inv = EfromRCOS(r, cos) / (2.* M_PI * (r*r) * pow(pMax(), 3));
            return 1./inv;
        };

        // ---------------------------------------------------------------------------
        // CARTESIAN VARIABLES (x, y)

        // Energy of X in (x, y)
        inline double EfromXY(double x, double y)
        {
            double e2 = _mX2 + pMax()*pMax() * (x*x + y*y);
            return sqrt(e2);
        };

        // Jacobian in (x, y)
        inline double jacobianXY(double x, double y)
        {
            double inv = EfromXY(x, y) / (2.* M_PI * y * pow(pMax(), 3));
            return 1./inv;
        };

        // Energy of X in (x, y2)
        inline double EfromXY2(double x, double y2)
        {
            double e2 = _mX2 + pMax()*pMax() * (x*x + y2);
            return sqrt(e2);
        };

        // Jacobian in (x, y2)
        inline double jacobianXY2(double x, double y2)
        {
            double inv = EfromXY2(x, y2) / (M_PI * pow(pMax(), 3));
            return 1./inv;
        };

        // ---------------------------------------------------------------------------
        // INVARIANT variables (t, M2)

        // momentum of produced at a fixed missing mass
        inline double pXfromM2(double M2)
        {
            return sqrt(Kallen(_s, _mX2, M2)) / (2. * sqrt(_s));
        };  

        // Jacobian in (t, M2)
        inline double jacobianTM2(double t, double M2)
        {
            double inv = 2. * sqrt(_s) * qGamma() / M_PI;
            return 1./inv;
        };

        // t from cos and M2 
        inline double TfromM2COS(double M2, double cos)
        {
            double t = _mX2 - (_s - _mT2) * (_s - M2 + _mX2) / (2. * _s) + 2.* qGamma() * pXfromM2(M2) * cos;
            return t;
        };

        // Bounds of integration in t at fixed M2
        inline double TMINfromM2(double M2) { return TfromM2COS(M2,  1.); };
        inline double TMAXfromM2(double M2) { return TfromM2COS(M2, -1.); };

        // Bounds of integration in M2 at fixed t
        inline double M2MINfromT(double t) { return _minM2; };
        inline double M2MAXfromT(double t) 
        {
            return (_mT2 + _mX2 - _s - t) * (_mT2 * _mX2 - _s * t) / (_mT2 - _s) / (_mX2 - t);
        };

        // ---------------------------------------------------------------------------
        // MIXED variables (t, x)

        // Missing mass from the cartesian XY
        inline double M2fromXY(double x, double y)
        {
            return _s + _mX2 - 2. * sqrt(_s) * EfromXY(x, y);
        }; 
        inline double M2fromXY2(double x, double y2)
        {
            double y = sqrt(y2);
            return _s + _mX2 - 2. * sqrt(_s) * EfromXY(x, y);
        }; 


        // Similarly, momentum transfer from XY
        inline double TfromXY(double x, double y)
        {
            return _mX2 - 2.*qGamma() * (EfromXY(x, y) - pMax() * x);
        };
        inline double TfromXY2(double x, double y2){ return TfromXY(x, sqrt(y2)); };

        // Bounds of integration for T at fixed X
        inline double TMINfromX(double x) { return TfromXY(x, sqrt(1. - x*x)); };
        inline double TMAXfromX(double x) { return TfromXY(x, 0.); };

        // Jacobian for mixed variables
        inline double jacobianTX(double t, double x)
        {
            double inv = qGamma() / M_PI / pMax();
            return 1./inv;
        };

        // FIXME: Bounds of integration for X at fixed T
        inline double XMINfromT(double t)
        {
            return 0.;
        };
        inline double XMAXfromT(double t)
        {
            return 0.;
        };
    };  
};

#endif