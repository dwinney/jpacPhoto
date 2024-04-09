// Extention of the amplitude structure but for inclusive reactions,
// here things are much easier since everything needs to be implements
// at the cross-section level not the amplitude level
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------

#ifndef SEMI_INCLUSIVE_HPP
#define SEMI_INCLUSIVE_HPP

#include "constants.hpp"
#include "key.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"

#include "Math/IntegratorMultiDim.h"
#include "Math/GSLIntegrator.h"
#include "Math/IntegrationTypes.h"
#include "Math/Functor.h"

#include <memory>
#include <functional>
#include <vector>

namespace jpacPhoto
{
    // Forward declare for typedef 
    class raw_semi_inclusive;

    // Similar to amplitudes, we only want our inclusive objects to exist as pointers 
    using semi_inclusive = std::shared_ptr<raw_semi_inclusive>;

    //----------------------------------------------------------------------------------------------
    // Use these functions as our constructor
    template<class A>
    inline semi_inclusive new_semi_inclusive(kinematics mX, std::string id)
    {
        auto ptr = std::make_shared<A>(key(), mX, id);
        return std::static_pointer_cast<raw_semi_inclusive>(ptr);
    };

    template<class A, class B>
    inline semi_inclusive new_semi_inclusive(kinematics mX, B parameter, std::string id)
    {
        auto ptr = std::make_shared<A>(key(), mX, parameter, id);
        return std::static_pointer_cast<raw_semi_inclusive>(ptr);
    };
    
    //----------------------------------------------------------------------------------------------
    // Methods for adding terms together

    bool are_compatible(semi_inclusive a, semi_inclusive b);
    bool are_compatible(semi_inclusive a, amplitude b);
    bool are_compatible(amplitude a, semi_inclusive b);

    semi_inclusive operator+(semi_inclusive a, semi_inclusive b);
    semi_inclusive operator+(semi_inclusive a, amplitude b);
    inline semi_inclusive operator+(amplitude a, semi_inclusive b){ return b + a; };
    void operator+=(semi_inclusive a, amplitude b);

    //----------------------------------------------------------------------------------------------

    class raw_semi_inclusive
    {
        public: 

        // Set both observed particle and target masses
        raw_semi_inclusive(key k, kinematics kinem, std::string id)
        : _kinematics(kinem), _mX2(_kinematics->get_meson_mass()*_kinematics->get_meson_mass()), 
          _id(id)
        {};

        raw_semi_inclusive(key k, kinematics kinem, std::vector<semi_inclusive> x, std::vector<amplitude> y, std::string id)
        : _kinematics(kinem), _mX2(_kinematics->get_meson_mass()*_kinematics->get_meson_mass()),  
          _inclusives(x), _exclusives(y),
          _id(id)
        {};
        
        // Access internals
        inline double get_observed_mass(){ return sqrt(_mX2); };
        inline double get_target_mass()  { return sqrt(_mT2); };
        inline std::string id(){ return _id; };

        // Set internals
        inline void set_observed_mass(double m){ _mX2 = m*m; };
        inline void set_target_mass(  double m){ _mT2 = m*m; };
        inline void set_id(std::string id){ _id = id; };

        // Every implementation needs to set the minimum missing mass 
        virtual double minimum_M2(){ return _inclusives[0]->minimum_M2(); };

        // Number of free parameters
        inline int N_pars(){ return _N_pars; };

        // If we are a sum
        inline bool is_sum(){ return (_inclusives.size() > 0) || (_exclusives.size() > 0); };

        // This function is what a user actually calls
        // It wraps the protected vitual method allocate_parameters() with checks of correct size and caching
        void set_parameters( std::vector<double> x );
        inline void set_parameters( double x )
        { 
            std::vector<double> pars = {x};
            set_parameters(pars); 
        };

        // Given a vector of double of appropriate length, allocate free parameters to model
        // By default we do nothing
        virtual void allocate_parameters( std::vector<double> x ){ return; };

        // Pass a flag and make the appropriate changes, defaults to do nothing excpet save the flag
        virtual inline void set_option( int opt ){ _option = opt; };

        // Specify whether our cross section is reggeized
        virtual inline void reggeized(bool x){ _regge = x; };
        
        // ----------------------------------------------------------------------
        // Kinematics 

        kinematics get_kinematics(){ return _kinematics; };

        // Momentum of the photon
        double qGamma(double s);

        // Max momentum of the produced particle (exclusive limit)
        double pMax(double s);

        // POLAR VARIABLES (r, cos)
        double EfromRCOS(   double s, double r, double cos);
        double jacobianRCOS(double s, double r, double cos);

        // CARTESIAN VARIABLES (x, y)
        double EfromXY(    double s, double x, double y);
        double jacobianXY( double s, double x, double y);
        double EfromXY2(   double s, double x, double y2);
        double jacobianXY2(double s, double x, double y2);
        double XfromRCOS(  double s, double r, double cos);
        double XfromTM2(   double s, double t, double M2);

        // INVARIANT variables (t, M2)
        double pXfromM2(   double s, double M2);
        double COSfromTM2( double s, double t,  double M2);
        double RfromTM2(   double s, double t,  double M2);
        double jacobianTM2(double s, double t,  double M2);
        double TfromM2COS( double s, double M2, double cos);
        double M2fromRCOS( double s, double r,  double cos);
        double TfromRCOS(  double s, double r,  double cos);
        double TMINfromM2( double s, double M2);
        double TMAXfromM2( double s, double M2);
        double M2MINfromT( double s, double t);
        double M2MAXfromT( double s, double t);

        // MIXED variables (t, x)
        double M2fromXY(  double s, double x, double y);
        double M2fromXY2( double s, double x, double y2);
        double TfromXY(   double s, double x, double y);
        double TfromXY2(  double s, double x, double y2);
        double jacobianTX(double s, double t, double x);
        double M2fromTX(  double s, double t, double x);
        double TMINfromX( double s, double x);
        double TMAXfromX( double s, double x);

        //--------------------------------------------------------------------
        // d3sigma/d3p (invariant cross-section)
        // These need to be specified by a specific parameterization
        // The third argument (double mm) can be either M2 or x
        // which is assumed to correspond to the _useTX member above

        // Alias function to rename d3sigma_d3p to more human name
        virtual double invariant_xsection(double s, double t, double mm);

        // Doubly differential 
        double dsigma_dtdM2(double s, double t, double M2);
        double dsigma_dtdx( double s, double t, double x);
        double dsigma_dxdy2(double s, double x, double y2);
        inline double dsigma_dxdy(double s, double x, double y){ return dsigma_dxdy2(s, x, y*y) / (2.*y); };

        // Integrated cross-sections

        // (t, M2)
        double dsigma_dt( double s, double t);    // integrated over M2
        double dsigma_dM2(double s, double M2);   // integrated over t

        // (x, y2)
        double dsigma_dy2(double s,  double y2);   // integrated over x
        double dsigma_dx( double s,  double x);    // integrated over pT2
        inline double dsigma_dy(double  s,  double y ){ return dsigma_dy2(s, y*y) / (2.*y); };

        // Fully integrated
        double integrated_xsection(double s, double cos_min = -1.);   
        
        protected:

        // Instance of a kinematics object
        kinematics _kinematics;

        // Same as amplitude save state variables to be able to pass around
        // This also makes sure masses are synced with the kinematics object
        inline void store(double s, double t, double xm2)
        {
            _s = s; _t = t; 
            // Since xm2 can either be x or M2, save the same value twice to avoid confusion
            _x = xm2; _m2 = xm2;

            // Sync masses
            _mX2 = _kinematics->get_meson_mass()*_kinematics->get_meson_mass();
            _mT2 = _kinematics->get_target_mass()*_kinematics->get_target_mass();
        };
        // Total center-of-mass energy
        double _s = 0., _t = 0., _x = 0., _m2 = 0;

        // Int flag used to signal desired changes to amplitude
        int _option = 0;

        // Whether or not the inclusive process is reggeized
        bool _regge = false;
        
        // Mass of the observed and target particles
        double _mX2 = 0, _mT2 = M2_PROTON;

        // Different parameterizations may use different variables which make things tricky
        // when integrating. So I include this flag (set to default as false):
        // False -> assume independent variables are t and M2 
        // True  -> assume independent variables are t and x
        inline virtual bool use_TX(){ return false; };   

        // String ID 
        std::string _id = "semi_inclusive";

        // Number of parameters
        inline void set_N_pars(int Npars){ _N_pars = Npars; };
        int  _N_pars = 0;

        // Simple check that a given vector is of the expected size
        bool correct_size(std::vector<double> pars);

        //---------------------------------------------------------------------------
        // Methods for summing terms together

        friend semi_inclusive operator+(semi_inclusive a, semi_inclusive b);
        friend semi_inclusive operator+(semi_inclusive a, amplitude b);
        friend void operator+=(semi_inclusive a, amplitude b);

        // Sub diagrams / processes
        std::vector<semi_inclusive> _inclusives;
        std::vector<amplitude>      _exclusives;
    };
};

#endif