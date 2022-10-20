// Amplitude assuming a simple effective-range expansion near-threshold
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef EFF_RANGE
#define EFF_RANGE

#include "amplitude.hpp"

namespace jpacPhoto
{
    class effective_range : public amplitude
    {
        // ---------------------------------------------------------------------------
        
        public:

        effective_range(reaction_kinematics * xkinem, int ellmax, std::string id = "effective_range")
        : amplitude(xkinem, "effective_range", id), _ellMax(ellmax)
        {
            set_nParams(3);
            check_JP(xkinem);
        };

        // Set the free parameters
        inline void set_params(std::vector<double> params)
        {
            check_nParams(params);

            _N = params[0];
            _a = params[1];
            _r = params[2];

            _extra_couplings.clear();
            for (int n_th = 0; n_th < _extra_thresholds.size(); n_th++)
            {
                _extra_couplings.push_back(params[3+n_th]);
            };
        };

        inline void add_threshold(double m1, double m2)
        {
            // Save the masses of the new threshold
            _extra_thresholds.push_back({m1, m2});

            // and add one to the number of free parameters
            set_nParams( get_nParams() + 1 );
        };  

        // Assemble the helicity amplitude by summing with appropriate angular funxtions
        std::complex<double> helicity_amplitude(std::array<int, 4> helicities, double xs, double xt);

        // Expansion is directly in the S channel so helicities should always be defined in this frame
        helicity_channel helicity_CM_frame(){ return helicity_channel::S; };

        // only axial-vector, vector, and pseudo-scalar available
        inline std::vector<std::array<int,2>> allowed_meson_JP(){  return { {1, -1} }; };
        inline std::vector<std::array<int,2>> allowed_baryon_JP(){ return { {1, +1} }; };

        // ---------------------------------------------------------------------------

        protected:

        // Number of partial-wave terms to consider 
        int _ellMax = 1;

        // Free parameters
        double _a;         // Scattering length (assumed same for all waves)
        double _r;         // Effective range (assumed same for all waves)
        double _N;         // Overall normalization

        // A vector containing all open thresholds
        std::vector<std::array<double,2>> _extra_thresholds; 
        std::vector<double> _extra_couplings;

        // Legendre functions
        double P_l(int l, double z);

        // Partial wave amplitude
        std::complex<double> f_l(int l);
        
        // Redefine momenta here instead of using the reaction_kinematics versions
        // so that they can be appropriately analytically continued below threshold
        inline std::complex<double> barrier_factor(int l)
        {         
            std::complex<double> pq = sqrt( Kallen(_s * XR, _mB*_mB * XR, _mT*_mT * XR)*Kallen(_s * XR, _mX*_mX * XR, _mR*_mR * XR) ) / (4.*_s);           
            return pow( _r* pq * XR, double(l) );
        };

        // Phase-space factor, this is defined with respect to a threshold to allow multiple thresholds to be considered
        inline std::complex<double> rho()
        {
            return PI * sqrt( Kallen(_s * XR, _mX*_mX * XR, _mR*_mR * XR) ) / (4.*_s);
        };

        inline std::complex<double> rho_inel(int i)
        {
            if (i > _extra_thresholds.size()) return 0.;
            double m1 = _extra_thresholds[i][0], m2 = _extra_thresholds[i][1];
            return PI * sqrt( Kallen(_s * XR, m1*m1 * XR, m2*m2 * XR) ) / (4.*_s);
        };
        
    };
};

#endif