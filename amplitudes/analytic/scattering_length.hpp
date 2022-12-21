// Implementation of a PWA in the scattering-length approximation
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef SCATTERING_LENGTH_HPP
#define SCATTERING_LENGTH_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"
#include "partial_wave.hpp"

namespace jpacPhoto
{
    namespace analytic
    {
        class scattering_length : public raw_partial_wave
        {
            public: 

            scattering_length(amplitude_key key, kinematics xkinem, int J, std::string id = "scattering_length")
            : raw_partial_wave(key, xkinem, J, "scattering_length", id)
            {
                set_N_pars(2);
                check_QNs(xkinem);
            };

            // -----------------------------------------------------------------------
            // Virtuals 

            // These are projections onto the orbital angular momentum and therefore
            // helicity independent
            inline complex helicity_amplitude(std::array<int,4> helicities, double s, double t)
            {
                // Arbitrarily pick one of the helicities to evaluate
                if (helicities != _kinematics->helicities(0)) return 0;

                // Save inputes
                store(helicities, s, t);

                // Normalization here to get rid of helicity dependence in amplitude::probability_distribution
                // First a 2 removes the factor 1/4 when averaging over initial helicities
                // then a 1/sqrt(2) removes the factor of 2 from the parity relation in amplitude::update_cache
                return sqrt(2) * legendre(_J, cos(_theta)) * evaluate(helicities, s);
            };

            // We can have any quantum numbers
            // but for now explicitly put only pseudo-scalar, vector, and axial-vector
            inline std::vector<std::array<int,2>> allowed_meson_JP()
            {
                return { PSUEDOSCALAR, VECTOR, AXIALVECTOR };
            };

            // and either parity spin-1/2
            inline std::vector<std::array<int,2>> allowed_baryon_JP()
            {
                return { HALFPLUS, HALFMINUS };
            };

            // Parameter names are a[J] and b[J] for scattering length and normalization respectively
            inline std::vector<std::string> parameter_labels()
            {
                std::string J = std::to_string(_J);
                return{ "a[" + J + "]", "b[" + J + "]" };
            };

            // PW is the unitarized K-matrix form
            inline complex evaluate(std::array<int,4> helicities, double s)
            {
                complex K = pow(q2(), _J) * _a;
                complex B = pow(pq(), _J) * _b;
                complex T = K / (1 - rhoCM()*K);

                return B * (1 + rhoCM()*T);
            };

            protected:

            inline void allocate_parameters(std::vector<double> pars)
            {
                _a = pars[0];
                _b = pars[1];
                return;
            };

            // -----------------------------------------------------------------------
            private:

            // couplings are the overall normalization (b) and scattering length (a)
            double _a = 0, _b = 0; 

            // Chew-Mandelstam phase-space
            inline complex rhoCM(double m1, double m2)
            {
                complex rho, xi;
                complex result;

                rho    = csqrt(Kallen(_s, m1*m1, m2*m2)) / _s;
                xi     = 1 - (m1+m2)*(m1+m2)/_s;
                result = ( rho*log((xi + rho) / (xi - rho)) - xi*(m2-m1)/(m2+m1)*log(m2/m1) ) / PI;
                return result;
            };
            inline complex rhoCM(){ return rhoCM(_mX, _mR); };

            // Product of momenta for the photoproduction process
            inline complex pq()
            {
                return _kinematics->initial_momentum(_s) * _kinematics->final_momentum(_s);
            };
            
            // Product of momenta of the hadronic rescattering process
            inline complex q2()
            {
                return pow(_kinematics->final_momentum(_s), 2);
            };
        };
    };
};

#endif