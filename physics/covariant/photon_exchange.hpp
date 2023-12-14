// Axial vector photoproduction via virtual photon exchange. Bottom coupling
// given in terms of the Sach's form factors of the nucleon based on [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------
// References:
// [1] - https://arxiv.org/abs/1707.09063
// ------------------------------------------------------------------------------

#ifndef PHOTON_EXCHANGE_HPP
#define PHOTON_EXCHANGE_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"
#include "elementwise.hpp"

namespace jpacPhoto
{
    namespace covariant
    {
        class photon_exchange : public raw_amplitude
        {
            public: 

            photon_exchange(key k, kinematics xkinem, std::string id = "photon_exchange")
            : raw_amplitude(k, xkinem, id), _mEx(0)
            {
                initialize(3);
            };

            photon_exchange(key k, kinematics xkinem, double mEx, std::string id = "photon_exchange")
            : raw_amplitude(k, xkinem, id), _mEx(mEx)
            {
                initialize(3);
            };

            // -----------------------------------------------------------------------
            // REQUIRED VIRTUAL FUNCTIONS

            inline complex helicity_amplitude(std::array<int,4> helicities, double s, double t)
            {
                store(helicities, s, t);
                _covariants->update(helicities, s, t);

                // Contract indices
                complex result = contract( top_coupling(), bottom_coupling() );

                // Divide by photon propagator 
                result /= (t - _mEx*_mEx);

                if (!is_zero(_mEx) && !is_zero(_lam)) 
                {
                    result *= beta(t - _kinematics->t_min(s));
                }
                return result;
            };

            // Covariants are always natively s-channel helicity amplitudes
            inline helicity_frame native_helicity_frame(){  return S_CHANNEL; };

            // Specify which final state particles amplitude can acommodate
            inline std::vector<quantum_numbers> allowed_mesons() { return { AXIALVECTOR }; };
            inline std::vector<quantum_numbers> allowed_baryons(){ return { HALFPLUS    }; };

            inline void allocate_parameters(std::vector<double> pars)
            {
                _gTop = pars[0];
                _eta  = pars[1];
                _lam  = pars[2];
            };

            // Assign each parameter a name, useful for fitting utlities
            inline std::vector<std::string> parameter_labels()
            {
                return {"gTop", "etaEx", "mEx"};
            };

            // Options choose either proton or neutron target
            static const int kProton  = 0;
            static const int kNeutron = 1;
            static const int kDefault = kProton;
            inline void set_option (int opt)
            {
                switch (opt)
                {
                    case kProton:   _option = kProton;  break;
                    case kNeutron:  _option = kNeutron; break;
                    default: return;
                };
            }

            // Sach's electric form factor
            inline double G_E(double Q2)
            {
                double z  = z_conformal(Q2);

                auto pars = (_option == kProton) ? _GEp_pars : _GEn_pars;

                double sum = 0;
                for (int i = 0; i < pars.size(); i++)
                {
                    sum += pars[i] * pow(z, double(i));
                };
                return sum;
            };

            // Sach's magnetic form factor 
            // Normalization is divided by magnetic moment
            inline double G_M(double Q2)
            {
                double z  = z_conformal(Q2);

                auto pars = (_option == kProton) ? _GMp_pars : _GMn_pars;
                double mu = (_option == kProton) ? _mup : _mun;

                double sum = 0;
                for (int i = 0; i < pars.size(); i++)
                {
                    sum += pars[i] * pow(z, double(i));
                };
                return mu*sum;
            }

            // Simple dipole form factor
            inline double G_D(double Q2)
            {
                return 1. / pow(1+Q2/0.71, 2);
            };

            // Ratio of form factors for a massive vector exchange
            inline double beta(double t)
            {
                return exp(t/_lam/_lam) / G_D(-t);
            }

            // -----------------------------------------------------------------------
            // Internal data members 

            protected:

            // A -- gamma -- gamma* coupling
            double _gTop = 0;
            double _eta  = 1;
            double _mEx  = 0;
            double _lam  = 1;

             // Top coupling refers to the beam-gamma-meson interaction
            inline lorentz_tensor<complex,1> top_coupling()
            {
                // Beam
                auto k     = _covariants->q();
                auto eps   = _covariants->eps();
                
                // Outgoing meson
                auto eps_p = _covariants->eps_prime();

                // Exchange momentum
                auto q = _covariants->k_t();

                // Coupling function depends on
                // the quantum numbers of the produced meson
                lorentz_tensor<complex,1> T = _t*levi_civita(k, eps, eps_p) - levi_civita(k, eps, q, eps_p) * q; 
                return _eta*_gTop/_mX/_mX * T;
            };

             // Bottom coupling refers to the gamma^* NNbar interation
            inline lorentz_tensor<complex,1> bottom_coupling()
            {
                // Recalculate form_factors 
                double GE  =  G_E(-_t), GM = G_M(-_t);
                double tau =  -_t/(4*M2_PROTON);

                auto u     = _covariants->u();    // Target spinor
                auto ubar  = _covariants->ubar(); // Recoil spinor
                auto q     = _covariants->k_t();  // t-channel exchange momentum

                // This is also dependent on baryon quantum numbers
                // We have two pieces for the vector and tensor currents
                lorentz_tensor<complex,1> vector, tensor;

                auto vector_current =  gamma_vector();
                auto tensor_current = (gamma_vector()*slash(q) - slash(q)* gamma_vector())/2;

                vector = bilinear(ubar, vector_current, u); 
                tensor = bilinear(ubar, tensor_current, u);

                return _eta * E*((GE+tau*GM)*vector + (GM-GE)/(2*M_PROTON)*tensor)/(1+tau);
            };

            //---------------------------------------------------------------------
            // Rest is related to the form factors G_E and G_M

            // Conformal expansion variable
            inline double z_conformal(double Q2)
            {
                double tcut = 4*M2_PION,       t0 = -0.7;
                double x    = sqrt(tcut + Q2),  y = sqrt(tcut - t0);

                return (x - y) / (x + y);
            };

            static constexpr double _mup =  2.79*0.161; // Proton  magnetic moment in GeV-1
            static constexpr double _mun = -1.92*0.161; // Nucleon magnetic moment in GeV-1

            // Polynomial coefficients
            std::array<double,13> _GEp_pars = 
            {
                 0.239163298067,  -1.109858574410,  1.444380813060,  0.479569465603,
                -2.286894741870,   1.126632984980,  1.250619843540, -3.631020471590,
                 4.082217023790,   0.504097346499, -5.085120460510,  3.967742543950, -0.981529071103
            };
            std::array<double,13> _GEn_pars = 
            {
                 0.048919981379, -0.064525053912, -0.240825897382, 0.392108744873,
                 0.300445258602, -0.661888687179, -0.175639769687, 0.624691724461,
                -0.077684299367, -0.236003975259,  0.090401973470
            };

            std::array<double,13> _GMp_pars = 
            {
                  0.264142994136, -1.095306122120,  1.218553781780,  0.661136493537,
                 -1.405678925030, -1.356418438880,  1.447029155340,  4.235669735900,
                 -5.334045653410, -2.916300520960,  8.707403067570, -5.706999943750,  1.280814375890
            };
            std::array<double,13> _GMn_pars =
            {
                 0.257758326959, -1.079540642058,  1.182183812195, 0.711015085833, 
                -1.348080936796, -1.662444025208,  2.624354426029, 1.751234494568,
                -4.922300878888,  3.197892727312, -0.712072389946 
            };
        };
    };
};

#endif