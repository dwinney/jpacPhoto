// Implementation of a PWA in the scattering-length approximation with up to three
// coupled channels.
//
// This amplitude explicitly uses Eigen C++ to manipulate complex matrices,
// make sure an environment variable EIGEN points to the top level directory
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef SPIN_INDEPENDENT_HPP
#define SPIN_INDEPENDENT_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "utilities.hpp"
#include "partial_wave.hpp"

#include <Eigen/Dense>


namespace jpacPhoto
{
    namespace kmatrix
    {
        struct arguments 
        {
            arguments(uint j, std::array<uint,3>exp={1,1,1}) 
            : _spin(j), 
            _production_expansion(exp[0]),
            _diagonal_elastic_expansion(exp[1]),
            _off_diagonal_elastic_expansion(exp[2])
            {};
            
            void add_coupled_channel(double x, double y){ _coupled_channels.push_back({x,y}); };
            
            // Spin of the partial wave
            uint  _spin;
            // How many terms to consider in the production vector
            uint  _production_expansion  = 1; 
            // How many terms to include in the K-matrix for the diagonal and off-diagonal entries
            uint  _diagonal_elastic_expansion = 1, _off_diagonal_elastic_expansion = 1;
            // If this is coupled channel or not
            std::vector<std::array<double,2>> _coupled_channels;
        };

        class spin_independent : public raw_partial_wave
        {
            public: 

            // Single channel K-matrix
            spin_independent(key k, kinematics xkinem, arguments args)
            : raw_partial_wave(k, xkinem, args._spin, "kmatrix::spin_independent")
            {
                // Populate the thresholds
                _thresholds.push_back({xkinem->get_meson_mass(), xkinem->get_recoil_mass()});
                for (auto extra_threshold : args._coupled_channels)
                {
                    _thresholds.push_back({extra_threshold[0], extra_threshold[1]});
                };

                double nchan = _thresholds.size();
                int offdiags = (nchan-1)*nchan/2;
                _n_prod    = args._production_expansion;
                _n_diag    = args._diagonal_elastic_expansion;
                _n_offdiag = args._off_diagonal_elastic_expansion;
                
                // Total number of parameters
                int npars = (_n_prod+_n_diag)*nchan + _n_offdiag*offdiags;
                initialize(npars);
            };

            // -----------------------------------------------------------------------
            // Virtuals 

            // We can have any quantum numbers
            inline std::vector<quantum_numbers> allowed_mesons() { return {ANY}; };
            inline std::vector<quantum_numbers> allowed_baryons(){ return {ANY}; };
            // And helicity independent
            inline helicity_frame native_helicity_frame(){ return HELICITY_INDEPENDENT; };

            // Partial wave comes from K-matrix unitarized form
            inline complex partial_wave(std::array<int,4> helicities, double s)
            {
                store(helicities, s, _t);
                int N  = _thresholds.size();
                int Np = N*_n_diag; // Number of elastic parameters

                // Set up Q-vector
                Eigen::VectorXcd Q = Eigen::VectorXcd::Zero(N);
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < _n_prod; j++) Q(i) += _production_pars[i*_n_prod+j]*pow(p(0)*q(i) , _J+j);
                };

                // Set up K-matrix
                Eigen::MatrixXcd K   = Eigen::MatrixXcd::Zero(N,N), G = Eigen::MatrixXcd::Zero(N,N);
                Eigen::MatrixXcd One = Eigen::MatrixXcd::Identity(N,N);
                
                // Populate diagonals
                for (int i = 0; i < N; i++)
                {
                    G(i,i) = i_rho(i);
                    for (int j = 0; j < _n_diag; j++) K(i,i) += _elastic_pars[i*_n_diag+j]*pow(q(i)*q(i), _J+j);
                };
                // Populate off-diagonals
                for (int i = 0; i < N; i++)
                {
                    for (int j = i+1; j < N-i; j++)
                    {
                        for (int k = 0; k < _n_offdiag; k++) K(i,j) += _elastic_pars[Np+i*_n_offdiag+k]*pow(q(i)*q(j), _J+k);
                        K(j,i) = K(i,j); // Symmetrize
                    }
                };
                
                auto T = K*(One-G*K).inverse();
                auto F = (One+G*T)*Q;

                return F(0);
            };

            inline void allocate_parameters(std::vector<double> pars)
            {
                _production_pars.clear(); _elastic_pars.clear();
                for (int i = 0; i < pars.size(); i++)
                {
                    if (i < _n_prod*_thresholds.size()) _production_pars.push_back(pars[i]); 
                    else                                _elastic_pars.push_back(pars[i]);
                };
            };

            protected:

            int _n_prod, _n_diag, _n_offdiag;

            // Save the different parameters
            std::vector<double> _production_pars, _elastic_pars;

            // Mass of intermediate coupled channels
            // we can have up to two additional channels
            std::vector<std::array<double,2>> _thresholds;
            
            inline complex i_rho(unsigned i){ return i_rho(_thresholds[i][0], _thresholds[i][1]); };

            // Incoming break-up momentum (define it with a threshold index but we only need i=0)
            inline complex p(int i){ return (i==0) ? _kinematics->initial_momentum(_s) : 0; };
            // Outgoing break-up momentum
            inline complex q(double m1, double m2)
            {
                return csqrt(kallen(_s, m1*m1, m2*m2)) / csqrt(4.*_s);
            };
            inline complex q(unsigned i){ return q(_thresholds[i][0], _thresholds[i][1]); };
        };
    };
};

#endif