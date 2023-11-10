// Class to incorporate a toy experiment (beam spectrum, acceptance, etc)
// into the jpacPhoto::EventGenerator
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU),
//               University of Bonn (HISKP)
// Email:        winney@hiskp.uni-bonn.de
// ------------------------------------------------------------------------------

#ifndef EXPERIMENT_HPP
#define EXPERIMENT_HPP

#include "IUAmpTools/Kinematics.h"
#include <string>

namespace jpacPhoto
{
    class Experiment
    {
        public: 
        Experiment(std::string id = "")
        : _id(id)
        {};

        // Sample a beam energy from whatever spectrum we want
        virtual double beam_energy(){ return _constantEbeam; };

        // Get the degree of polarization at a given energy
        virtual double polarization(double E){ return 0.; };

        // Add a toy acceptance in accept/reject
        virtual bool acceptance(Kinematics * event){ return true; };

        // String id to identify the setup
        std::string id(){ return _id; };

        double _constantEbeam = 9.0;
        std::string _id;
    };
};

#endif