// Class to populate a data fileswith simulated events 
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU),
//               University of Bonn (HISKP)
// Email:        winney@hiskp.uni-bonn.de
// ------------------------------------------------------------------------------

#ifndef EVENTGENERATOR_HPP
#define EVENTGENERATOR_HPP

#include <iostream>
#include <string>
#include <vector>
#include <time.h>
#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/report.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/AmpToolsInterface.h"

#include "debug.hpp"
#include "print.hpp"

namespace jpacPhoto
{
    // Need to specify the number of particles in the final state N
    // and the writer type to format the generated data (Writer::DataWriter<N> should match N above)
    template<class Writer, int N>
    class EventGenerator
    {
        public: 

        EventGenerator( std::array<double, N> masses, std::array<std::string,N> particle_labels )
        : _labels(particle_labels)
        {
            // Convert the input masses into a C-string
            for (int i = 0; i < N; i++)
            {
                _masses[i] = masses[i];
            };
        };

        // Set a fixed beam energy to generate phase space / physics with
        // This can later be replaced by a function of beam specturm
        inline void setBeamEnergy(double Egam) 
        {
            double M_PROTON = 0.938272;
            TLorentzVector W = TLorentzVector(0., 0., Egam, Egam + M_PROTON);
            _generator.SetDecay(W, N, _masses);
            ENERGY_SET = true;
        };

        // Generate nEvents phase space points given a fixed beam energy
        // Follows the generatePhaseSpace AmpTools tutorial
        // https://github.com/mashephe/AmpTools/blob/master/Tutorials/Dalitz/DalitzExe/generatePhaseSpace.cc
        inline void generatePhaseSpace(int nEvents, std::string outfile = "out.root")
        {
            // Set random seeds
            gRandom = new TRandom3(0);

            Writer writer(_labels, outfile);
         
            if ( !ENERGY_SET )
            {
                jpacPhoto::warning("jpacTools::EventGenerator", "Beam energy not set!");
                return;
            };

            // Generate events
            double maxWeight = _generator.GetWtMax();
            for (int i = 0; i < nEvents; i++)
            {
                double weight = _generator.Generate();
                if (weight / maxWeight < gRandom->Rndm())
                {
                    i--; continue;
                };
                
                std::vector<TLorentzVector> fvecs;
                for (int n = 0; n < N; n++)
                {
                    fvecs.push_back( TLorentzVector( *_generator.GetDecay(n) ) );
                };
                Kinematics kin(fvecs);
                writer.writeEvent(kin);
            };
        };

        // Generate nEvents phasespace points but additionally save weights corresponding
        // to the intensity specified in configfile
        template<class A>
        inline void generateWeightedPhaseSpace(int nEvents, std::string configfile, std::string outfile = "out.root")
        {
            // Set random seeds
            gRandom = new TRandom3(0);

            ConfigFileParser parser(configfile);
            ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
            cfgInfo->display();
            ReactionInfo* reaction = cfgInfo->reactionList()[0];

            // AmpToolsInterface
            AmpToolsInterface::registerAmplitude( A() );
            AmpToolsInterface ATI(cfgInfo, AmpToolsInterface::kMCGeneration);

            // DataWriter
            Writer writer(_labels, outfile);

            // Generate a phase-space data set
            double maxWeight = _generator.GetWtMax();
            for (int i = 0; i < nEvents; i++)
            {
                double weight = _generator.Generate();
                if (weight / maxWeight < gRandom->Rndm())
                {
                    i--; continue;
                };
                
                std::vector<TLorentzVector> fvecs;
                for (int n = 0; n < N; n++)
                {
                    fvecs.push_back( TLorentzVector( *_generator.GetDecay(n) ) );
                };
                Kinematics* kin = new Kinematics(fvecs);

                // Here instead of load load to the AmpToolsInterface to calculate intensity
                ATI.loadEvent(kin, i, nEvents);
                delete kin;
            };

            // Instead of accept reject simply save with weight
            double maxIntensity = ATI.processEvents(reaction->reactionName());
            for (int i = 0; i < nEvents; i++)
            {
                Kinematics* kin = ATI.kinematics(i);
                kin->setWeight(ATI.intensity(i) / maxIntensity);
                writer.writeEvent(*kin);
                delete kin;
            };
        };

        // Generate nEvents using hit-or-miss with a specified model intensity
        template<class A>
        inline void generatePhysics(int nEvents, std::string configfile, std::string outfile = "out.root")
        {
            // Set random seeds
            gRandom = new TRandom3(0);

            ConfigFileParser parser(configfile);
            ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
            cfgInfo->display();
            ReactionInfo* reaction = cfgInfo->reactionList()[0];

            // AmpToolsInterface
            AmpToolsInterface::registerAmplitude( A() );
            AmpToolsInterface ATI(cfgInfo, AmpToolsInterface::kMCGeneration);

            // DataWriter
            Writer writer(_labels, outfile);

            // Generate events
            line();
            divider();
            print("Generating " + std::to_string(nEvents) + " with amplitude: " + A().name());
            print("Configuration file: " + configfile);
            print("Output file: " + outfile);
            divider();
            int nPS = 5E5; // Generate phase space in chunks

            int nGenerated = 0;
            int Passes = 1, nFailed = 0;
            while (nGenerated < nEvents)
            {
                if (nFailed > 3)
                {
                    warning("EventGenerator::generatePhysics", "Three passes with no events generated, possible infinite loop? Exiting...");
                    exit(1);
                };
                int nPass = 0;
                double maxWeight = _generator.GetWtMax();
                for (int i = 0; i < nPS; i++)
                {
                    double weight = _generator.Generate();
                    if (weight / maxWeight < gRandom->Rndm())
                    {
                        i--; continue;
                    };
                    
                    std::vector<TLorentzVector> fvecs;
                    for (int n = 0; n < N; n++)
                    {
                        fvecs.push_back( TLorentzVector( *_generator.GetDecay(n) ) );
                    };
                    Kinematics* kin = new Kinematics(fvecs);

                    // Here instead of loading to writer, load to the AmpToolsInterface
                    ATI.loadEvent(kin, i, nPS);
                    delete kin;
                };

                // DO accept / reject
                double maxIntensity = ATI.processEvents(reaction->reactionName());
                for (int i = 0; i < nPS; i++)
                {
                    double Intensity = ATI.intensity(i); 
                    if ((nGenerated + nPass) < nEvents && Intensity / maxIntensity > gRandom->Rndm())
                    {
                        nPass++;
                        Kinematics* kin = ATI.kinematics(i);
                        writer.writeEvent(*kin);
                        delete kin;
                    }
                };
                print("-- Pass " + std::to_string(Passes) + ": Generated " + std::to_string(nPass) + " events (" + std::to_string(nGenerated + nPass) + "/" + std::to_string(nEvents) + ")...");
                nGenerated += nPass;
                Passes++;
                if (nPass == 0) nFailed++;
            };
            divider();
            line();
        };

        protected:

        // Final state masses
        double _masses[N];

        // Final state particle labels
        std::array<std::string,N> _labels;

        // Related to TGenPhaseSpace
        bool ENERGY_SET = false;
        TGenPhaseSpace _generator;
    };

};

#endif