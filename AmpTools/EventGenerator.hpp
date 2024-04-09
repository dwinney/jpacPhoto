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
#include "constants.hpp"
#include "plotter.hpp"
#include "Experiment.hpp"

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

        EventGenerator( std::array<double, N> masses, std::array<std::string,N> particle_labels, Experiment * x)
        : _labels(particle_labels), _exp(x)
        {
            // Convert the input masses into a C-string
            for (int i = 0; i < N; i++)
            {
                _masses[i] = masses[i];
            };
        };

        inline void setExperiment(Experiment * exp)
        {
            _exp = exp;
        };

        inline void setSeed(int x)
        {
            _seed = x;
        };

        // ------------------------------------------------------------------------------
        // These methods are for generating WEIGHTED events
        // These generate and save the raw 4-vectors from TGenPhaseSpace
        // and save the weights (either of pure phase-space or an amplitude )

        // Generate nEvents phase space points given a fixed beam energy
        // Follows the generatePhaseSpace AmpTools tutorial
        // https://github.com/mashephe/AmpTools/blob/master/Tutorials/Dalitz/DalitzExe/generatePhaseSpace.cc
        inline void generateWeightedPhaseSpace(int nEvents, std::string outfile = "out.root")
        {
            if (_exp == nullptr) 
            {
                warning("generateWeightedPhaseSpace", "No experiment setup set!");
                return;
            }

            // Set random seeds
            gRandom = new TRandom3(_seed);

            Writer writer(_labels, outfile);
         
            // Instead of accept reject just save weight
            for (int i = 0; i < nEvents; i++)
            {
                double egam = _exp->beam_energy();
                TLorentzVector W = TLorentzVector(0., 0., egam, egam + M_PROTON);
                _generator.SetDecay(W, N, _masses);
                double weight = _generator.Generate();

                // Unpack 4-vectors 
                std::vector<TLorentzVector> fvecs;
                for (int n = 0; n < N; n++)
                {
                    fvecs.push_back( TLorentzVector( *_generator.GetDecay(n) ) );
                };

                Kinematics kin = Kinematics(fvecs);
                kin.setWeight(weight); // Set weight and save all events
                writer.writeEvent(kin);
            };
        };

        // Generate nEvents phasespace points but additionally save weights corresponding
        // to the intensity specified in configfile
        template<class A>
        inline void generateWeightedPhysics(int nEvents, std::string configfile, std::string outfile = "out.root")
        {
            if (_exp == nullptr) 
            {
                warning("generateWeightedPhaseSpace", "No experiment setup set!");
                return;
            }

            // Set random seeds
            gRandom = new TRandom3(_seed);

            // Generate events
            line();
            divider();
            print("Generating " + std::to_string(nEvents) + " with amplitude: " + A().name());
            print("Method:", "Weighted");
            print("Beam setup:", _exp->id());
            print("cfg file:", configfile);
            print("Output file:", outfile);
            divider();

            ConfigFileParser parser(configfile);
            ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
            // cfgInfo->display();
            ReactionInfo* reaction = cfgInfo->reactionList()[0];

            // AmpToolsInterface
            AmpToolsInterface::registerAmplitude( A() );
            AmpToolsInterface ATI(cfgInfo, AmpToolsInterface::kMCGeneration);

            // DataWriter
            Writer writer(_labels, outfile);

            // Generate a phase-space data set
            for (int i = 0; i < nEvents; i++)
            {
                double egam = _exp->beam_energy();
                TLorentzVector W = TLorentzVector(0., 0., egam, egam + M_PROTON);
                _generator.SetDecay(W, N, _masses);
                double weight = _generator.Generate();
                std::vector<TLorentzVector> fvecs;
                for (int n = 0; n < N; n++)
                {
                    fvecs.push_back( TLorentzVector( *_generator.GetDecay(n) ) );
                };
                Kinematics* kin = new Kinematics(fvecs);
                kin->setWeight(weight);

                // Here instead of load load to the AmpToolsInterface to calculate intensity
                ATI.loadEvent(kin, i, nEvents);
                delete kin;
            };

            // Instead of accept reject simply save with weight
            double maxIntensity = ATI.processEvents(reaction->reactionName());
            for (int i = 0; i < nEvents; i++)
            {
                Kinematics* kin = ATI.kinematics(i);
                kin->setWeight(ATI.intensity(i) * kin->weight());
                writer.writeEvent(*kin);
                delete kin;
            };

            line();
        };

        // ------------------------------------------------------------------------------
        // These methods are for generating UNWEIGHTED events using hit-or-miss
        // and therefore have weights of 1

        inline std::vector<Kinematics> generateBatch()
        {
            std::vector<Kinematics> out;
            for (int i = 0; i < _nBatch; i++)
            {
                double egam = _exp->beam_energy();
                TLorentzVector W = TLorentzVector(0., 0., egam, egam + M_PROTON);
                _generator.SetDecay(W, N, _masses);
                double weight = _generator.Generate(); 
                std::vector<TLorentzVector> fvecs;
                for (int n = 0; n < N; n++)
                {
                    fvecs.push_back( TLorentzVector( *_generator.GetDecay(n) ) );
                };
                Kinematics kin = Kinematics(fvecs);
                kin.setWeight(weight);
                out.push_back(kin);
            };
            return out;
        };

        inline double getMaxWeight( std::vector<Kinematics> in)
        {
            double max = in[0].weight();
            for (auto event : in)
            {
                if (event.weight() > max) max = event.weight();
            };
            return max;
        };

        inline void generatePhaseSpace(int nEvents, std::string outfile = "out.root")
        {
            if (_exp == nullptr) 
            {
                warning("generateWeightedPhaseSpace", "No experiment setup set!");
                return;
            }

            // Set random seeds
            gRandom = new TRandom3(_seed);

            Writer writer(_labels, outfile);
        
            int nGenerated = 0; // Total points generated so far
            int nFailed = 0;    // Number of passes made with no generated events
            while (nGenerated < nEvents)
            {
                if (nFailed > 3)
                {
                    warning("EventGenerator::generatePhasespace", "Three passes with no events generated, possible infinite loop? Exiting...");
                    exit(1);
                };

                auto batch       = generateBatch();
                double maxWeight = getMaxWeight(batch);
                
                // Do accept / reject
                int nPass = 0; // Number of events generated this pass
                for (auto event : batch)
                {
                    // Terminate early if generated enough events
                    if ((nGenerated + nPass) >= nEvents) break;

                    double weight = event.weight(); 
                    if (weight / maxWeight > gRandom->Rndm())
                    {
                        event.setWeight(1.); // Reset the weight to 1
                        if (_exp->acceptance(&event)) writer.writeEvent(event);
                        nPass++;
                    }
                };
                nGenerated += nPass;
                if (nPass == 0) nFailed++;
            };
        };

        // Generate nEvents using hit-or-miss with a specified model intensity
        template<class A>
        inline void generatePhysics(int nEvents, std::string configfile, std::string outfile = "out.root", double maxwt = -1)
        {
            if (_exp == nullptr) 
            {
                warning("generateWeightedPhaseSpace", "No experiment setup set!");
                return;
            }

            // Set random seeds
            gRandom = new TRandom3(_seed);

            // Generate events
            line();
            divider();
            print("Generating " + std::to_string(nEvents) + " with amplitude: " + A().name());
            print("Method:", "Hit-or-Miss");
            print("Beam setup:", _exp->id());
            print("cfg file:", configfile);
            print("Output file:", outfile);
            divider();

            ConfigFileParser parser(configfile);
            ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
            // cfgInfo->display();
            ReactionInfo* reaction = cfgInfo->reactionList()[0];

            // AmpToolsInterface
            AmpToolsInterface::registerAmplitude( A() );
            AmpToolsInterface ATI(cfgInfo, AmpToolsInterface::kMCGeneration);
            line();

            line();
            divider();
            
            // DataWriter
            Writer writer(_labels, outfile);

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
                // First generate a batch of phase space
                auto batch       = generateBatch();
                double maxWeight = getMaxWeight(batch);
                int nPS = 0;  // Number of generated PS points

                // Do hit-or-miss to generate phasespace
                for (auto event : batch)
                {
                    double weight = event.weight(); 
                    if (weight / maxWeight > gRandom->Rndm())
                    {
                        event.setWeight(1.); // Reset the weight to 1
                        ATI.loadEvent(&event, nPS, batch.size()); // Instead of writing, load to AmpTools
                        nPS++;
                    };
                };

                // Loop over the loaded events and do a second hit-or-miss
                // DO accept / reject
                double maxIntensity = ATI.processEvents(reaction->reactionName());
                if (maxwt > 0) maxIntensity = maxwt;
                for (int i = 0; i < nPS; i++)
                {
                    if ((nGenerated + nPass) >= nEvents) break;

                    double Intensity = ATI.intensity(i); 
                    if (Intensity / maxIntensity > gRandom->Rndm())
                    {
                        if (Intensity > maxIntensity) continue;
                        Kinematics* kin = ATI.kinematics(i);
                        if (_exp->acceptance(kin)) writer.writeEvent(*kin);
                        delete kin;
                        nPass++;
                    }
                };
                print("-- Pass " + std::to_string(Passes) + ": Generated " + std::to_string(nPass) + " events (" + std::to_string(nGenerated + nPass) + "/" + std::to_string(nEvents) + ")...");
                nGenerated += nPass;
                Passes++;
                if (nPass == 0) nFailed++;
                ATI.clearEvents();
            };
            divider();
            line();
        };

        // template<class A>
        // inline void intensityProfile(std::string configfile, std::string outfile, int nEvents = 1E5)
        // {
        //     if (_exp == nullptr) 
        //     {
        //         warning("generateWeightedPhaseSpace", "No experiment setup set!");
        //         return;
        //     }

        //     // Set random seeds
        //     gRandom = new TRandom3(_seed);

        //     ConfigFileParser parser(configfile);
        //     ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
        //     // cfgInfo->display();
        //     ReactionInfo* reaction = cfgInfo->reactionList()[0];

        //     // AmpToolsInterface
        //     AmpToolsInterface::registerAmplitude( A() );
        //     AmpToolsInterface ATI(cfgInfo, AmpToolsInterface::kMCGeneration);
        //     line();

        //     TH1D * hw = new TH1D("weight", "weight", 150, 1, -1);
        //     plotter plotter;
        //     histogram_1D h = plotter.new_histogram_1D();
        //     h.set_labels("Intensity", "Cumulative % of Events");

        //     int nGenerated = 0;
        //     int Passes = 1, nFailed = 0;
        //     while (nGenerated < nEvents)
        //     {
        //         if (nFailed > 3)
        //         {
        //             warning("EventGenerator::weightProfile", "Three passes with no events generated, possible infinite loop? Exiting...");
        //             exit(1);
        //         };

        //         int nPass = 0;
                
        //         // First generate a batch of phase space
        //         auto batch       = generateBatch();
        //         double maxWeight = getMaxWeight(batch);
        //         int nPS = 0;  // Number of generated PS points

        //         // Do hit-or-miss to generate phasespace
        //         for (auto event : batch)
        //         {
        //             double weight = event.weight(); 
        //             if (weight / maxWeight > gRandom->Rndm())
        //             {
        //                 event.setWeight(1.); // Reset the weight to 1
        //                 ATI.loadEvent(&event, nPS, batch.size()); // Instead of writing, load to AmpTools
        //                 nPS++;
        //             };
        //         };

        //         // Instead of hit-or-miss we load all weights into a histogram
        //         double maxIntensity = ATI.processEvents(reaction->reactionName());
        //         for (int i = 0; i < nPS; i++)
        //         {
        //             if ((nGenerated + nPass) >= nEvents) break;

        //             double Intensity = ATI.intensity(i); 
        //             h.fill(Intensity);
        //             nPass++;
        //         };
        //         nGenerated += nPass;
        //         Passes++;
        //         if (nPass == 0) nFailed++;
        //         ATI.clearEvents();
        //     };

        //     h.save(outfile, true);
        //     return;
        // };

        protected:

        // Energy dependence
        Experiment * _exp = nullptr;
        double _constant_E = 9.0;

        // Final state masses
        double _masses[N];

        // Final state particle labels
        std::array<std::string,N> _labels;

        // For hit-or-miss we always generate phase-spaces in chunks of this size
        int _nBatch = 1E6;

        // Related to TGenPhaseSpace
        int _seed = 0;
        TGenPhaseSpace _generator;
    };

};

#endif