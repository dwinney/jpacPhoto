// Abstract data writer class for an N particle final state
// This will read the the 4-vectors and weights of each of the particles
// Additional quantities can be read by overriding the setupExtras() and
// calculateExtras() in a derived class
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU),
//               University of Bonn (HISKP)
// Email:        winney@hiskp.uni-bonn.de
// ------------------------------------------------------------------------------

#ifndef DATA_WRITER_HPP
#define DATA_WRITER_HPP

#include <array>
#include <cassert>
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TSystem.h"
#include "IUAmpTools/Kinematics.h"
#include "constants.hpp"
#include "print.hpp"

namespace jpacPhoto
{
    template<int N>
    class DataWriter 
    {
        public:
        
        // Constructor doesnt run initialize root structure here
        // Instead initialize() call needed in derived class constructor
        DataWriter<N>( std::array<std::string,N> labels )
        : _labels(labels)
        {};

        // If this is the actuall class to be used, then call initialize
        DataWriter<N>( std::array<std::string,N> labels, const std::string & outfile )
        : _labels(labels)
        {
            initialize(outfile);
        };

        // Destructor writes to file
        ~DataWriter<N>()
        {
            _tree->AutoSave();
            delete _file;
        };

        inline void writeEvent( const Kinematics& kin)
        {
            std::vector<TLorentzVector> particleList = kin.particleList();
            
            TLorentzVector q;
            for (int i = 0; i < N; i++)
            {
                _E[i]  = particleList[i].E();
                _Px[i] = particleList[i].Px();
                _Py[i] = particleList[i].Py();
                _Pz[i] = particleList[i].Pz();

                q += particleList[i];
            };

            _Ebeam  = jpacPhoto::E_beam( q.M() );
            _weight = kin.weight();
            calculateExtras(particleList);

            _tree->Fill();
            _nEvents++;
        };

        inline const int eventCounter(){ return _nEvents; };

        protected:

        // Set up tree and default branches
        inline void initialize(const std::string & outfile)
        {
            TH1::AddDirectory( kFALSE );
            gSystem->Load( "libTree" );

            _file = new TFile( outfile.c_str(), "recreate" );
            _tree = new TTree( "nt", "nt" );
            _nEvents = 0;

            _tree->Branch( "E_beam", &_Ebeam,  "E_beam/D" );
            _tree->Branch( "weight", &_weight, "weight/D" );

            // Every writer will save the elements of the four vector automatically
            for (int i = 0; i < N; i++)
            {
                _tree->Branch( ("E_" +_labels[i]).c_str(), &_E[i],  ("E_" +_labels[i]+"/D").c_str());
                _tree->Branch( ("Px_"+_labels[i]).c_str(), &_Px[i], ("Px_"+_labels[i]+"/D").c_str());
                _tree->Branch( ("Py_"+_labels[i]).c_str(), &_Py[i], ("Py_"+_labels[i]+"/D").c_str());
                _tree->Branch( ("Pz_"+_labels[i]).c_str(), &_Pz[i], ("Pz_"+_labels[i]+"/D").c_str());
            }

            // Any additional saved info must be user defined through overriding this function
            setupExtras();
        };

        // Each implentation can allow additional branches to be saves to the file
        // This function is called at the end of initialize and does nothing by default
        virtual void setupExtras(){};

        // At each write event any extra quantities not the raw 4-vectors
        // need to be calcualted through this function
        virtual void calculateExtras(std::vector<TLorentzVector> fvecs){};

        // Particle labels to appear in branches
        std::array<std::string, N> _labels;

        // (Final state) Four vector elements
        std::array<double, N> _E, _Px, _Py, _Pz;

        // Beam energy
        double _Ebeam;

        // Incase using weighted events
        double _weight;

        // Root structures
        TFile* _file;
        TTree* _tree;
        
        int _nEvents = 0;
    };
}

#endif