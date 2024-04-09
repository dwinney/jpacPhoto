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

    // Version with two particles int the final state
    class TwoBodyWriter : public DataWriter<2>
    {
        public:

        TwoBodyWriter( std::array<std::string,2> labels, const std::string & outfile )
        : DataWriter(labels)
        {
            initialize(outfile);
        };

        protected: 

        // Only extra variables are invariants
        inline void setupExtras()
        {
            _tree->Branch( "s", &_s, "s/D");
            _tree->Branch( "t", &_t, "t/D");
        };

        inline void calculateExtras(std::vector<TLorentzVector> fvecs)
        {
            // Invariant total energy
            _s = (fvecs[0] + fvecs[1]).M2();

            // Invariant momentum transfer
            TLorentzVector q(0., 0., _Ebeam, _Ebeam);
            _t = (q - fvecs[0]).M2();
        };

        // Mandelstam Variables
        double _s, _t;
    };

    // Same thing but with three
    class ThreeBodyWriter : public DataWriter<3>
    {
        public:

        ThreeBodyWriter( std::array<std::string,3> labels, const std::string & outfile )
        : DataWriter(labels)
        {
            initialize(outfile);
        };

        protected:

        // Calculate five invariants
        inline void setupExtras()
        {
            _tree->Branch( "s",    &_s, "s/D" );
            _tree->Branch( "t",    &_t, "t/D" );
            _tree->Branch( subscript("s", 0, 1).c_str(), &_s12, subscriptD("s", 0, 1).c_str());
            _tree->Branch( subscript("s", 0, 2).c_str(), &_s1,  subscriptD("s", 0, 2).c_str());
            _tree->Branch( subscript("s", 1, 2).c_str(), &_s2,  subscriptD("s", 1, 2).c_str());
            _tree->Branch( subscript("t", 0).c_str(), &_t1,  subscriptD("t", 0).c_str());
            _tree->Branch( subscript("t", 1).c_str(), &_t2,  subscriptD("t", 1).c_str());

            _tree->Branch( "cosGJ",   &_cosGJ,   "cosGJ/D");
            _tree->Branch( "thetaGJ", &_thetaGJ, "thetaGJ/D");
            _tree->Branch( "phiGJ",   &_phiGJ,   "phiGJ/D");
        };

        inline void calculateExtras(std::vector<TLorentzVector> fvecs)
        {
            // Inital state reconstructed from the internally calcualted _Ebeam
            TLorentzVector beam(  0., 0., _Ebeam, _Ebeam);
            TLorentzVector target(0., 0., 0.,   M_PROTON);
            TLorentzVector m1(fvecs[0]);
            TLorentzVector m2(fvecs[1]);
            TLorentzVector recoil(fvecs[2]);

            _s   = (beam + target).M2();
            _t   = (recoil - target).M2();
            _s12 = (m1 + m2).M2();

            // Subsystem invariant masses
            _s1  = (m1 + recoil).M2();
            _s2  = (m2 + recoil).M2();

            // Subsystem momentum transfers
            _t1 = (beam - m1).M2();
            _t2 = (beam - m2).M2();

            // Boost only the initialstate and particle 1 to GJ frame to calculate angles
            beam.Boost(  - (m1 + m2).BoostVector());
            recoil.Boost(- (m1 + m2).BoostVector());
            m1.Boost(    - (m1 + m2).BoostVector());

            TVector3 z =  beam.Vect().Unit();
            TVector3 y = (beam.Vect().Cross(recoil.Vect())).Unit();
            TVector3 x = y.Cross(z);

            TVector3 GJ( m1.Vect().Dot(x), m1.Vect().Dot(y), m1.Vect().Dot(z));
            _thetaGJ = GJ.Theta(); _phiGJ = GJ.Phi(); _cosGJ = cos(_thetaGJ);
        };

        inline std::string subscript(std::string x, int i, int j = -1)
        {
            return  (j == -1) ? x + "_" + _labels[i] : x + "_" + _labels[i] + _labels[j];
        };
        inline std::string subscriptD(std::string x, int i, int j = -1)
        {
            return subscript(x, i, j) + "/D";
        };

        double _s, _t, _s12, _s1, _s2, _t1, _t2, _cosGJ, _thetaGJ, _phiGJ;
    };
}

#endif