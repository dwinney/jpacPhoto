#include "Experiment.hpp"

#include <string>
#include <cstdlib>
#include "TRandom3.h"

#ifndef EXPERIMENT_GLUEX_HPP
#define EXPERIMENT_GLUEX_HPP

namespace jpacPhoto
{
    struct GlueX : public Experiment
    {
        GlueX()
        : Experiment("GlueX")
        {
            initialize();
        };

        GlueX(double x, double y)
        : Experiment("GlueX"), _emin(x), _emax(y)
        {
            initialize();
        };

        inline void set_energy_range(double x, double y)
        {
            if (_emax > 12. || _emin < 8.2) return;
            _emin = x; _emax = y;
        };

        inline double beam_energy()
        {
            double sampledE = -1;
            while(sampledE < _emin || sampledE > _emax)
            {
                double x = r3.Uniform(0,max_c);
                sampledE = gr->Eval(x);
            }

            return sampledE;
        };

        inline void initialize()
        {
            // Get beam distribution
            enefile = TFile::Open("beamE_v4.root");

            enehist = (TH1*)enefile->Get("ene");
            if (!enehist)
            {
                warning("Cannot open beam spectrum file");
            };

            hc = enehist->GetCumulative();
            int n = hc->GetNbinsX();
            double x_c[n], y_x[n];
            max_c = hc->GetMaximum();

            for(int i=0;i<n;i++) {
            x_c[i] = hc->GetBinContent(i);
            y_x[i] = hc->GetBinCenter(i);
            }
            gr = new TGraph(n,x_c,y_x);

            int j = enehist->GetNbinsX();
            double x_e[j], y_e[j];
            for (int i = 0; i < j; i++)
            {
                x_e[i] = enehist->GetBinCenter(i);
                y_e[i] = enehist->GetBinContent(i);
            }
            weightE = new TGraph(j, x_e, y_e); 
        };

        // ROOT structures
        TFile * enefile = NULL;
        TH1 * enehist = NULL;
        TH1 * hc = NULL; //cumulative
        TGraph * gr = NULL;
        TGraph * weightE = NULL;
        TRandom3 r3;
        int max_c;
        double _emin = 10.0, _emax = 10.2;
    };
};

#endif