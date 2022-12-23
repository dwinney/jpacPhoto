// Reworking of the jpacStyle as internal library special for jpacPhoto
// This defines a single plot object, which is what actually prints out the files
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#include "plot.hpp"
#include "colors.hpp"

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    // Outputting function, generates plot and saves it to file
    void plot::save()
    {
        if (_entries.empty() ) 
        {
            warning("plot::save()", "No entries added! Returning...");
            return;
        };
        
        // Logscale settings are global canvas settings so do that first
        _canvas->SetLogx(_xlog); _canvas->SetLogy(_ylog);

        // Set up the axes by grabbing them from the first entry
        auto first = _entries[0]._graph;
        first->SetTitle("");

        TAxis *xaxis = first->GetXaxis();
        xaxis->SetTitle(_xlabel.c_str());
        xaxis->CenterTitle(true);

        TAxis *yaxis = first->GetYaxis();
        yaxis->SetTitle(_ylabel.c_str());
        yaxis->CenterTitle(true);
        if (_customranges)
        {
            xaxis->SetLimits(   _xbounds[0], _xbounds[1]);
            yaxis->SetRangeUser(_ybounds[0], _ybounds[1]);
        };
        
        if (_entries[0]._isdata) { first->Draw("AP"); }
        else                     { first->Draw("AP"); };

        // Draw the first entry
        // need to parse if its a curve or data points

        // Set up legend
        _legendxoffset  = 0.3;
        _legendyoffset  = 0.036*(_Nlegend + _addheader);
        auto legend = new TLegend(_legendxcord, _legendycord, 
                                  _legendxcord + _legendxoffset, 
                                  _legendycord + _legendyoffset);
        legend->SetFillStyle(0); // Make legend transparent

        // IF theres a custom header add it
        if (_addheader) legend->SetHeader(("    " + _header).c_str());

        for (auto entry : _entries)
        {
            std::string legend_option;
            if (entry._isdata) { entry._graph->Draw("P");    legend_option = "P"; }
            else               { entry._graph->Draw("same"); legend_option = "L"; };

            if (entry._style._add_to_legend)
            {
                legend->AddEntry(entry._graph, entry._style._label.c_str(), legend_option.c_str());
            };
        };

        add_logo();

        if (_addlegend) legend->Draw();

        // Draw the canvas
        _canvas->Draw();

        // and print to file
        _canvas->Print(_filename.c_str());
    };

    // ---------------------------------------------------------------------------
    // Convert data_set and amplitude easily into plot_entries

    void plot::add_data(data_set data)
    {
        double *x, *y, *xl, *xh, *yl, *yh;
        switch (data._type)
        {
            case integrated_data: 
            {
                x  = &(data._w[0]);         y  = &(data._obs[0]);
                xl = &(data._werr[0][0]);   xh = &(data._werr[1][0]);
                yl = &(data._obserr[0][0]); yh = &(data._obserr[1][0]);
                break;
            };
            case differential_data: 
            {
                x  = &(data._t[0]);         y  = &(data._obs[0]);
                xl = &(data._terr[0][0]);   xh = &(data._terr[1][0]);
                yl = &(data._obserr[0][0]); yh = &(data._obserr[1][0]);
                break;
            };
            default: return;
        };

        TGraph *graph = new TGraphAsymmErrors(data._N, x, y, xl, xh, yl, yh);

        entry_style style;
        style._label = data._id;
        style._style = 20 + _Ndata;
        style._color = jpacColor::DarkGrey;

        _Ndata++;
        _Nlegend++;

        _entries.push_back(plot_entry(graph, style, true));
    };

    // -----------------------------------------------------------------------
    // Add different curves from amplitudes / functions

    // Add a curve from raw vectors of (x, fx) pairs
    // this is the most reminicent of jpacStyle
    // color is automatically cycled through the jpacColors
    void plot::add_curve(std::vector<double> x, std::vector<double> fx, entry_style style)
    {
        if (style._add_to_legend) _Nlegend++;
        TGraph *g = new TGraph(x.size(), &(x[0]), &(fx[0]));
        _entries.push_back(plot_entry(g, style, false));
    };

    void plot::add_curve(std::vector<double> x, std::vector<double> fx, std::string id)
    {
        _Ncurve++;
        entry_style style;
        style._color = JPACCOLORS[_Ncurve];
        style._style = kSolid;
        style._label = id;
        add_curve(x, fx, style);
    };

    // Take in a lambda an evaluation range to get the vectors
    void plot::add_curve(int N, std::array<double,2> bounds, std::function<double(double)> F, entry_style style)
    {
        double step = (bounds[1] - bounds[0]) / double(N);

        std::vector<double> x, fx;
        for (int n = 0; n < N; n++)
        {
            double xs  = bounds[0] + double(n) * (bounds[1] - bounds[0]) / double(N-1);
            double fxs = F(xs);

            x.push_back(xs);
            fx.push_back(fxs);
        };

        add_curve(x, fx, style);
    };

    void plot::add_curve(int N, std::array<double,2> bounds, std::function<double(double)> F, std::string id)
    {
        entry_style style;
        style._color = JPACCOLORS[_Ncurve];
        style._style = kSolid;
        style._label = id;
        style._add_to_legend = true;

        add_curve(N, bounds, F, style);
    };
    
    void plot::add_dashed(std::vector<double> x, std::vector<double> fx)
    {
        entry_style style;
        style._color = JPACCOLORS[_Ncurve];
        style._style = kDashed;
        style._add_to_legend = false;

        TGraph *g = new TGraph(x.size(), &(x[0]), &(fx[0]));
        _entries.push_back(plot_entry(g, style, false));
    };

    void plot::add_dashed(int N, std::array<double,2> bounds, std::function<double(double)> F)
    {
        double step = (bounds[1] - bounds[0]) / double(N);

        std::vector<double> x, fx;
        for (int n = 0; n < N; n++)
        {
            double xs  = bounds[0] + double(n) * (bounds[1] - bounds[0]) / double(N-1);
            double fxs = F(xs);

            x.push_back(xs);
            fx.push_back(fxs);
        };

        add_dashed(x, fx);
    };
};