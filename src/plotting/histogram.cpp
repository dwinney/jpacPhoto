// Class to generate histograms with the JPAC style
// 
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               University of Bonn (HISKP)
// Email:        winney@hiskp.uni-bonn.de
// ------------------------------------------------------------------------------

#include "histogram.hpp"

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    // Draw the histogram onto whatver is the active canvas
    void histogram_1D::draw()
    {       
        // Remove title and add axis labels
        std::string labels = ";" + _xlabel + ";" + _ylabel;
        _hist->SetTitle(labels.c_str());
        
        // Center the labels
        _hist->GetXaxis()->CenterTitle(true);
        _hist->GetYaxis()->CenterTitle(true);
        _hist->GetYaxis()->SetMaxDigits(3);

        // Increase plot range by 20% to fit the logo lol
        double ymax = _hist->GetMaximum();
        _hist->SetMaximum(ymax * 1.15);

        // Draw
        _hist->Draw();

        if (_customranges)
        {
            _hist->GetXaxis()->SetRangeUser(   _xbounds[0], _xbounds[1]);
            gPad->Modified();
        };

        add_logo();
    };

    void histogram_1D::combine_draw(double x)
    {
        // Apply global settings to the pad
        gPad->UseCurrentStyle();
        gPad->SetFixedAspectRatio();
        gPad->SetTopMargin(0.05);
        gPad->SetRightMargin(0.03);
        gPad->SetLeftMargin(0.11);
        gPad->SetBottomMargin(0.12);

        // Apply the linewidth 
        draw();
    };

    // ---------------------------------------------------------------------------
    // Outputting function, generates plot and saves it to file
    void histogram_1D::save(std::string filename)
    {       
        _canvas->cd();
        // Draw the canvas
        draw();
        _canvas->Draw();

        // and print to file
        _canvas->Print(filename.c_str());
    };

    // ---------------------------------------------------------------------------
    // Draw the histogram onto whatver is the active canvas
    void histogram_2D::draw()
    {       
        // Remove title and add axis labels
        std::string labels = ";" + _xlabel + ";" + _ylabel;
        _hist->SetTitle(labels.c_str());
        
        // Center the labels
        _hist->GetXaxis()->CenterTitle(true);
        _hist->GetYaxis()->CenterTitle(true);

        // Draw
        _hist->Draw("COLZ");

        if (_customranges)
        {
            _hist->GetXaxis()->SetRangeUser(   _xbounds[0], _xbounds[1]);
            _hist->GetYaxis()->SetRangeUser(   _ybounds[0], _ybounds[1]);
            gPad->Modified();
        };
    };

    void histogram_2D::combine_draw(double x)
    {
        // Apply global settings to the pad
        _hist->SetTitleOffset(1., "y");
        gPad->UseCurrentStyle();
        gPad->SetFixedAspectRatio();
        gPad->SetTopMargin(0.05);
        gPad->SetRightMargin(0.14);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);

        draw();
    };

    // ---------------------------------------------------------------------------
    // Outputting function, generates plot and saves it to file
    void histogram_2D::save(std::string filename)
    {       
        _canvas->cd();
        // Draw the canvas
        draw();
        _canvas->Draw();

        // and print to file
        _canvas->Print(filename.c_str());
    };
};