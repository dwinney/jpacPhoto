// Reworking of the jpacStyle as internal library special for jpacPhoto
// This defines the interface object which creates plots and defines
// the common settings for the JPC style
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef PLOTTER_HPP
#define PLOTTER_HPP

#include "plot.hpp"
#include "histogram.hpp"
#include "colors.hpp"

#include <array>

#include <TROOT.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TStyle.h>
#include <TError.h>
#include <TLatex.h>
#include <string>

namespace jpacPhoto
{    
    // The plotter object is the basic structure for setting up a plot in ROOT 
    // using the JPAC style.
    class plotter
    {
        public: 

        // Basic constructor
        // initializes all the global things (colors and style settings)
        plotter()
        {
            initialize_style();
            initialize_colors();
        };

        // Create a new plot!
        inline plot new_plot(std::string file = "")
        {
            // Make it the global default style
            gROOT->SetStyle("jpacStyle");

            // Set up a new canvas 
            _Nplots++;
            std::string name = "c" + std::to_string(_Nplots);
            TCanvas *canvas = new TCanvas(name.c_str(), name.c_str(), 600, 600);
            canvas->UseCurrentStyle();
            canvas->SetTopMargin(0.05);
            canvas->SetRightMargin(0.03);
            canvas->SetLeftMargin(0.16);
            canvas->SetBottomMargin(0.12);
            canvas->SetFixedAspectRatio();

            // if no filename is passes, choose a default one
            if (file == "") file += "plot" + std::to_string(_Nplots) + ".pdf";
            
            // Pass this to the plot which will populate it with everything else
            return plot(canvas, file);
        };

        // Create a new histogram!
        inline histogram_1D new_histogram_1D(std::string label = "")
        {
            // Make it the global default style
            gROOT->SetStyle("jpacStyle");

            // Set up a new canvas 
            _Nplots++;
            std::string name = "c" + std::to_string(_Nplots);
            TCanvas *canvas = new TCanvas(name.c_str(), name.c_str(), 600, 600);
            canvas->UseCurrentStyle();
            canvas->SetTopMargin(0.05);
            canvas->SetRightMargin(0.03);
            canvas->SetLeftMargin(0.12);
            canvas->SetBottomMargin(0.12);
            canvas->SetFixedAspectRatio();
            
            // Pass this to the plot which will populate it with everything else
            return histogram_1D(canvas, label);
        };

        // Create a new histogram!
        inline histogram_2D new_histogram_2D(std::string labelx = "", std::string labely = "")
        {
            // Make it the global default style
            gROOT->SetStyle("jpacStyle");

            // Set up a new canvas 
            _Nplots++;
            std::string name = "c" + std::to_string(_Nplots);
            TCanvas *canvas = new TCanvas(name.c_str(), name.c_str(), 600, 600);
            canvas->UseCurrentStyle();
            canvas->SetTopMargin(0.05);
            canvas->SetRightMargin(0.14);
            canvas->SetLeftMargin(0.14);
            canvas->SetBottomMargin(0.12);
            canvas->SetFixedAspectRatio();
            
            // Pass this to the plot which will populate it with everything else
            return histogram_2D(canvas, {labelx, labely});
        };

        // Combine a vector of plots into a grid of specified dimensions
        static void combine(std::array<int,2> dims, std::vector<plot> plots, std::string filename);
        static void combine(std::array<int,2> dims, std::vector<combinable*> hists, std::string filename);

        static void stack(std::vector<plot> plots, std::string filename);

        private:

        // "Global" style settings for all plots
        TStyle *_style = new TStyle("jpacStyle", "JPAC Style");

        // JPAC Logo to put on corner
        TLatex *_logo       = NULL; 
        const int kjpacFont = 132;

        // Initialize each color with a ROOT colorID
        void initialize_colors();

        // Actual color objects to define a palette
        TColor *jpacBlue,   *jpacRed,    *jpacGreen;
        TColor *jpacOrange, *jpacPurple, *jpacBrown;
        TColor *jpacPink,   *jpacGold,   *jpacAqua;
        TColor *jpacGrey,   *jpacDarkGrey;

        // Number of plots created
        int _Nplots = 0;

        // Initialize all the appropriate settings into _style
        void initialize_style();
    };
};

#endif