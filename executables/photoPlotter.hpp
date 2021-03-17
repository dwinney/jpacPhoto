// ---------------------------------------------------------------------------
// Wrapper for the jpacGraph1D class to get rid of copying the same code snippets
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _PHOTOPLOT_
#define _PHOTOPLOT_

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include "amplitudes/reaction_kinematics.hpp"
#include "amplitudes/amplitude.hpp"

#include <vector>
#include <string>

namespace jpacPhoto
{
    class photoPlotter : public jpacGraph1D
    {
        public:
        
        // default constructor
        photoPlotter( std::vector<amplitude*> amps_)
        : amps(amps_)
        {};
        
        int line_width = 4;
        
        int N = 20;
        bool PRINT_TO_COMMANDLINE = true;
        bool LAB_ENERGY = false;
        double xmin, xmax, ymin, ymax;

        bool SHOW_LEGEND = false;
        double xlegend, ylegend;

        std::string xlabel, ylabel, filename;

        void Plot(std::string observable, double theta = 0.)
        {
            int obs = translate(observable);
            if ((obs < 0) || (obs > 7))
            {
                std::cout << "Error! Invalid string \"" << observable << "\" passed to photoPlotter::Plot()!";
                return;
            }
        
            for (int n = 0; n < amps.size(); n++)
            {
                std::cout << std::endl << "Printing amplitude: " << amps[n]->_identifier << "\n";
    
                double th;
                (LAB_ENERGY) ? (th = E_beam(amps[n]->_kinematics->Wth())) : (th = amps[n]->_kinematics->Wth());

                auto F = [&](double x)
                {
                    double W;
                    (LAB_ENERGY) ? (W = W_cm(x)) : (W = x);

                    double s = W*W;               
                    double t = amps[n]->_kinematics->t_man(s, theta * DEG2RAD);

                    switch (obs) 
                    {
                        case 0:
                        {
                            return amps[n]->probability_distribution(s, t);
                        }
                        case 1:
                        {
                            
                            return amps[n]->integrated_xsection(s);
                        }
                        case 2:
                        {
                            return amps[n]->differential_xsection(s, t);
                        }
                        case 3:
                        {
                            return amps[n]->A_LL(s, t);
                        }
                        case 4:
                        {
                            return amps[n]->K_LL(s, t);
                        }
                        case 5: 
                        {
                            return amps[n]->beam_asymmetry_4pi(s, t);
                        }
                        case 6: 
                        {
                            return amps[n]->beam_asymmetry_y(s, t);
                        }
                        case 7:
                        {
                            return amps[n]->parity_asymmetry(s, t);
                        }
                    };

                    return 0.;
                };

                std::array<std::vector<double>, 2> x_fx;
                if (xmin < th)
                {
                    x_fx = vec_fill(N, F, th + EPS, xmax, PRINT_TO_COMMANDLINE);
                }
                else
                {
                    x_fx = vec_fill(N, F, xmin, xmax, PRINT_TO_COMMANDLINE);
                }

                AddEntry(x_fx[0], x_fx[1], amps[n]->_identifier);
            }

            SetXaxis(xlabel, xmin, xmax);
            SetYaxis(ylabel, ymin, ymax);

            if (SHOW_LEGEND == false)
            {
                SetLegend(false);
            }
            else
            {
                SetLegend(xlegend, ylegend);
            }

            // Output to file
            jpacGraph1D::Plot(filename);
        };

        private:
        std::vector<amplitude*> amps;

        int translate(std::string observable)
        {
            if      (observable == "probability_distribution")  return 0;
            else if (observable == "integrated_xsection")       return 1;
            else if (observable == "differential_xsection")     return 2;
            else if (observable == "A_LL")                      return 3;
            else if (observable == "K_LL")                      return 4;
            else if (observable == "beam_asymmetry_4pi")        return 5;
            else if (observable == "beam_asymmetry_y")          return 6;
            else if (observable == "parity_asymmetry")          return 7;
            else return -1;
        };
    }; 
};

#endif