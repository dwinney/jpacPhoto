// Pseudo-code template for setting up a fitter for an analysis
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2026)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Instituto de Ciencias Nucleares (ICN),
//               Universidad Nacional Autonoma de México (UNAM)
// Email:        daniel.winney@nucleares.unam.mx
// ------------------------------------------------------------------------------

#ifndef FITTER_TEMPLATE_HPP
#define FITTER_TEMPLATE_HPP

#include "amplitude.hpp"
#include "utilities.hpp"

namespace jpacPhoto
{
    // Its often useful to have a namespace for the specific reaction to avoid name clashes
    namespace my_reaction
    {
        // data_set contains an member _type which is used to identify what type of data
        // it is and used by fitters to filter what to calculate 
        // The types is just an int but for human readability its useful to define constants
        const int kType1 = 0;
        const int kType2 = 1;

        // This is the actual fitter implementation
        // All the functions need to be static for it to work and thus we arent allowed nonstatic member data
        struct my_fitter
        {
            // This method takes the data_set._type and returns a string identifier of what data it is
            static std::string data_type(int i)
            {
                switch (i)
                {
                    case kType1: return "Type 1 - x vs z";
                    case kType2: return "Type 2 - y vs z";
                    default: return "ERROR";
                };
                return "ERROR";
            };

            // This method is what is actually being minimized, i.e. chi2 or log likelihood or whatever
            static double fcn(const std::vector<data_set> &data_sets, amplitude to_fit)
            {
                double sum = 0;
                for (auto data : data_sets)
                {
                    // Calculate fcn for each data, filtered by data._type if needed
                };
            };

            // These last two methods are for easily handling fit results by exporting/importing parameter
            // values from a file

            static std::vector<double> import_parameters(std::string filepath)
            {
                // Given a path to a file parse the file and arrange them in the correct order
                // expected by our amplitude
            };

            static void export_parameters(std::vector<double> pars, std::string filepath, std::string description = "")
            {
                // Given a vector of best-fit pars, export them into a file in the format that can be imported with
                // my_fitter::import_parameters above. Typically we also want to put some extra description related
                // to e.g. the fit results or the circumstances of the data, etc., in a preamble to easily identify 
                // different results. 
            };
        };
    };
};

#endif