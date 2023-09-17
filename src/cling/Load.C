// Load script that links all the required libraries.
// Adapted from the installation of elSpectro
// [https://github.com/dglazier/elSpectro]
// by Derek Glazier and others
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@alumni.iu.edu
// -----------------------------------------------------------------------------

void Load()
{
    TString lib_ext = gSystem->GetSoExt();

    //----------------------------------------------------------------------
    // Core physics library

    TString main_dir  = gSystem->Getenv("JPACPHOTO");

    // Load the main library files
    TString core        = main_dir + "/src/core"; 
    TString main_lib    = main_dir + "/lib/libJPACPHOTO." + lib_ext;

    // Supplementary header files 
    TString physics  = main_dir + "/physics";
    TString data        = main_dir + "/data";

    if (!gSystem->AccessPathName(main_lib.Data()))
    {
        Int_t pholib = gSystem->Load( main_lib.Data());
        gInterpreter->AddIncludePath( core.Data());
        gInterpreter->AddIncludePath( physics.Data());
        gInterpreter->AddIncludePath( data.Data());
    }
    else
    {
        Warning("jpacPhoto::Load", "jpacPhoto library not found! Path given: %s", main_lib.Data());
    }
}