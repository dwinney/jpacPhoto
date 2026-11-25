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
    TString main_lib    = main_dir + "/lib/libJPACPHOTO." + lib_ext;
    TString core        = main_dir + "/src"; 

    // Supplementary header files 
    TString physics  = main_dir + "/physics";
    TString data     = main_dir + "/analysis";
    
    if (!gSystem->AccessPathName(main_lib.Data()))
    {
        gSystem->Load(main_lib.Data());
        gInterpreter->AddIncludePath(core.Data());
        gInterpreter->AddIncludePath(physics.Data());
        gInterpreter->AddIncludePath(data.Data());
    }
    else
    {
        Warning("jpacPhoto::Load", "jpacPhoto library not found! Path given: %s", main_lib.Data());
    };

    TString plot_dir = gSystem->Getenv("JPACPLOTTER");
    TString plot_lib = plot_dir+"/lib/libJPACPLOTTER."+lib_ext;
    TString plot_inc = plot_dir+"/src";

    if (!gSystem->AccessPathName(plot_lib.Data()))
    {
        gSystem->Load(plot_lib.Data());
        gInterpreter->AddIncludePath(plot_inc.Data());
    };

    //----------------------------------------------------------------------
    // If Eigen is available add it
    TString eigen = gSystem->Getenv("EIGEN");
    gInterpreter->AddIncludePath(eigen.Data());

    //----------------------------------------------------------------------
    // IF dynamic AmpTools is found in lib, assume we want to load it too

    TString AmpTools_dir  = gSystem->Getenv("AMPTOOLS");
    TString AmpTools_lib  = main_dir + "/lib/libAMPTOOLS." + lib_ext;
    TString AmpTools_inc = main_dir + "/AmpTools";

    if (AmpTools_dir != "" && !gSystem->AccessPathName(AmpTools_lib.Data()))
    {
        gSystem->Load( AmpTools_lib.Data() );
        gInterpreter->AddIncludePath(AmpTools_inc.Data() );
        gInterpreter->AddIncludePath(AmpTools_dir.Data() );
    };
    
    //----------------------------------------------------------------------
    // Often we want to add amplitude files in another directory
    // Here we'll simply look for the environment variable 'JPACPHOTO_WORKING_DIR'
    TString extra = gSystem->Getenv("JPACPHOTO_WORKING_DIR");
    gInterpreter->AddIncludePath(extra.Data());
};