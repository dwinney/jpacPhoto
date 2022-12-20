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
    TString fitter      = main_dir + "/src/fitter";
    TString main_lib    = main_dir + "/lib/libJPACPHOTO." + lib_ext;

    // Supplementary header files 
    TString amplitudes  = main_dir + "/amplitudes";
    TString data        = main_dir + "/data";

    if (!gSystem->AccessPathName(main_lib.Data()))
    {
        Int_t pholib = gSystem->Load( main_lib.Data());
        gInterpreter->AddIncludePath( core.Data());
        gInterpreter->AddIncludePath( fitter.Data());
        gInterpreter->AddIncludePath( amplitudes.Data());
        gInterpreter->AddIncludePath( data.Data());
    }
    else
    {
        Warning("jpacPhoto::Load", "jpacPhoto library not found! Path given: %s", main_lib.Data());
    }

    //----------------------------------------------------------------------
    // Plotting library

    TString JPACSTYLE_DIR  = gSystem->Getenv("JPACSTYLE");

    TString JPACSTYLE_INC  = JPACSTYLE_DIR;
            JPACSTYLE_INC += "/include/";
            
    TString JPACSTYLE_LIB  = JPACSTYLE_DIR;
            JPACSTYLE_LIB += "/lib/libjpacStyle.";
            JPACSTYLE_LIB += lib_ext;

    if (!gSystem->AccessPathName(JPACSTYLE_LIB.Data()))
    {
        gInterpreter->AddIncludePath( JPACSTYLE_INC.Data());
        Int_t stylib = gSystem->Load( JPACSTYLE_LIB.Data());
    }
    else
    {
        Warning("jpacPhoto::Load", "jpacStyle library not found! Path given: %s", JPACSTYLE_LIB.Data());
    }


    // //----------------------------------------------------------------------
    // // Non-essential libraries

    // TString INCLUSIVE_INC  = JPACPHOTO_DIR;
    //         INCLUSIVE_INC += "/include/inclusive";
    // TString INCLUSIVE_LIB  = JPACPHOTO_DIR;
    //         INCLUSIVE_LIB += "/lib/libjpacInclusive.";
    //         INCLUSIVE_LIB += LIB_EXT;

    // if (!gSystem->AccessPathName(INCLUSIVE_LIB.Data()))
    // {
    //     gInterpreter->AddIncludePath( INCLUSIVE_INC.Data());
    //     Int_t inclib = gSystem->Load( INCLUSIVE_LIB.Data());
    // }

    // TString BOX_INC  = JPACPHOTO_DIR;
    //         BOX_INC += "/include/box";
    // TString BOX_LIB  = JPACPHOTO_DIR;
    //         BOX_LIB += "/lib/libjpacBox.";
    //         BOX_LIB += LIB_EXT;

    // if (!gSystem->AccessPathName(BOX_LIB.Data()))
    // {
    //     gInterpreter->AddIncludePath( BOX_INC.Data());
    //     Int_t boxlib = gSystem->Load( BOX_LIB.Data());
    // }
}