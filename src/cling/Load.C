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
    TString LIB_EXT = gSystem->GetSoExt();

    //----------------------------------------------------------------------
    // Core physics library

    TString JPACPHOTO_DIR  = gSystem->Getenv("JPACPHOTO");

    TString JPACPHOTO_INCCORE  = JPACPHOTO_DIR;
            JPACPHOTO_INCCORE += "/src/core";
    TString JPACPHOTO_INCAMP  = JPACPHOTO_DIR;
            JPACPHOTO_INCAMP += "/amplitudes";
            
    TString JPACPHOTO_LIB  = JPACPHOTO_DIR;
            JPACPHOTO_LIB += "/lib/libJPACPHOTO.";
            JPACPHOTO_LIB += LIB_EXT;

    if (!gSystem->AccessPathName(JPACPHOTO_LIB.Data()))
    {
        gInterpreter->AddIncludePath( JPACPHOTO_INCCORE.Data());
        gInterpreter->AddIncludePath( JPACPHOTO_INCAMP.Data());

        Int_t pholib = gSystem->Load( JPACPHOTO_LIB.Data());
    }
    else
    {
        Warning("jpacPhoto::Load", "jpacPhoto library not found! Path given: %s", JPACPHOTO_LIB.Data());
    }

    //----------------------------------------------------------------------
    // Plotting library

    TString JPACSTYLE_DIR  = gSystem->Getenv("JPACSTYLE");

    TString JPACSTYLE_INC  = JPACSTYLE_DIR;
            JPACSTYLE_INC += "/include/";
            
    TString JPACSTYLE_LIB  = JPACSTYLE_DIR;
            JPACSTYLE_LIB += "/lib/libjpacStyle.";
            JPACSTYLE_LIB += LIB_EXT;

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