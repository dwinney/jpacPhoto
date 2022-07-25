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
    TString JPACPHOTO = gSystem->Getenv("JPACPHOTO");
    TString JPACSTYLE = gSystem->Getenv("JPACSTYLE");

    // look for and load jpacstyle
    gInterpreter->AddIncludePath(JPACSTYLE + "/include/");
    auto stylib = gSystem->Load(JPACSTYLE + "/lib/libjpacStyle." + gSystem->GetSoExt());
    if (stylib != 0) Fatal("jpacPhoto::Load","libjpacStyle not found");

    // do the same for jpacphoto
    gInterpreter->AddIncludePath(JPACPHOTO + "/include/core");
    auto photolib = gSystem->Load(JPACPHOTO + "/lib/libjpacPhoto." + gSystem->GetSoExt());
    if (photolib != 0) Fatal("jpacPhoto::Load","libjpacPhoto not found");

    // Now check for the non-essential libraries
    gInterpreter->AddIncludePath(JPACPHOTO + "/include/inclusive");
    auto inclib = gSystem->Load(JPACPHOTO + "/lib/libjpacInclusive." + gSystem->GetSoExt());
    if (inclib != 0) Warning("jpacPhoto::Load","libjpacInclusive not found");

    gInterpreter->AddIncludePath(JPACPHOTO + "/include/box");
    auto boxlib = gSystem->Load(JPACPHOTO + "/lib/libjpacBox." + gSystem->GetSoExt());
    if (boxlib != 0) Warning("jpacPhoto::Load","libjpacBox not found");
}