{
    printf("Running AmpGen rootlogon script\n");
    gSystem->Load("~/Documents/HEP/AmpGen/build/lib/libAmpGen.so");
    gROOT->ProcessLine(".include ~/Documents/Hep/Ampgen");
}

