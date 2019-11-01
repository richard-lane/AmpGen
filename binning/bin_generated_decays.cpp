/*
 * bin_generated_decays.cpp
 * ROOT macro to bin AmpGen generated D -> K3pi decays into predefined bins
 */
#include <string>

#include "k3pi_binning.h"

// ---- Magic Numbers
// @@@ most of these can probably be removed now that we're using "real" generated data
// Decaying particle (Momentum, Energy) vector in GeV

// DCS and CF relative amplitude and phase
#define DCS_MAGNITUDE 0.0601387
#define DCS_PHASE 1.04827

/// Bin limits in phase, centred on zero by construction
#define NUM_BINS 5
#define BIN_LIMITS -39, 0, 43, 180

inline bool fileExists(const std::string &name)
{
    ifstream f(name.c_str());
    return f.good();
}

/*
 * Store the data from myBranch on myTree into myArray
 * Returns an array of doubles that must be freed by the caller.
 */
double *writeArray(TTree *myTree, const char *myBranchName)
{

    const long long numEntries{myTree->GetEntries()};
    double          myData{0.0};

    // Point the desired branch at the myData variable.
    myTree->SetBranchAddress(myBranchName, &myData);

    // Create an array of the appropriate size to store this data
    double *const myDataArray = new double[numEntries];

    for (int i = 0, N = myTree->GetEntries(); i < N; ++i) {
        myTree->GetEntry(i);
        myDataArray[i] = myData;
    }

    return myDataArray;
}

/*
 * Bin the decays modelled an AmpGen generated inputFile into phase bins as defined by $BIN_LIMITS
 *
 */
void bin_generated_decays(TFile *inputFile)
{
    // Apply scaling and rotation to the DCS amplitude such that we get dcs/cf amplitude ratio 'r' = 0.055
    // and the average relative strong-phase between the two amplitudes ~ 0.
    std::complex<double> dcs_offset = DCS_MAGNITUDE * exp(std::complex<double>(0, 1) * DCS_PHASE * M_PI / 180.);

    // Define the bins based on the form of the DCS and CF decays
    const std::string dcsFile{"binning/dcs.so"};
    const std::string cfFile{"binning/cf.so"};

    if (!fileExists(dcsFile)) {
        std::cout << "File '" << dcsFile << "' not found." << std::endl;
        throw 1;
    }
    if (!fileExists(cfFile)) {
        std::cout << "File '" << cfFile << "' not found." << std::endl;
        throw 1;
    }

    k3pi_binning::binning bins(dcsFile, cfFile, dcs_offset, {BIN_LIMITS});

    /// calculate global hadronic parameters, and parameters in each of the bins.
    std::complex<double>              z(0, 0);
    double                            n_cf(0);
    double                            n_dcs(0);
    std::vector<std::complex<double>> z_binned(NUM_BINS, std::complex<double>(0, 0));
    std::vector<double>               n_cf_binned(NUM_BINS, 0);
    std::vector<double>               n_dcs_binned(NUM_BINS, 0);

    // Read in the tree and branches from the provided ROOT file
    // @@@ for now just read in the k energy
    TTree *     myTree   = nullptr;
    const char *treeName = "DalitzEventList"; // The name of the tree of interest in the ROOT file
    inputFile->GetObject(treeName, myTree);

    double      kE{0.0};
    double      kpx{0.0};
    double      kpy{0.0};
    double      kpz{0.0};

    // Could improve this to not loop over each array separately, but i dont want to do that yet
    const char *kEnergyBranchName = "_1_K~_E";
    const char *kpxBranchName     = "_1_K~_Px";
    const char *kpyBranchName     = "_1_K~_Py";
    const char *kpzBranchName     = "_1_K~_Pz";
    myTree->SetBranchAddress(kEnergyBranchName, &kE);
    myTree->SetBranchAddress(kpxBranchName, &kpx);
    myTree->SetBranchAddress(kpyBranchName, &kpy);
    myTree->SetBranchAddress(kpzBranchName, &kpz);

    const char *pi1EnergyBranchName = "_1_K~_E";
    const char *pi1pxBranchName = "_1_K~_Px";
    const char *pi1pyBranchName = "_1_K~_Py";
    const char *pi1pzBranchName = "_1_K~_Pz";
    myTree->SetBranchAddress(pi1EnergyBranchName, ); // @@@ fill these in 
    myTree->SetBranchAddress(pi1pxBranchName, );
    myTree->SetBranchAddress(pi1pyBranchName, );
    myTree->SetBranchAddress(pi1pzBranchName, );

    double *kEArray  = nullptr;
    double *kpxArray = nullptr;
    double *kpyArray = nullptr;
    double *kpzArray = nullptr;
    kEArray          = writeArray(myTree, kEnergyBranchName);
    kpxArray         = writeArray(myTree, kpxBranchName);
    kpyArray         = writeArray(myTree, kpyBranchName);
    kpzArray         = writeArray(myTree, kpzBranchName);

    for (int i = 0; i < myTree->GetEntries(); ++i) {
        TLorentzVector kLorentzVector{};
        TLorentzVector pi1LorentzVector{};
        TLorentzVector pi2LorentzVector{};
        TLorentzVector pi3LorentzVector{}; //@@@ initialise these with values from the above allocated arrays

        // Use order K+ pi- pi- pi+
        std::vector<double> event =
            k3pi_binning::eventFromVectors(kLorentzVector, pi1LorentzVector, pi2LorentzVector, pi3LorentzVector);
        //
        //      Create a vector of TLorentzVectors for this event (K+, pi-, pi-, pi+)
        //      run eventFromVectors on it to get an event vector (or could just create the even vector straight away)

        //      @@@ Work out the CF and DCS amplitudes of this event, using the bins object created above
        //      @@@ Update global hadronic parameters
        //      @@@ Find which bin the event belongs in then update its hadronic parameters
    }

    delete[] kEArray;
    delete[] kpxArray;
    delete[] kpyArray;
    delete[] kpzArray;

    std::cout << "==== Global: =====================" << std::endl;
    std::cout << "==================================" << std::endl;
}
