/*
 * bin_generated_decays.cpp
 * ROOT macro to bin AmpGen generated D -> K3pi decays into predefined bins
 *
 * This script is gross and needs a lot of clearing up; only exists to test my understanding
 */
#include <cstring>
#include <iostream>
#include <string>

#include "k3pi_binning.h"

// ---- Magic Numbers
// @@@ most of these can probably be removed now that we're using "real" generated data
// @@@ some of them probably need to be changed

// Length of DalitzEventList branch names, +1 for null terminator
#define NAME_LENGTH 9

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
 * Store the px, py, pz and E data of a particle named particle Name from myTree into an array of arrays
 * Allocates memory to this array which must be freed by the caller
 *
 */
double **writeArrays(TTree *myTree, const char *particleName)
{
    double **const arrays = new double *[4] { nullptr, nullptr, nullptr, nullptr };
    double *       energies{nullptr};
    double *       xMomenta{nullptr};
    double *       yMomenta{nullptr};
    double *       zMomenta{nullptr};
    const char *   px = "_Px";
    const char *   py = "_Py";
    const char *   pz = "_Pz";
    const char *   E  = "_E";

    // Generate names of the particle's momenta and energy
    char particlePx[NAME_LENGTH];
    strcpy(particlePx, particleName);
    strcat(particlePx, px);

    char particlePy[NAME_LENGTH];
    strcpy(particlePy, particleName);
    strcat(particlePy, py);

    char particlePz[NAME_LENGTH];
    strcpy(particlePz, particleName);
    strcat(particlePz, pz);

    char particleE[NAME_LENGTH];
    strcpy(particleE, particleName);
    strcat(particleE, E);

    // Create arrays of this particle's momentum and energy
    energies = writeArray(myTree, particleE);
    xMomenta = writeArray(myTree, particlePx);
    yMomenta = writeArray(myTree, particlePy);
    zMomenta = writeArray(myTree, particlePz);

    arrays[0] = xMomenta;
    arrays[1] = yMomenta;
    arrays[2] = zMomenta;
    arrays[3] = energies;

    return arrays;
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

    double kE{0.0};
    double kpx{0.0};
    double kpy{0.0};
    double kpz{0.0};

    // Create an array of arrays pointing to the K data
    double **   kArrays{nullptr};
    const char *kParticleName = "_1_K~";
    kArrays                   = writeArrays(myTree, kParticleName);

    for (int i = 0; i < myTree->GetEntries(); ++i) {
        TLorentzVector kLorentzVector{};
        TLorentzVector pi1LorentzVector{};
        TLorentzVector pi2LorentzVector{};
        TLorentzVector pi3LorentzVector{}; //@@@ initialise these with values from the above allocated arrays

        // Use order K+ pi- pi- pi+

        //
        //      Create a vector of TLorentzVectors for this event (K+, pi-, pi-, pi+)
        //      run eventFromVectors on it to get an event vector (or could just create the even vector straight away)

        //      @@@ Work out the CF and DCS amplitudes of this event, using the bins object created above
        //      @@@ Update global hadronic parameters
        //      @@@ Find which bin the event belongs in then update its hadronic parameters
    }
    auto tmpCanvas = new TCanvas("c", "c", 600, 600);
    TH1* foo = new TH1D("foo", "bar", 100, 0.4, 1);
    foo->FillN(1000, kArrays[3], 0);
    foo->Draw();

    delete[] kArrays[0];
    delete[] kArrays[1];
    delete[] kArrays[2];
    delete[] kArrays[3];
    delete[] kArrays;

    std::cout << "==== Global: =====================" << std::endl;
    std::cout << "==================================" << std::endl;
}
