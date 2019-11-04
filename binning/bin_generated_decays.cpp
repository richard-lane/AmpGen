/*
 * bin_generated_decays.cpp
 * ROOT macro to bin AmpGen generated D -> K3pi decays into predefined bins
 *
 * This script is gross and needs a lot of clearing up; only exists to test my understanding
 */
#include <cmath>
#include <cstring>
#include <iostream>
#include <string>

#include "k3pi_binning.h"

// ---- Magic Numbers
// @@@ most of these can probably be removed now that we're using "real" generated data
// @@@ some of them probably need to be changed

// Length of DalitzEventList branch names, +1 for null terminator
#define NAME_LENGTH 10

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
 * CoM energy squared of the ab system
 * Takes an array of arrays of kinematic particle data
 *
 * Assumes particle data arrays take the form *(Px, Py, Pz, E)
 *
 * Allocates memory to the array of s values which must be freed by the caller
 *
 */
double *s(double **particleA, double **particleB, const unsigned int length)
{
    // Initialise an array of CoM energies
    double *const sValues = new double[length];

    for (int i = 0; i < length; ++i) {
        sValues[i] = std::pow(particleA[3][i], 2) - std::pow(particleA[0][i], 2) - std::pow(particleA[1][i], 2) -
                     std::pow(particleA[2][i], 2) + std::pow(particleB[3][i], 2) - std::pow(particleB[0][i], 2) -
                     std::pow(particleB[1][i], 2) - std::pow(particleB[2][i], 2) +
                     2 * particleA[3][i] * particleB[3][i] - 2 * particleA[0][i] * particleB[0][i] -
                     2 * particleA[1][i] * particleB[1][i] - 2 * particleA[2][i] * particleB[2][i];
    }

    return sValues;
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
    const char *   pxSuffix     = "_Px";
    const char *   pySuffix     = "_Py";
    const char *   pzSuffix     = "_Pz";
    const char *   energySuffix = "_E";

    // Generate names of the particle's momenta and energy
    char particlePx[NAME_LENGTH];
    strcpy(particlePx, particleName);
    strcat(particlePx, pxSuffix);

    char particlePy[NAME_LENGTH];
    strcpy(particlePy, particleName);
    strcat(particlePy, pySuffix);

    char particlePz[NAME_LENGTH];
    strcpy(particlePz, particleName);
    strcat(particlePz, pzSuffix);

    char particleE[NAME_LENGTH];
    strcpy(particleE, particleName);
    strcat(particleE, energySuffix);

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

    // Number of datapoints
    const unsigned int length = myTree->GetEntries();

    // Create an array of arrays pointing to the K data
    double **   kArrays{nullptr};
    const char *kParticleName = "_1_K~";
    kArrays                   = writeArrays(myTree, kParticleName);

    // Do the same for the pi- data
    double **   pi1Arrays{nullptr};
    const char *pi1ParticleName = "_2_pi#";
    pi1Arrays                   = writeArrays(myTree, pi1ParticleName);

    double **   pi2Arrays{nullptr};
    const char *pi2ParticleName = "_3_pi#";
    pi2Arrays                   = writeArrays(myTree, pi2ParticleName);

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

    double *s01Values{nullptr};
    double *s02Values{nullptr};
    s01Values = s(kArrays, pi1Arrays, length);
    s02Values = s(kArrays, pi2Arrays, length);

    auto tmpCanvas = new TCanvas("c", "c", 600, 600);

    // 1d histogram of CoM energy
    // @@@ works for s01; all s02 are the same as pi2-(pz) for some reason
    // @@@ sometimes it works as expected if you restart root WHAT THE heck
    // @@@ possibly memory issue
    TH1D *foo = new TH1D("foo", "bar", 100, 0.35, 2.6);
    foo->FillN(length, s01Values, 0);
    foo->Draw();

    //auto  tmpCanvas2 = new TCanvas("d", "d", 600, 600);
    //TH2D *baz        = new TH2D("baz", "zoop", 100, 0.35, 2.6, 100, 0.35, 2.6);
    //baz->FillN(length, s01Values, s02Values, 0);
    //baz->Draw();

    std::cout << "==== Global: =====================" << std::endl;
    std::cout << "==================================" << std::endl;

    delete s01Values;
    delete s02Values;

    for (int i = 0; i < 4; ++i) {
        delete kArrays[i];
        delete pi1Arrays[i];
        delete pi2Arrays[i];
    }
    delete kArrays;
    delete pi1Arrays;
    delete pi2Arrays;
}

