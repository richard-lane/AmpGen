/*
 * bin_generated_decays.cpp
 * ROOT macro to bin AmpGen generated D -> K3pi decays into predefined bins
 *
 * NOTE: if this fails to build with error "dcs.so: No such file or directory", try restarting ROOT and building again.
 *
 * This script is gross and needs a lot of clearing up; only exists to test my understanding
 */
#include <cmath>
#include <complex>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TRandom.h"
#include "TTree.h"

#include "k3pi_binning.h"

// ---- Magic Numbers
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
 * Takes two arrays of arrays of kinematic particle data
 *
 * Assumes particle data arrays take the form *(Px, Py, Pz, E)
 *
 * Allocates memory to the array of s values which must be freed by the caller
 *
 */
double *s(const double **particleA, const double **particleB, const unsigned int length)
{
    // Initialise an array of CoM energies
    double *const sValues = new double[length];

    for (size_t i = 0; i < length; ++i) {
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
 */
void writeArray(TTree *myTree, const char *myBranchName, double *myArray)
{
    const long long numEntries{myTree->GetEntries()};
    double          myData{0.0};

    // Point the desired branch at the myData variable.
    myTree->SetBranchAddress(myBranchName, &myData);

    // Create an array of the appropriate size to store this data
    for (Long64_t i = 0; i < myTree->GetEntries(); ++i) {
        myTree->GetEntry(i);
        myArray[i] = myData;
    }

    // Reset all branch addresses to avoid a bug where repeatedly calling this function would set an array to the wrong
    // value
    myTree->ResetBranchAddresses();
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

    // Create bins based on DCS and CF amplitudes
    k3pi_binning::binning bins(dcsFile, cfFile, dcs_offset, {BIN_LIMITS});

    /// initialise global hadronic parameters, and parameters in each of the bins.
    std::complex<double>              z(0, 0);
    double                            n_cf(0);
    double                            n_dcs(0);
    std::vector<std::complex<double>> z_binned(NUM_BINS, std::complex<double>(0, 0));
    std::vector<double>               n_cf_binned(NUM_BINS, 0);
    std::vector<double>               n_dcs_binned(NUM_BINS, 0);

    // Read in the tree and branches from the provided ROOT file
    // The tree of interest is DalitzEventList as this contains the decay products' kinematic data
    TTree *myTree = nullptr;
    inputFile->GetObject("DalitzEventList", myTree);

    // Number of datapoints
    const unsigned int length = myTree->GetEntries();

    // Create an array of arrays pointing to the K data
    double *const kPxArray   = new double[length];
    double *const kPyArray   = new double[length];
    double *const kPzArray   = new double[length];
    double *const kEArray    = new double[length];
    const double *kArrays[4] = {kPxArray, kPyArray, kPzArray, kEArray};

    writeArray(myTree, "_1_K~_E", kEArray);
    writeArray(myTree, "_1_K~_Px", kPxArray);
    writeArray(myTree, "_1_K~_Py", kPyArray);
    writeArray(myTree, "_1_K~_Pz", kPzArray);

    // Create arrays pointing to the pi- data
    double *const pi1PxArray   = new double[length];
    double *const pi1PyArray   = new double[length];
    double *const pi1PzArray   = new double[length];
    double *const pi1EArray    = new double[length];
    const double *pi1Arrays[4] = {pi1PxArray, pi1PyArray, pi1PzArray, pi1EArray};

    double *const pi2PxArray   = new double[length];
    double *const pi2PyArray   = new double[length];
    double *const pi2PzArray   = new double[length];
    double *const pi2EArray    = new double[length];
    const double *pi2Arrays[4] = {pi2PxArray, pi2PyArray, pi2PzArray, pi2EArray};

    writeArray(myTree, "_2_pi#_E", pi1EArray);
    writeArray(myTree, "_2_pi#_Px", pi1PxArray);
    writeArray(myTree, "_2_pi#_Py", pi1PyArray);
    writeArray(myTree, "_2_pi#_Pz", pi1PzArray);

    writeArray(myTree, "_3_pi#_E", pi2EArray);
    writeArray(myTree, "_3_pi#_Px", pi2PxArray);
    writeArray(myTree, "_3_pi#_Py", pi2PyArray);
    writeArray(myTree, "_3_pi#_Pz", pi2PzArray);

    // Create arrays for pi+ data
    double *const pi3PxArray   = new double[length];
    double *const pi3PyArray   = new double[length];
    double *const pi3PzArray   = new double[length];
    double *const pi3EArray    = new double[length];
    const double *pi3Arrays[4] = {pi3PxArray, pi3PyArray, pi3PzArray, pi3EArray};

    writeArray(myTree, "_4_pi~_E", pi3EArray);
    writeArray(myTree, "_4_pi~_Px", pi3PxArray);
    writeArray(myTree, "_4_pi~_Py", pi3PyArray);
    writeArray(myTree, "_4_pi~_Pz", pi3PzArray);

    //@@@ allocate memory to the 4 k arrays
    //@@@ refactor writeArray to take in an allocated pointer and the branch name
    //@@@ fill them with the data using the writeArray function

    // Do the same for the pi- data
    // double *   pi1Arrays[4];
    // const char *pi1ParticleName = "_2_pi#";
    // pi1Arrays                   = writeArrays(myTree, pi1ParticleName);

    // double *   pi2Arrays[4];
    // const char *pi2ParticleName = "_3_pi#";
    // pi2Arrays                   = writeArrays(myTree, pi2ParticleName);

    for (int i = 0; i < myTree->GetEntries(); ++i) {
        TLorentzVector kLorentzVector{};
        TLorentzVector pi1LorentzVector{};
        TLorentzVector pi2LorentzVector{};
        TLorentzVector pi3LorentzVector{}; //@@@ initialise these with values from the above allocated arrays

        // @@@ todo:
        //      Create a vector of TLorentzVectors for this event (K+, pi-, pi-, pi+)
        //      run eventFromVectors on it to get an event vector (or could just create the event vector straight away)

        //      @@@ Work out the CF and DCS amplitudes of this event, using the bins object created above
        //      @@@ Update global hadronic parameters
        //      @@@ Find which bin the event belongs in then update its hadronic parameters
    }

    //@@@bugfix- some sort of memory error when writing to these arrays; sometimes s is calculated incorrectly
    double *s01Values{nullptr};
    double *s02Values{nullptr};
    s01Values = s(kArrays, pi1Arrays, length);
    s02Values = s(kArrays, pi2Arrays, length);

    // auto tmpCanvas = new TCanvas("c", "c", 600, 600);

    // TH1D *foo = new TH1D("foo", "bar", 100, 0, 1);
    // foo->FillN(length, kArrays[2], 0);
    // foo->Draw();

    // Attempt to reproduce the s01 vs s02 plot to check consistency
    auto    tmpCanvas2 = new TCanvas("d", "d", 600, 600);
    TGraph *baz        = new TGraph(length, s01Values, s02Values);
    baz->Draw("AP");

    // Placeholder for outputting hadronic parameters
    std::cout << "==== Global: =====================" << std::endl;
    std::cout << "==================================" << std::endl;

    // free arrays
    delete[] kPxArray;
    delete[] kPyArray;
    delete[] kPzArray;
    delete[] kEArray;

    delete[] pi1PxArray;
    delete[] pi1PyArray;
    delete[] pi1PzArray;
    delete[] pi1EArray;

    delete[] pi2PxArray;
    delete[] pi2PyArray;
    delete[] pi2PzArray;
    delete[] pi2EArray;

    delete[] pi3PxArray;
    delete[] pi3PyArray;
    delete[] pi3PzArray;
    delete[] pi3EArray;
}
