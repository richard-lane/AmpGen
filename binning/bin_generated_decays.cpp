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
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TRandom.h"
#include "TTree.h"

#include "k3pi_binning.h"

// ---- Magic Numbers
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
double *s(double **particleA, double **particleB, const unsigned int length)
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
 * Plot CoM energies s01 vs s02
 *
 */
void plotS01VsS02(const char *title, double **kArrays, double **pi1Arrays, double **pi2Arrays, size_t length)
{
    const double *s01Values{nullptr};
    const double *s02Values{nullptr};

    s01Values = s(kArrays, pi1Arrays, length);
    s02Values = s(kArrays, pi2Arrays, length);

    auto    myCanvas = new TCanvas("CoM Energies", "CoM Energies", 600, 600);
    TGraph *myGraph  = new TGraph(length, s01Values, s02Values);
    myGraph->Draw("AP");

    delete s01Values;
    delete s02Values;
}

/*
 * Plot an array
 */
void plotArray(const char *  title,
               const double *array,
               size_t        length,
               size_t        nBins = 100,
               double        xMin  = 0.0,
               double        xMax  = 1.0)
{
    auto  myCanvas = new TCanvas(title, title, 600, 600);
    TH1D *hist     = new TH1D(title, title, nBins, xMin, xMax);

    hist->FillN(length, array, 0);
    hist->Draw();
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

    // Read in the tree and branches from the provided ROOT file
    // The tree of interest is DalitzEventList as this contains the decay products' kinematic data
    TTree *myTree = nullptr;
    inputFile->GetObject("DalitzEventList", myTree);

    // Number of datapoints
    const unsigned int length = myTree->GetEntries();

    // Create an array of arrays pointing to the K data
    // Lots of reused code here but the effort to refactor won't be worth the trouble
    double *const kPxArray   = new double[length];
    double *const kPyArray   = new double[length];
    double *const kPzArray   = new double[length];
    double *const kEArray    = new double[length];
    double *kArrays[4] = {kPxArray, kPyArray, kPzArray, kEArray};

    writeArray(myTree, "_1_K~_E", kEArray);
    writeArray(myTree, "_1_K~_Px", kPxArray);
    writeArray(myTree, "_1_K~_Py", kPyArray);
    writeArray(myTree, "_1_K~_Pz", kPzArray);

    // Create arrays pointing to the pi- data
    double *const pi1PxArray   = new double[length];
    double *const pi1PyArray   = new double[length];
    double *const pi1PzArray   = new double[length];
    double *const pi1EArray    = new double[length];
    double *pi1Arrays[4] = {pi1PxArray, pi1PyArray, pi1PzArray, pi1EArray};

    double *const pi2PxArray   = new double[length];
    double *const pi2PyArray   = new double[length];
    double *const pi2PzArray   = new double[length];
    double *const pi2EArray    = new double[length];
    double *pi2Arrays[4] = {pi2PxArray, pi2PyArray, pi2PzArray, pi2EArray};

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
    double *pi3Arrays[4] = {pi3PxArray, pi3PyArray, pi3PzArray, pi3EArray};

    writeArray(myTree, "_4_pi~_E", pi3EArray);
    writeArray(myTree, "_4_pi~_Px", pi3PxArray);
    writeArray(myTree, "_4_pi~_Py", pi3PyArray);
    writeArray(myTree, "_4_pi~_Pz", pi3PzArray);


    /// initialise global hadronic parameters, and parameters in each of the bins.
    std::complex<double>              z(0, 0);
    double                            n_cf(0);
    double                            n_dcs(0);
    std::vector<std::complex<double>> z_binned(NUM_BINS, std::complex<double>(0, 0));
    std::vector<double>               n_cf_binned(NUM_BINS, 0);
    std::vector<double>               n_dcs_binned(NUM_BINS, 0);

    for (int i = 0; i < myTree->GetEntries(); ++i) {
        TLorentzVector kLorentzVector{kArrays[0][i], kArrays[1][i], kArrays[2][i], kArrays[3][i]};
        TLorentzVector pi1LorentzVector{pi1Arrays[0][i], pi1Arrays[1][i], pi1Arrays[2][i], pi1Arrays[3][i]};
        TLorentzVector pi2LorentzVector{pi2Arrays[0][i], pi2Arrays[1][i], pi2Arrays[2][i], pi2Arrays[3][i]};
        TLorentzVector pi3LorentzVector{pi3Arrays[0][i], pi3Arrays[1][i], pi3Arrays[2][i], pi3Arrays[3][i]};

        // @@@ todo:
        //      Create a vector of TLorentzVectors for this event (K+, pi-, pi-, pi+)
        std::vector<TLorentzVector> eventVector{kLorentzVector, pi1LorentzVector, pi2LorentzVector, pi3LorentzVector};
        auto event = k3pi_binning::eventFromVectors(eventVector);

        //      @@@ Work out the CF and DCS amplitudes of this event, using the bins object created above
        auto eval_cf  = bins.cf(event.data(), 1);
        auto eval_dcs = dcs_offset * bins.dcs(event.data(), 1);

        // Update the global hadronic parameters
        z += eval_cf * std::conj(eval_dcs);
        n_cf += std::norm(eval_cf);
        n_dcs += std::norm(eval_dcs);

        // Find which bin the event belongs in and update its hadronic parameters
        auto bin = bins.bin(eventVector, 1);
        z_binned[bin] += eval_cf * std::conj(eval_dcs);
        n_cf_binned[bin] += std::norm(eval_cf);
        n_dcs_binned[bin] += std::norm(eval_dcs);

    }

    //double *s01Values{nullptr};
    //double *s02Values{nullptr};
    //s01Values = s(kArrays, pi1Arrays, length);
    //s02Values = s(kArrays, pi2Arrays, length);

    //plotArray("K energies", kArrays[3], length);
    //plotS01VsS02("S01 vs S02", kArrays, pi1Arrays, pi2Arrays, length);

    // Output hadronic parameters
    std::cout << "==== Global: =====================" << std::endl;
    std::cout << "R = " << std::abs(z) / sqrt(n_cf * n_dcs) << std::endl; /// should be ~ 0.47
    std::cout << "d = " << std::arg(z) * 180 / M_PI << std::endl; /// should be ~ zero by construction (< 1 degree)
    std::cout << "r = " << sqrt(n_dcs / n_cf) << std::endl;       /// should be ~ 0.055

    for (int i = 0; i < NUM_BINS; ++i) {
        std::cout << "==== Bin " << i + 1 << ": ======================" << std::endl;
        std::cout << "R[" << i + 1 << "] = " << std::abs(z_binned[i]) / sqrt(n_cf_binned[i] * n_dcs_binned[i])
                  << std::endl;
        std::cout << "d[" << i + 1 << "] = " << std::arg(z_binned[i]) * 180 / M_PI << std::endl;
        std::cout << "K[" << i + 1 << "] = " << n_dcs_binned[i] / n_dcs << std::endl;
        std::cout << "K'[" << i + 1 << "] = " << n_cf_binned[i] / n_cf << std::endl;
    }
    std::cout << "==================================" << std::endl;


    // free arrays
    for (int i = 0; i < 4; ++i) {
        delete[] kArrays[i];
        delete[] pi1Arrays[i];
        delete[] pi2Arrays[i];
        delete[] pi3Arrays[i];
    }
}
