/*
 * bin_generated_decays.cpp
 * ROOT macro to bin AmpGen generated D -> K3pi decays into predefined bins
 *
 * NOTE: if this fails to build with error "dcs.so: No such file or directory", try restarting ROOT and building again.
 *
 * This script contains no error handling, arg verification or anything to make it work nicely
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
// Max length of a branch name in the ROOT file of interest
#define NAME_LENGTH 10

// DCS and CF relative amplitude and phase
#define DCS_MAGNITUDE 0.0445
#define DCS_PHASE -3.04

/// Bin limits in phase, centred on zero by construction
#define NUM_BINS 5
#define BIN_LIMITS -39, 0, 43, 180 // not sure what to set these to

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
 * Plot the K energies and make a plot of s01 vs s02 to check consistency with the ROOT TBrowser
 *
 */
void plot_things(double **kArrays, double **pi1Arrays, double **pi2Arrays, size_t length)
{

    // Plot K energies
    auto  kCanvas = new TCanvas("K energies", "K energies", 600, 600);
    TH1D *hist    = new TH1D("K energies", "K energies", 100, 0.45, 1);
    hist->FillN(length, kArrays[3], 0);
    hist->Draw();

    // Plot CoM energies on a new canvas
    auto    comCanvas = new TCanvas("CoM Energies", "CoM Energies", 600, 600);
    double *s01       = s(kArrays, pi1Arrays, length);
    double *s02       = s(kArrays, pi2Arrays, length);
    TGraph *myGraph   = new TGraph(length, s(kArrays, pi1Arrays, length), s(kArrays, pi2Arrays, length));
    myGraph->Draw("AP");

    delete[] s01;
    delete[] s02;
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
    // values
    myTree->ResetBranchAddresses();
}

/*
 * Write an array of arrays for storing particle data
 *
 * Allocates memory for the data arrays which must be freed by the caller.
 * Also allocates memory to the array of arrays which must be freed
 */
double **writeArrays(TTree *myTree, const char *name, size_t length)
{

    // Generate names for the arrays
    // This is a stupid way to initialise an array but i can't think of a better one
    char xArrayName[NAME_LENGTH]{
        name[0], name[1], name[2], name[3], name[4], name[5], name[6], name[7], name[8], name[9]};
    char yArrayName[NAME_LENGTH]{
        name[0], name[1], name[2], name[3], name[4], name[5], name[6], name[7], name[8], name[9]};
    char zArrayName[NAME_LENGTH]{
        name[0], name[1], name[2], name[3], name[4], name[5], name[6], name[7], name[8], name[9]};
    char eArrayName[NAME_LENGTH]{
        name[0], name[1], name[2], name[3], name[4], name[5], name[6], name[7], name[8], name[9]};

    strcat(xArrayName, "_Px");
    strcat(yArrayName, "_Py");
    strcat(zArrayName, "_Pz");
    strcat(eArrayName, "_E");

    // Allocate memory for data arrays
    double *const xArray = new double[length];
    double *const yArray = new double[length];
    double *const zArray = new double[length];
    double *const eArray = new double[length];

    // Create an array of arrays
    double **arrays = new double *[4] { xArray, yArray, zArray, eArray };

    // Write the data to the arrays
    writeArray(myTree, xArrayName, xArray);
    writeArray(myTree, yArrayName, yArray);
    writeArray(myTree, zArrayName, zArray);
    writeArray(myTree, eArrayName, eArray);

    return arrays;
}

/*
 * Bin the decays modelled an AmpGen generated inputFile into phase bins as defined by $BIN_LIMITS
 *
 */
void bin_generated_decays(TFile *inputFile)
{
    // Read in the tree and branches from the provided ROOT file
    // The tree of interest is DalitzEventList as this contains the decay products' kinematic data
    TTree *myTree = nullptr;
    inputFile->GetObject("DalitzEventList", myTree);

    // Number of datapoints
    const unsigned int length = myTree->GetEntries();

    // Create an arrays of arrays pointing to the particle data
    double **kArrays   = writeArrays(myTree, "_1_K~", length);
    double **pi1Arrays = writeArrays(myTree, "_2_pi#", length);
    double **pi2Arrays = writeArrays(myTree, "_3_pi#", length);
    double **pi3Arrays = writeArrays(myTree, "_4_pi~", length);

    // Create vectors of particle data
    std::vector<double> kPxVector = std::vector<double>(length);
    double              myData{0.0};

    // Point the desired branch at the myData variable.
    myTree->SetBranchAddress("_1_K~_Px", &myData);

    // Create an array of the appropriate size to store this data
    for (Long64_t i = 0; i < myTree->GetEntries(); ++i) {
        myTree->GetEntry(i);
        kPxVector[i] = myData;
    }
    for (size_t i = 0; i < length; ++i) {
        std::cout << kPxVector[i] << std::endl;
    }
    myTree->ResetBranchAddresses();

    // Apply scaling and rotation to the DCS amplitude such that we get dcs/cf amplitude ratio 'r' = 0.055
    // and the average relative strong-phase between the two amplitudes ~ 0.
    std::complex<double> dcs_offset = DCS_MAGNITUDE * exp(std::complex<double>(0, 1) * DCS_PHASE * M_PI / 180.);

    // Define the bins based on the form of the DCS and CF decays
    const std::string     dcsFile{"binning/dcs.so"};
    const std::string     cfFile{"binning/cf.so"};
    k3pi_binning::binning bins(dcsFile, cfFile, dcs_offset, {BIN_LIMITS});

    /// initialise global hadronic parameters, and parameters in each of the bins.
    std::complex<double>              z(0, 0);
    double                            n_cf(0);
    double                            n_dcs(0);
    std::vector<std::complex<double>> z_binned(NUM_BINS, std::complex<double>(0, 0));
    std::vector<double>               n_cf_binned(NUM_BINS, 0);
    std::vector<double>               n_dcs_binned(NUM_BINS, 0);

    for (int i = 0; i < myTree->GetEntries(); ++i) {
        // Create a vector of TLorentzVectors for this event (K+, pi-, pi-, pi+)
        TLorentzVector kLorentzVector{kArrays[0][i], kArrays[1][i], kArrays[2][i], kArrays[3][i]};
        TLorentzVector pi1LorentzVector{pi1Arrays[0][i], pi1Arrays[1][i], pi1Arrays[2][i], pi1Arrays[3][i]};
        TLorentzVector pi2LorentzVector{pi2Arrays[0][i], pi2Arrays[1][i], pi2Arrays[2][i], pi2Arrays[3][i]};
        TLorentzVector pi3LorentzVector{pi3Arrays[0][i], pi3Arrays[1][i], pi3Arrays[2][i], pi3Arrays[3][i]};

        std::vector<TLorentzVector> eventVector{kLorentzVector, pi1LorentzVector, pi2LorentzVector, pi3LorentzVector};
        auto                        event = k3pi_binning::eventFromVectors(eventVector);

        // Work out the CF and DCS amplitudes of this event, using the bins object created above
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

    // Make some plots to check that the data from ROOT has been read in correctly
    // plot_things(kArrays, pi1Arrays, pi2Arrays, length);

    // Output hadronic parameters
    std::cout << "==== Global: =====================" << std::endl;
    std::cout << "R = " << std::abs(z) / sqrt(n_cf * n_dcs) << "\t// Should be ~0.47" << std::endl;
    std::cout << "d = " << std::arg(z) * 180 / M_PI << "\t// Should be ~ zero by construction (< 1 degree)"
              << std::endl;
    std::cout << "r = " << sqrt(n_dcs / n_cf) << "\t// Should be ~ 0.055" << std::endl;

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
    delete[] kArrays;
    delete[] pi1Arrays;
    delete[] pi2Arrays;
    delete[] pi3Arrays;
}
