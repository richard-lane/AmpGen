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
// DCS and CF relative amplitude and phase
#define DCS_MAGNITUDE 0.0445
#define DCS_PHASE -3.04

/// Bin limits in phase, centred on zero by construction
#define NUM_BINS 5
#define BIN_LIMITS -39, 0, 43, 180 // not sure what to set these to

/*
 * From a vector of TLorentzVectors and the desired index (0,1,2,3), find a C-style array of data
 *
 * Allocates memory to the array which must be freed by the caller.
 *
 * e.g. vector2Array(myVector, 0) for x-momentum
 */
double *vector2Array(const std::vector<TLorentzVector> &particleVector, const size_t index)
{
    size_t  length   = particleVector.size();
    double *outArray = new double[length];

    for (size_t i = 0; i < length; ++i) {
        outArray[i] = particleVector[i][index];
    }

    return outArray;
}

/*
 * CoM energy squared of the ab system
 * Takes two vectors of TLorentzVectors as particle data
 *
 * Assumes each entry of the vector is of the form (Px, Py, Pz, E)
 *
 */
const std::vector<double> s(const std::vector<TLorentzVector> &particleA, const std::vector<TLorentzVector> &particleB)
{
    size_t              length = particleA.size();
    std::vector<double> sValues(particleA.size());

    for (size_t i = 0; i < length; ++i) {

        sValues[i] = std::pow(particleA[i][3], 2) - std::pow(particleA[i][0], 2) - std::pow(particleA[i][1], 2) -
                     std::pow(particleA[i][2], 2) + std::pow(particleB[i][3], 2) - std::pow(particleB[i][0], 2) -
                     std::pow(particleB[i][1], 2) - std::pow(particleB[i][2], 2) +
                     2 * particleA[i][3] * particleB[i][3] - 2 * particleA[i][0] * particleB[i][0] -
                     2 * particleA[i][1] * particleB[i][1] - 2 * particleA[i][2] * particleB[i][2];
    }

    return sValues;
}

/*
 * Plot a 100-bin histogram from an array
 */
void plot_hist(const std::string &title, double *myData, size_t length, const float xmin, const float xmax)
{

    const char *titleStr = title.c_str();
    auto        kCanvas  = new TCanvas(titleStr, titleStr, 600, 600);
    TH1D *      hist     = new TH1D(titleStr, titleStr, 100, xmin, xmax);
    hist->FillN(length, myData, 0);
    hist->Draw();
}

/*
 * Plot the K energies and make a plot of s01 vs s02 to check consistency with the ROOT TBrowser
 *
 */
void plot_things(const std::vector<TLorentzVector> &kVectors,
                 const std::vector<TLorentzVector> &pi1Vectors,
                 const std::vector<TLorentzVector> &pi2Vectors)
{
    size_t length = kVectors.size();

    // Plot K energies
    plot_hist("K energies", vector2Array(kVectors, 3), length, 0.45, 1);

    // Plot CoM energies on a new canvas
    auto                      comCanvas = new TCanvas("CoM Energies", "CoM Energies", 600, 600);
    const std::vector<double> s01       = s(kVectors, pi1Vectors);
    const std::vector<double> s02       = s(kVectors, pi2Vectors);
    TGraph *                  myGraph   = new TGraph(length, s01.data(), s02.data());
    myGraph->Draw("AP");
}

/*
 * Write the data on branchName to each TLorentzVector in myVector.
 */
inline void writeData(TTree &myTree, const std::string &branchName, std::vector<double> &myVector)
{
    double myData{0.0};

    myTree.SetBranchAddress(branchName.c_str(), &myData);
    for (Long64_t i = 0; i < myTree.GetEntries(); ++i) {
        myTree.GetEntry(i);
        myVector[i] = myData;
    }
    // Reset all branch addresses to avoid a bug where repeatedly calling this function would set an array to the wrong
    // values
    myTree.ResetBranchAddresses();
}

/*
 * Write the data on branchName to the index'th position of each TLorentzVector in myVector.
 * e.g. to write x-momenta of a particle described by ROOT branch foo_Px, call writeData("foo_Px", myVector, 0)
 *
 * The TLorentzVector should be of the form (Px, Py, Pz, E).
 */
inline void
writeData(TTree &myTree, const std::string &branchName, std::vector<TLorentzVector> &myVector, const size_t &index)
{
    double myData{0.0};

    myTree.SetBranchAddress(branchName.c_str(), &myData);
    for (Long64_t i = 0; i < myTree.GetEntries(); ++i) {
        myTree.GetEntry(i);
        myVector[i][index] = myData;
    }
    // Reset all branch addresses to avoid a bug where repeatedly calling this function would set an array to the wrong
    // values
    myTree.ResetBranchAddresses();
}

/*
 * Write a vector of TLorentzVectors containing the data for the particle branchName
 */
std::vector<TLorentzVector> writeVector(TTree &myTree, const std::string &particleName)
{
    std::vector<TLorentzVector> myVector = std::vector<TLorentzVector>(myTree.GetEntries());

    writeData(myTree, particleName + "_Px", myVector, 0);
    writeData(myTree, particleName + "_Py", myVector, 1);
    writeData(myTree, particleName + "_Pz", myVector, 2);
    writeData(myTree, particleName + "_E", myVector, 3);

    return myVector;
}

/*
 * Bin the decays modelled in an AmpGen generated inputFile into phase bins as defined by $BIN_LIMITS
 *
 */
void bin_generated_decays(TFile *inputFile)
{
    // ---- Parameters
    // Scaling and rotation to the DCS amplitude such that we get dcs/cf amplitude ratio 'r' = 0.055
    // and the average relative strong-phase between the two amplitudes ~ 0.
    const std::complex<double> dcs_offset = DCS_MAGNITUDE * exp(std::complex<double>(0, 1) * DCS_PHASE * M_PI / 180.);

    // Define the bins based on the form of the DCS and CF decays
    const std::string     dcsFile{"binning/dcs.so"};
    const std::string     cfFile{"binning/cf.so"};
    k3pi_binning::binning bins(dcsFile, cfFile, dcs_offset, {BIN_LIMITS});

    // ---- Calculations
    // Read in the tree and branches from the provided ROOT file
    // The tree of interest is DalitzEventList as this contains the decay products' kinematic data
    TTree *myTree = nullptr;
    inputFile->GetObject("DalitzEventList", myTree);
    long long length{myTree->GetEntries()};

    // Read in vectors of particle data from the ROOT file
    const std::vector<TLorentzVector> kVectors   = writeVector(*myTree, "_1_K~");
    const std::vector<TLorentzVector> pi1Vectors = writeVector(*myTree, "_2_pi#");
    const std::vector<TLorentzVector> pi2Vectors = writeVector(*myTree, "_3_pi#");
    const std::vector<TLorentzVector> pi3Vectors = writeVector(*myTree, "_4_pi~");

    // Create a vector for each bin holding pairs of event DCS/CF amplitudes and the time
    // Create vectors for each bin holding DCS/CF amplitudes and the time
    std::vector<std::vector<double>> binRatio(NUM_BINS, std::vector<double>(length, 0.0));
    std::vector<std::vector<double>> binTimes(NUM_BINS, std::vector<double>(length, -1.0));

    std::vector<double> times(length, -1);
    writeData(*myTree, "D_decayTime", times);

    for (int i = 0; i < length; ++i) {
        // Create a vector of TLorentzVectors for this event (K+, pi-, pi-, pi+)
        std::vector<TLorentzVector> eventVector{kVectors[i], pi1Vectors[i], pi2Vectors[i], pi3Vectors[i]};
        auto                        event = k3pi_binning::eventFromVectors(eventVector);

        // Work out the CF and DCS amplitudes of this event, using the bins object created above
        auto eval_cf  = bins.cf(event.data(), 1);
        auto eval_dcs = dcs_offset * bins.dcs(event.data(), 1);

        // Find which bin the event belongs in
        int bin = bins.bin(eventVector, 1);

        // Log the time and ratio of the event in this bin
        double dcsCfRatio = abs(eval_dcs / eval_cf);
        binRatio[bin].push_back(dcsCfRatio);
        binTimes[bin].push_back(times[i]);
    }

    // Make some plots to check that the data from ROOT has been read in correctly
    // plot_things(kVectors, pi1Vectors, pi2Vectors);
}
