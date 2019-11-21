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
#include <numeric>
#include <string>
#include <vector>

#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"

#include "k3pi_binning.h"
#include "plottingHelpers.cpp"

// ---- Magic Numbers
// DCS and CF relative amplitude and phase
#define DCS_MAGNITUDE 0.0445
#define DCS_PHASE -3.04

/// Bin limits in phase, centred on zero by construction
#define NUM_BINS 5
#define BIN_LIMITS -39, 0, 43, 180 // not sure what to set these to

/*
 * Write the data on branchName to each TLorentzVector in myVector.
 */
void writeData(TTree &myTree, const std::string &branchName, std::vector<double> &myVector)
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
void
writeDataToLorentzVectors(TTree &myTree, const std::string &branchName, std::vector<TLorentzVector> &myVector, const size_t &index)
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

    writeDataToLorentzVectors(myTree, particleName + "_Px", myVector, 0);
    writeDataToLorentzVectors(myTree, particleName + "_Py", myVector, 1);
    writeDataToLorentzVectors(myTree, particleName + "_Pz", myVector, 2);
    writeDataToLorentzVectors(myTree, particleName + "_E", myVector, 3);

    return myVector;
}

/*
 * Split a vector of N vectors into a a vector containing N vectors of vectors, each of which has a maximum of chunkSize
 * elements
 */
std::vector<std::vector<std::vector<double>>> splitVectors(const std::vector<std::vector<double>> &myVector,
                                                           size_t                                  chunkSize)
{
    size_t                                        N{myVector.size()};
    std::vector<std::vector<std::vector<double>>> outVectors(N);

    for (size_t i = 0; i < N; ++i) {
        // Find the number of subvectors to split this vector into
        size_t vectorSize = myVector[i].size();
        size_t numChunks{vectorSize / chunkSize + (vectorSize % chunkSize != 0)};

        // Append numChunks empty vectors to outVectors[i]
        for (size_t chunk = 0; chunk < numChunks; ++chunk) {
            std::vector<double> emptySubvector;
            outVectors[i].push_back(
                emptySubvector); // this looks wrong, shouldn't outVectors[i] be a vector of vectors?
        }

        // Loop over all the data in the vector by index
        // Find which chunk the data belongs in, push it back
        for (size_t j = 0; j < vectorSize; ++j) {
            size_t chunkNumber = j / chunkSize;
            outVectors[i][chunkNumber].push_back(myVector[i][j]);
        }
    }

    return outVectors;
}

/*
 * Find average of vector
 */
double vectorAvg(const std::vector<double> &vector)
{
    return std::accumulate(std::begin(vector), std::end(vector), 0.0) / vector.size();
}

/*
 * Find std dev of a vector
 * Could improve implementation so that the avg is passed into this fcn but i dont want to do that
 */
double vectorStdDev(const std::vector<double> &vector)
{
    size_t              size = vector.size();
    double              mean = vectorAvg(vector);
    std::vector<double> diff(size);

    std::transform(vector.begin(), vector.end(), diff.begin(), std::bind2nd(std::minus<double>(), mean));
    double sqSum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);

    return std::sqrt(sqSum / size);
}

/*
 * Given a vector of vector of pairs, sort each vector of pairs into ascending order based on the second element in the
 * pair.
 *
 */
void sortVectorOfPairs(std::vector<std::vector<std::pair<double, double>>> &myVector)
{
    // Sort each vector based on time
    for (size_t i = 0; i < myVector.size(); ++i) {
        std::sort(
            myVector[i].begin(), myVector[i].end(), [](auto &left, auto &right) { return left.second < right.second; });
    }
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

    // Read the time data from the ROOT file into a vector
    std::vector<double> times(length, -1);
    writeData(*myTree, "D_decayTime", times);

    // Create a vector of pairs of event ratio/time for each bin
    // Store these in a vector for convenience
    std::vector<std::vector<std::pair<double, double>>> binData(NUM_BINS, std::vector<std::pair<double, double>>());

    // Make some plots to check that the data from ROOT has been read in correctly
    plot_things(kVectors, pi1Vectors, pi2Vectors);

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
        // This might be slow because it's dynamically resizing the vector, but should be ok for our purposes.
        double dcsCfRatio = abs(eval_dcs / eval_cf);
        binData[bin].push_back(std::make_pair(dcsCfRatio, times[i]));
    }

    // Find the number of points in each bin
    std::vector<size_t> binSizes(NUM_BINS);
    for (size_t bin = 0; bin < NUM_BINS; ++bin) {
        binSizes[bin] = binData[bin].size();
        std::cout << "points in bin " << bin << ": " << binSizes[bin] << std::endl;
    }

    // Sort the data in each bin in increasing time order
    sortVectorOfPairs(binData);

    // Create two vectors for each bin holding ratio and time data
    std::vector<std::vector<double>> binRatios(NUM_BINS);
    std::vector<std::vector<double>> binTimes(NUM_BINS);
    for (size_t i = 0; i < NUM_BINS; ++i) {
        for (size_t j = 0; j < binData[i].size(); ++j) {
            binRatios[i].push_back(binData[i][j].first);
            binTimes[i].push_back(binData[i][j].second);
        }
    }

    // Split vectors of ratios and times into subvectors
    size_t                                        binSize{100};
    std::vector<std::vector<std::vector<double>>> splitBinRatios = splitVectors(binRatios, binSize);
    std::vector<std::vector<std::vector<double>>> splitBinTimes  = splitVectors(binTimes, binSize);

    // Find average and std dev of ratios and times in each subvector
    std::vector<std::vector<double>> binRatioAverage(NUM_BINS);
    std::vector<std::vector<double>> binTimesAverage(NUM_BINS);
    std::vector<std::vector<double>> binRatioStdDev(NUM_BINS);
    std::vector<std::vector<double>> binTimesStdDev(NUM_BINS);

    // Create vectors of the right length to hold the average and std devs for each bin
    for (size_t bin = 0; bin < splitBinTimes.size(); ++bin) {
        binRatioAverage[bin] = std::vector<double>(splitBinRatios[bin].size());
        binTimesAverage[bin] = std::vector<double>(splitBinRatios[bin].size());
        binRatioStdDev[bin]  = std::vector<double>(splitBinRatios[bin].size());
        binTimesStdDev[bin]  = std::vector<double>(splitBinRatios[bin].size());
    }

    for (size_t bin = 0; bin < NUM_BINS; ++bin) {
        for (size_t i = 0; i < splitBinTimes[bin].size(); ++i) {
            binRatioAverage[bin][i] = vectorAvg(splitBinRatios[bin][i]);
            binTimesAverage[bin][i] = vectorAvg(splitBinTimes[bin][i]);

            binRatioStdDev[bin][i] = vectorStdDev(splitBinRatios[bin][i]);
            binTimesStdDev[bin][i] = vectorStdDev(splitBinTimes[bin][i]);
        }
    }
}
