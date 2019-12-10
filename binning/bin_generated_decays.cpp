/*
 * bin_generated_decays.cpp
 * ROOT macro to bin AmpGen generated D -> K3pi decays into predefined bins
 *
 * NOTE: if this fails to build with error "dcs.so: No such file or directory", try restarting ROOT and building again.
 *
 * This script contains no error handling, arg verification or anything to make it work nicely
 */

#include <algorithm>
#include <assert.h>
#include <cmath>
#include <complex>
#include <cstring>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "TFile.h"
#include "TGraphErrors.h"
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
#define BIN_LIMITS -39, 0, 43, 180

/*
 * Write the data on branchName to myVector.
 */
void saveBranchToVector(TTree &myTree, const std::string &branchName, std::vector<double> &myVector)
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
 * e.g. to write x-momenta of a particle described by ROOT branch foo_Px, call saveBranchToVector("foo_Px", myVector, 0)
 *
 * The TLorentzVector should be of the form (Px, Py, Pz, E).
 */
void writeBranchToLorentzVectors(TTree &                      myTree,
                                 const std::string &          branchName,
                                 std::vector<TLorentzVector> &myVector,
                                 const size_t &               index)
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

    writeBranchToLorentzVectors(myTree, particleName + "_Px", myVector, 0);
    writeBranchToLorentzVectors(myTree, particleName + "_Py", myVector, 1);
    writeBranchToLorentzVectors(myTree, particleName + "_Pz", myVector, 2);
    writeBranchToLorentzVectors(myTree, particleName + "_E", myVector, 3);

    return myVector;
}

/*
 * Split a vector of N vectors into a a vector containing N vectors of vectors, each of which has a maximum of chunkSize
 * elements
 */
std::vector<std::vector<std::vector<double>>> splitVectorsEqualChunks(const std::vector<std::vector<double>> &myVector,
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
 * Given a vector of vectors of pairs, split it into a vector of vectors of vectors of pairs.
 * Each vector is split based on the second element in the pair according to binLimits.
 *
 * Each subvector should be sorted in order of increasing second pair element.
 */
std::vector<std::vector<std::vector<std::pair<double, double>>>>
splitVectorsWithLimits(const std::vector<std::vector<std::pair<double, double>>> &myVectors,
                       std::vector<double>                                        binLimits)
{
    // Bin limits should be sorted
    if (!std::is_sorted(binLimits.begin(), binLimits.end())) {
        std::cout << "Bad time bin limits; should be sorted" << std::endl;
        throw;
    }

    // Bins for values up to binLimits[0], values between two limits and values >binLimits[N]
    size_t numBins = binLimits.size() + 1;

    // Create our vector of vectors of vectors of pairs
    std::vector<std::vector<std::vector<std::pair<double, double>>>> splitVectors(
        NUM_BINS, std::vector<std::vector<std::pair<double, double>>>(numBins));

    // In each bin, fill in our vector of vectors according to the bin limits
    for (size_t bin = 0; bin < NUM_BINS; ++bin) {
        size_t currentTimeBin = 0;
        for (size_t i = 0; i < myVectors[bin].size(); ++i) {
            // If this time is more than the bin limit, we want to put subsequent points in a higher bin
            // Unless we are already putting points in the highest bin
            while (myVectors[bin][i].second > binLimits[currentTimeBin] && currentTimeBin < numBins - 1) {
                currentTimeBin += 1;
            }
            splitVectors[bin][currentTimeBin].push_back(myVectors[bin][i]);
        }
    }
    return splitVectors;
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
void sortVectorsOfPairs(std::vector<std::vector<std::pair<double, double>>> &myVector)
{
    // Sort each vector based on time
    for (size_t i = 0; i < myVector.size(); ++i) {
        std::sort(
            myVector[i].begin(), myVector[i].end(), [](auto &left, auto &right) { return left.second < right.second; });
    }
}

/*
 * In each phase-space bin, bin the data by time into bins containing an equal number (binSize) of points.
 *
 * binRatioAverage etc. are out args that are modified by this function
 *
 * binRatios and binTimes should be sorted in order of increasing time
 *
 */
void binDataEqualSizeBins(std::vector<std::vector<double>> &binRatioAverage,
                          std::vector<std::vector<double>> &binTimesAverage,
                          std::vector<std::vector<double>> &binRatioStdDev,
                          std::vector<std::vector<double>> &binTimesStdDev,
                          std::vector<std::vector<double>> &binRatios,
                          std::vector<std::vector<double>> &binTimes,
                          size_t                            binSize)
{
    // Split vectors of ratios and times into subvectors
    std::vector<std::vector<std::vector<double>>> splitBinRatios = splitVectorsEqualChunks(binRatios, binSize);
    std::vector<std::vector<std::vector<double>>> splitBinTimes  = splitVectorsEqualChunks(binTimes, binSize);

    // Create vectors of the right length to hold the average and std devs for each bin
    for (size_t bin = 0; bin < splitBinTimes.size(); ++bin) {
        // TODO stop calling .size() so much
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

/*
 * Convert a vector of pairs to two vectors
 *
 */
void vectorOfPairs2vectors(std::vector<double> &                   firstOutVector,
                           std::vector<double> &                   secondOutVector,
                           std::vector<std::pair<double, double>> &vectorOfPairs)
{
    // Check that our two out args are empty vectors
    if (!firstOutVector.empty() || !secondOutVector.empty()) {
        std::cout << "Bad out args passed to vectorOfPairs2vectors" << std::endl;
        throw;
    }

    size_t numPairs = vectorOfPairs.size();
    for (size_t i = 0; i < numPairs; ++i) {
        firstOutVector.push_back(vectorOfPairs[i].first);
        secondOutVector.push_back(vectorOfPairs[i].second);
    }
}

/*
 * In each phase-space bin, bin the data by time into bins defined by timeBinLimits
 *
 * binRatioAverage etc. are out args that are modified by this function
 *
 * binRatios and binTimes should be sorted in order of increasing time
 *
 */
void binDataTimeBinLimits(std::vector<std::vector<double>> &                   binRatioAverage,
                          std::vector<std::vector<double>> &                   binTimesAverage,
                          std::vector<std::vector<double>> &                   binRatioStdDev,
                          std::vector<std::vector<double>> &                   binTimesStdDev,
                          std::vector<std::vector<std::pair<double, double>>> &binData,
                          std::vector<double> &                                timeBinLimits)
{
    size_t numTimeBins = timeBinLimits.size() + 1;

    // Create vectors of the right length to hold the average and std devs for each bin
    for (size_t bin = 0; bin < NUM_BINS; ++bin) {
        binRatioAverage[bin] = std::vector<double>(numTimeBins, -1);
        binTimesAverage[bin] = std::vector<double>(numTimeBins, -1);
        binRatioStdDev[bin]  = std::vector<double>(numTimeBins, -1);
        binTimesStdDev[bin]  = std::vector<double>(numTimeBins, -1);
    }

    // Split our vectors of pairs into vectors of vectors of pairs based on the bin limits
    std::vector<std::vector<std::vector<std::pair<double, double>>>> splitBinData =
        splitVectorsWithLimits(binData, timeBinLimits);

    for (size_t bin = 0; bin < NUM_BINS; ++bin) {
        for (size_t i = 0; i < numTimeBins; ++i) {
            // Create vectors for time and ratio in this phase space and time bin
            std::vector<double> ratios;
            std::vector<double> times;
            vectorOfPairs2vectors(ratios, times, splitBinData[bin][i]);

            binRatioAverage[bin][i] = vectorAvg(ratios);
            binTimesAverage[bin][i] = vectorAvg(times);

            binRatioStdDev[bin][i] = vectorStdDev(ratios);
            binTimesStdDev[bin][i] = vectorStdDev(times);
        }
    }
}

/*
 * Bin the decays modelled in an AmpGen generated inputFile into phase bins as defined by $BIN_LIMITS
 *
 */
void bin_generated_decays(TFile *inputFile)
{
    // Scaling and rotation of the DCS amplitude such that we get a small dcs/cf amplitude ratio strong phase.
    //
    // Ideally these should be 0.055 and 0 respectively, though at the moment there's no reason to believe that the
    // current values of DCS_MAGNITUDE and DCS_PHASE will do this for any AmpGen generated dataset.
    // NB: I don't actually think this value affects the DCS/CF ratio vs time plots
    const std::complex<double> dcs_offset = DCS_MAGNITUDE * exp(std::complex<double>(0, 1) * DCS_PHASE * M_PI / 180.);

    // Define the bins based on the form of the DCS and CF decays
    const std::string     dcsFile{"binning/dcs.so"};
    const std::string     cfFile{"binning/cf.so"};
    k3pi_binning::binning bins(dcsFile, cfFile, dcs_offset, {BIN_LIMITS});

    // Read in the tree and branches from the provided ROOT file
    // The tree of interest is DalitzEventList as this contains the decay products' kinematic data
    TTree *myTree = nullptr;
    inputFile->GetObject("DalitzEventList", myTree);
    long long numGeneratedEvents{myTree->GetEntries()};

    // Read in vectors of particle data from the ROOT file
    const std::vector<TLorentzVector> kVectors   = writeVector(*myTree, "_1_K~");
    const std::vector<TLorentzVector> pi1Vectors = writeVector(*myTree, "_2_pi#");
    const std::vector<TLorentzVector> pi2Vectors = writeVector(*myTree, "_3_pi#");
    const std::vector<TLorentzVector> pi3Vectors = writeVector(*myTree, "_4_pi~");

    // Read the time data from the ROOT file into a vector
    std::vector<double> times(numGeneratedEvents, -1);
    saveBranchToVector(*myTree, "D_decayTime", times);

    // Create a vector of pairs of event ratio;time for each bin
    // Store these in a vector for convenience
    std::vector<std::vector<std::pair<double, double>>> binData(NUM_BINS, std::vector<std::pair<double, double>>());

    // Make some plots to check that the data from ROOT has been read in correctly
    // plot_things(kVectors, pi1Vectors, pi2Vectors);

    for (int i = 0; i < numGeneratedEvents; ++i) {
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

    // Find the number of points in each bin and output to console
    std::vector<size_t> binSizes(NUM_BINS);
    for (size_t bin = 0; bin < NUM_BINS; ++bin) {
        binSizes[bin] = binData[bin].size();
        std::cout << "points in bin " << bin << ": " << binSizes[bin] << std::endl;
    }

    // Sort the data in each bin in increasing time order
    sortVectorsOfPairs(binData);

    // For each bin, create vectors holding ratio and time data
    std::vector<std::vector<double>> binRatios(NUM_BINS);
    std::vector<std::vector<double>> binTimes(NUM_BINS);
    for (size_t i = 0; i < NUM_BINS; ++i) {
        for (size_t j = 0; j < binSizes[i]; ++j) {
            // Dynamic resizing again... we could avoid this using the above number of points per bin but nah
            binRatios[i].push_back(binData[i][j].first);
            binTimes[i].push_back(binData[i][j].second);
        }
    }

    // Find average and std dev of ratios and times in each subvector
    std::vector<std::vector<double>> binRatioAverage(NUM_BINS);
    std::vector<std::vector<double>> binTimesAverage(NUM_BINS);
    std::vector<std::vector<double>> binRatioStdDev(NUM_BINS);
    std::vector<std::vector<double>> binTimesStdDev(NUM_BINS);

    // Bin the data into time bins of equal sizes
    // binDataEqualSizeBins(binRatioAverage, binTimesAverage, binRatioStdDev, binTimesStdDev, binRatios, binTimes, 10);

    // Bin data into time bins defined by a vector
    std::vector<double> timeBinLimits{};
    for (double i = 1; i < 300; ++i) {
        timeBinLimits.push_back(i / 100000);
    }

    binDataTimeBinLimits(binRatioAverage, binTimesAverage, binRatioStdDev, binTimesStdDev, binData, timeBinLimits);

    // Plot a graph of time against ratio in one of the bins to show jonas tomorrow
    TGraphErrors *plot = new TGraphErrors(binRatioAverage[0].size(),
                                          binTimesAverage[0].data(),
                                          binRatioAverage[0].data(),
                                          binTimesStdDev[0].data(),
                                          binRatioStdDev[0].data());
    plot->SetTitle("DCS/CF amplitude ratio; time; ratio");
    plot->Draw("*ap");
}
