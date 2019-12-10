/*
 * bin_generated_decays.cpp
 * ROOT macro to bin AmpGen generated D -> K3pi decays into predefined bins
 *
 * NOTE: if this fails to build with error "dcs.so: No such file or directory", try restarting ROOT and building again.
 *
 * This script contains no error handling, arg verification or anything to make it work nicely
 */

#include <algorithm>
#include <complex>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TGraphErrors.h"
#include "TRandom.h"

#include "binning_helpers.cpp"
#include "k3pi_binning.h"
#include "plottingHelpers.cpp"
#include "rootfile_helpers.cpp"

// ---- Magic Numbers
// DCS and CF relative amplitude and phase
#define DCS_MAGNITUDE 0.0445
#define DCS_PHASE -3.04

/// Bin limits in phase, centred on zero by construction
#define NUM_BINS 5
#define BIN_LIMITS -39, 0, 43, 180

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
 * This function is far too long but i think its probably ok for now
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

    // Find average and std dev of ratios and times in each subvector
    // These are passed to the binning functions below as OUT args
    std::vector<std::vector<double>> binRatioAverage(NUM_BINS);
    std::vector<std::vector<double>> binTimesAverage(NUM_BINS);
    std::vector<std::vector<double>> binRatioStdDev(NUM_BINS);
    std::vector<std::vector<double>> binTimesStdDev(NUM_BINS);

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
