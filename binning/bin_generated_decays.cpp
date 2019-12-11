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
#include "TH1D.h"
#include "TRandom.h"

#include "DecaysData.cpp"
#include "DecaysData.h"
#include "binning_helpers.cpp"
#include "plottingHelpers.cpp"

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
 * Bin the CF and Mixed decays modelled in an AmpGen generated inputFile into phase bins as defined by $BIN_LIMITS
 *
 * This function is far too long but i think its probably ok for now
 *
 */
void bin_generated_decays(TFile *mixedDecays)
{
    // Create a D2K3piData class instance to encapsulate the data in our ROOT file
    D2K3PiData MixedData = D2K3PiData(mixedDecays, "DalitzEventList");
    MixedData.populate("D_decayTime");

    // Make some plots to check that the data from ROOT has been read in correctly
    plot_things(MixedData.kVectors, MixedData.pi1Vectors, MixedData.pi2Vectors);

    // Perform binning; binned times will be set in MixedData.binnedTimes
    MixedData.binTimes();

    // Find the number of points in each bin and output to console
    std::vector<size_t> binSizes(NUM_BINS);
    for (size_t bin = 0; bin < NUM_BINS; ++bin) {
        binSizes[bin] = MixedData.binnedTimes[bin].size();
        std::cout << "points in bin " << bin << ": " << binSizes[bin] << std::endl;
    }

    // Sort the data in each bin in increasing time order
    sortVectorOfVectors(MixedData.binnedTimes);

    // Bin data into time bins defined by a vector
    std::vector<double> timeBinLimits{};
    for (double i = 0; i < 300; ++i) {
        timeBinLimits.push_back(i / 100000);
    }

    std::vector<std::vector<std::vector<double>>> timeBinnedData(NUM_BINS);
    for (size_t bin = 0; bin < NUM_BINS; ++bin) {
        timeBinnedData[bin] = splitVectorWithLimits(MixedData.binnedTimes[bin], timeBinLimits);
    }

    // Find how many points there are in each time bin
    std::vector<std::vector<size_t>> numPointsPerTimeBin(NUM_BINS);

    for (size_t bin = 0; bin < NUM_BINS; ++bin) {
        for (size_t i = 0; i < timeBinnedData[bin].size(); ++i) {
            numPointsPerTimeBin[bin].push_back(timeBinnedData[bin][i].size());
        }
    }

    TH1D *MyHist = new TH1D("foo", "foo", timeBinLimits.size() - 1, timeBinLimits.data());
    for (size_t i = 0; i < timeBinLimits.size() - 1; ++i) {
        MyHist->SetBinContent(i, numPointsPerTimeBin[0][i]);
    }
    MyHist->Draw();
}
