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
#include "plottingHelpers.cpp"

void printBinSizes(std::vector<std::vector<double>> binnedTimes)
{
    std::vector<size_t> binSizes(NUM_BINS);
    for (size_t bin = 0; bin < NUM_BINS; ++bin) {
        binSizes[bin] = binnedTimes[bin].size();
        std::cout << "points in bin " << bin << ": " << binSizes[bin] << std::endl;
    }
}

/*
 * Bin the CF and Mixed decays modelled in an AmpGen generated inputFile into phase bins as defined by $BIN_LIMITS
 *
 * This function is far too long but i think its probably ok for now
 *
 */
void bin_generated_decays(TFile *mixedDecays, TFile *favouredDecays)
{
    // Create a D2K3piData class instance to encapsulate the data in our ROOT file
    D2K3PiData MixedData = D2K3PiData(mixedDecays, "DalitzEventList");
    MixedData.populate("D_decayTime");

    D2K3PiData FavouredData = D2K3PiData(favouredDecays, "DalitzEventList");
    FavouredData.populate("Dbar0_decayTime");

    // Make some plots to check that the data from ROOT has been read in correctly
    // plot_things(MixedData.kVectors, MixedData.pi1Vectors, MixedData.pi2Vectors);

    // Perform binning; binned times will be set in MixedData.binnedTimes
    MixedData.binTimes();
    FavouredData.binTimes();

    // Find the number of points in each bin and output to console
    std::cout << "Mixed:" << std::endl;
    printBinSizes(MixedData.binnedTimes);
    std::cout << "Favoured:" << std::endl;
    printBinSizes(FavouredData.binnedTimes);

    // Sort the data in each bin in increasing time order and bin the times based on time bins defined in
    // D2K3PiData.setTimeBins()
    MixedData.sortBinnedTimes();
    MixedData.setTimeBins();
    MixedData.splitTimes();

    FavouredData.sortBinnedTimes();
    FavouredData.setTimeBins();
    FavouredData.splitTimes();

    // Find how many points there are in each time bin
    MixedData.setNumPointsPerTimeBin();
    FavouredData.setNumPointsPerTimeBin();

    // Plot the hist of times in each bin
    for (size_t bin = 0; bin < NUM_BINS; bin++) {
        MixedData.plotBinnedTimes(bin);
        FavouredData.plotBinnedTimes(bin);
    }
}
