/*
 * Helper functions for k3pi binning
 */
#ifndef BINNING_HELPERS_CPP
#define BINNING_HELPERS_CPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

/*
 * Split a (sorted) vector based on binLimits (it wouldn't be hard to make this work with unsorted vectors, though it's
 * slightly more efficient if we demand that vectors are sorted)
 *
 * Bin limits define lowest, highest and intermediate points
 * Bin limits should be sorted
 */
std::vector<std::vector<double>> splitVectorWithLimits(std::vector<double> &myVector, std::vector<double> binLimits)
{
    // Bin limits should be sorted
    if (!std::is_sorted(binLimits.begin(), binLimits.end())) {
        std::cerr << "Bad time bin limits; should be sorted" << std::endl;
        throw;
    }

    // Vector should be sorted
    if (!std::is_sorted(myVector.begin(), myVector.end())) {
        std::cerr << "Vector must be sorted in order to split" << std::endl;
        throw;
    }

    // First point in the vector should be in a bin; if it is smaller than the lowest bin emit a warning
    if (myVector[0] < binLimits[0]) {
        std::cerr << "[warning] Vector contains values smaller than lowest bin limit" << std::endl;
    }

    // Last point in the vector should be in a bin; throw if it is larger than the highest bin
    size_t vectorSize = myVector.size();
    size_t numBins    = binLimits.size() - 1;
    if (myVector[vectorSize - 1] > binLimits[numBins]) {
        std::cerr << "Vector extends over range larger than largest bin limit" << std::endl;
        throw;
    }

    // Split our vector into smaller vectors along the bin limits
    std::vector<std::vector<double>> splitVector(numBins);

    // Left edge of bin we are currently inserting points into
    size_t currentTimeBin = 0;
    for (size_t i = 0; i < vectorSize; ++i) {
        // If this time is more than the bin limit, we want to put subsequent points in a higher bin
        // Unless we are already putting points in the highest bin
        while (myVector[i] > binLimits[currentTimeBin + 1] && currentTimeBin < numBins - 1) {
            currentTimeBin += 1;

            // I think this code never gets called...?
            if (currentTimeBin > numBins) {
                std::cerr << "Attempted to insert points into a bin that does not exist" << std::endl;
                throw;
            }
        }
        splitVector[currentTimeBin].push_back(myVector[i]);
    }

    return splitVector;
}

#endif // BINNING_HELPERS_CPP
