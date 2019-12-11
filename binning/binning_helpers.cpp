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
 * Split a (sorted) vector based on binLimits
 * Bin limits define lowest, highest and intermediate points
 * Bin limits should be sorted
 */
std::vector<std::vector<double>> splitVectorWithLimits(std::vector<double> &myVector, std::vector<double> binLimits)
{
    // Bin limits should be sorted
    if (!std::is_sorted(binLimits.begin(), binLimits.end())) {
        std::cout << "Bad time bin limits; should be sorted" << std::endl;
        throw;
    }

    // Vector should be sorted
    if (!std::is_sorted(myVector.begin(), myVector.end())) {
        std::cout << "Vector must be sorted in order to split" << std::endl;
        throw;
    }

    // TODO create a consistent definition of what numbins are
    size_t                           numBins = binLimits.size() - 1;
    std::vector<std::vector<double>> splitVector(numBins);

    size_t currentTimeBin = 0;
    for (size_t i = 0; i < myVector.size(); ++i) {
        // If this time is more than the bin limit, we want to put subsequent points in a higher bin
        // Unless we are already putting points in the highest bin
        while (myVector[i] > binLimits[currentTimeBin] && currentTimeBin < numBins - 1) {
            currentTimeBin += 1;
        }
        splitVector[currentTimeBin].push_back(myVector[i]);
    }

    return splitVector;
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

#endif // BINNING_HELPERS_CPP
