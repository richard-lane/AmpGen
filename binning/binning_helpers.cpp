/*
 * Helper functions for k3pi binning
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

/*
 * Split a vector of N vectors into a a vector containing N vectors of vectors, each of which has a maximum of chunkSize
 * elements
 *
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
    size_t                                                           numVectors = myVectors.size();
    std::vector<std::vector<std::vector<std::pair<double, double>>>> splitVectors(
        numVectors, std::vector<std::vector<std::pair<double, double>>>(numBins));

    // In each bin, fill in our vector of vectors according to the bin limits
    for (size_t bin = 0; bin < numVectors; ++bin) {
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
