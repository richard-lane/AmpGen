#ifndef DATA_SETS_RATIO_CPP
#define DATA_SETS_RATIO_CPP

#include <iostream>
#include <vector>

#include "DataSetsRatio.h"

DataSetsRatio::DataSetsRatio(std::vector<double> &myNumeratorBinLimits,
                             std::vector<size_t> &myNumeratorData,
                             std::vector<double> &myDenominatorBinLimits,
                             std::vector<size_t> &myDenominatorData)
{
    verifyInputs(myNumeratorBinLimits, myNumeratorData, myDenominatorBinLimits, myDenominatorData);

    // Set bin limits; verifyInputs ensures that the numerator and denominator bin limits are the same so just set it to
    // one of these
    binLimits = myNumeratorBinLimits;

    numeratorData   = myNumeratorData;
    denominatorData = myDenominatorData;
}

/*
 * Check that our numerator and denominator bin limits are the same, and that our datasets have the right number of
 * points.
 */
void DataSetsRatio::verifyInputs(std::vector<double> &myNumeratorBinLimits,
                                 std::vector<size_t> &myNumeratorData,
                                 std::vector<double> &myDenominatorBinLimits,
                                 std::vector<size_t> &myDenominatorData)
{
    // Check that bin limits are the same for both data sets
    // Using == to compare doubles..?
    if (myNumeratorBinLimits != myDenominatorBinLimits) {
        std::cerr << "Cannot compare two datasets with different bin limits" << std::endl;
        throw;
    }

    // Check that the datasets have the right number of points
    // Our bin limits should give us the leftmost point of the lowest bin, all intermediate points and the rightmost
    // point of the highest bin
    // This means we have {length(binVector) - 1} bins
    numBins = myNumeratorBinLimits.size() - 1;
    if ((myNumeratorData.size() != numBins || myDenominatorData.size() != numBins)) {
        std::cerr << "Incompatible bins and datasets provided; must have the same number of bins and datapoints"
                  << std::endl;
        throw;
    }
}

#endif // DATA_SETS_RATIO_CPP
