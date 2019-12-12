#ifndef DATA_SETS_RATIO_CPP
#define DATA_SETS_RATIO_CPP

#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TH1D.h"

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
    // Using == to compare doubles?
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

/*
 * Divide two doubles, returning 0 if either the numerator or denominator are zero
 */
double DataSetsRatio::safeDivide(const double num, const double denom)
{
    if (num == 0 || denom == 0) {
        return 0;
    }
    return num / denom;
}

/*
 * Set the ratio of our numerator and denominator's points in each bin
 * The safeDivide function sets the ratio to 0 if either the numerator or denominator are zero
 */
void DataSetsRatio::_setBinRatios()
{
    // Reassign our vector of ratios to -1
    binRatios.assign(numBins, -1);
    binRatios.shrink_to_fit();

    // Unintelligently divide our elements
    // A good implementation would use std::transform but this is fine
    for (size_t i = 0; i < binRatios.size(); ++i) {
        binRatios[i] = safeDivide(numeratorData[i], denominatorData[i]);
    }
}

/*
 * Plot the ratios of numerator to denominator points in each bin
 */

void DataSetsRatio::plotBinRatios()
{
    TCanvas *c      = new TCanvas();
    TH1D *   MyHist = new TH1D("", "Ratios", numBins, binLimits.data());

    for (size_t i = 0; i < numBins; ++i) {
        MyHist->SetBinContent(i, binRatios[i]);
    }
    MyHist->Draw();
}

#endif // DATA_SETS_RATIO_CPP
