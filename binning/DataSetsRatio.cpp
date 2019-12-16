#ifndef DATA_SETS_RATIO_CPP
#define DATA_SETS_RATIO_CPP

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TGraphErrors.h"

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

    // Find also the errors in our ratios.
    _setBinRatioErrors();
}

/*
 * Assuming the error in our counts is sqrt(count) and that the fractional error in a ratio is given by adding count
 * errors in quadrature, find the error in a ratio.
 *
 * Sets errors for ratio=0 to 0
 */
double DataSetsRatio::ratioError(const double &ratio, const size_t &numeratorCounts, const size_t &denominatorCounts)
{
    // If our ratio is zero, our error should also be
    if (ratio == 0) {
        return 0;
    }

    // Otherwise our ratio is found by assuming delta(Counts) = sqrt(Counts)
    return std::sqrt(((double)numeratorCounts + (double)denominatorCounts) /
                     ((double)numeratorCounts * (double)denominatorCounts)) *
           ratio;
}

/*
 * Set the values of binRatioErrors
 */
void DataSetsRatio::_setBinRatioErrors()
{
    binRatioErrors.assign(numBins, 0);
    binRatioErrors.shrink_to_fit();

    for (size_t i = 0; i < numBins; ++i) {
        binRatioErrors[i] = ratioError(binRatios[i], numeratorData[i], denominatorData[i]);
    }
}

/*
 * Plot the ratios of numerator to denominator points in each bin
 * Doesn't actually draw the plot but assigns _ratioPlot to a new TGraphErrors object
 */
void DataSetsRatio::plotBinRatios()
{
    std::vector<double> timeBinCentres(numBins);
    std::vector<double> zeros(numBins, 0);

    for (size_t i = 0; i < numBins; ++i) {
        timeBinCentres[i] = 0.5 * (binLimits[i] + binLimits[i + 1]);
    }

    _ratioPlot =
        new TGraphErrors(numBins, timeBinCentres.data(), binRatios.data(), zeros.data(), binRatioErrors.data());
}

/*
 * Fit a second order polynomial to a plot of ratios in each bin
 *
 * @param draw: whether to draw the graph or just fit to its data.
 */
void DataSetsRatio::fitToData(bool draw)
{
    plotBinRatios();
    if (draw) {
        TCanvas *c = new TCanvas();
        _ratioPlot->Draw("*ap");
    }
    _ratioPlot->Fit("pol1");
}

#endif // DATA_SETS_RATIO_CPP
