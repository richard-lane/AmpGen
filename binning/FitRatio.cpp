#ifndef FITRATIO_CPP
#define FITRATIO_CPP

#include <iostream>
#include <vector>

#include "TGraphErrors.h"

#include "FitRatio.h"

FitRatio::FitRatio(std::vector<double> &binLimits, std::vector<double> &ratioData, std::vector<double> &ratioErrors)
{
    verifyInputs(binLimits, ratioData, ratioErrors);

    ratioData   = ratioData;
    ratioErrors = ratioErrors;

    timeData.assign(_numBins, 0);
    timeErrors.shrink_to_fit();
    for (size_t i = 0; i < _numBins; ++i) {
        timeData[i]   = 0.5 * (binLimits[i] + binLimits[i + 1]);
        timeErrors[i] = binLimits[i + 1] - binLimits[i];
    }
}

/*
 * Check that we have the right number of bins, data points and errors
 * BinLimits is a vector of the left edge of each bin plus the rightmost edge of the highest bin
 */
void FitRatio::verifyInputs(std::vector<double> &binLimits,
                            std::vector<double> &ratioData,
                            std::vector<double> &ratioErrors)
{
    // Check that we have the same number of bins, data points and errors
    _numBins         = binLimits.size() - 1;
    size_t numPoints = ratioData.size();
    size_t numErrors = ratioErrors.size();

    if (_numBins != numPoints || _numBins != numErrors) {
        std::cerr << "Number of bins must match number of datapoints and errors" << std::endl;
    }
}

/*
 * Fit the ratio
 */
void fitRatios()
{
    ;
}

#endif // FITRATIO_CPP