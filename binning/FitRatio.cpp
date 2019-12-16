#ifndef FITRATIO_CPP
#define FITRATIO_CPP

#include <iostream>
#include <vector>

#include "FitRatio.h"

FitRatio::FitRatio(std::vector<double> &binLimits, std::vector<double> &ratioData, std::vector<double> &ratioErrors)
{
    verifyInputs(binLimits, ratioData, ratioErrors);
}

void FitRatio::verifyInputs(std::vector<double> &binLimits,
                            std::vector<double> &ratioData,
                            std::vector<double> &ratioErrors)
{
    // Check that we have the same number of bins, data points and errors
    size_t numBins   = binLimits.size() - 1;
    size_t numPoints = ratioData.size();
    size_t numErrors = ratioErrors.size();

    if (numBins != numPoints || numBins != numErrors) {
        std::cerr << "Number of bins must match number of datapoints and errors" << std::endl;
    }
}

#endif // FITRATIO_CPP