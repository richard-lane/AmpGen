#ifndef DATA_SETS_RATIO_HPP
#define DATA_SETS_RATIO_HPP

#include <vector>

#include "TGraph.h"

/*
 * Class for calculating and representing the ratio of two binned datasets
 *
 * A ratio of the two datasets will be calculated, so one is referred to as "numerator" and one as "denominator"
 *
 */
class DataSetsRatio
{
  public:
    explicit DataSetsRatio(std::vector<double> &myNumeratorBinLimits,
                           std::vector<size_t> &myNumeratorData,
                           std::vector<double> &myDenominatorBinLimits,
                           std::vector<size_t> &myDenominatorData);

    void _setBinRatios();
    void fitToData(bool draw, std::string plotTitle);

    std::vector<double> binRatios{};
    std::vector<double> binRatioErrors{};
    std::vector<double> binLimits{};
    std::vector<double> binCentres{};
    std::vector<double> binErrors{};

  private:
    void verifyInputs(std::vector<double> &myNumeratorBinLimits,
                      std::vector<size_t> &myNumeratorData,
                      std::vector<double> &myDenominatorBinLimits,
                      std::vector<size_t> &myDenominatorData);
    void _setBinRatioErrors();
    void plotBinRatios();

    double ratioError(const double &ratio, const size_t &numeratorCounts, const size_t &denominatorCounts);

    TGraph *_ratioPlot = nullptr;

    size_t numBins{0};

    // Once we remove any bins that have inf/NaN/zeroes in their data, our number of data points will be smaller than
    // the number of bins.
    size_t numPoints{0};

    std::vector<size_t> numeratorData{};
    std::vector<size_t> denominatorData{};
};

#endif // DATA_SETS_RATIO_HPP