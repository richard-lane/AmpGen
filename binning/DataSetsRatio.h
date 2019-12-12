#ifndef DATA_SETS_RATIO_HPP
#define DATA_SETS_RATIO_HPP

#include <vector>

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
    void plotBinRatios();

    std::vector<double> binRatios{};

  private:
    void   verifyInputs(std::vector<double> &myNumeratorBinLimits,
                        std::vector<size_t> &myNumeratorData,
                        std::vector<double> &myDenominatorBinLimits,
                        std::vector<size_t> &myDenominatorData);
    double safeDivide(const double num, const double denom);

    size_t numBins{0};

    std::vector<size_t> numeratorData{};
    std::vector<size_t> denominatorData{};
    std::vector<double> binLimits{};
};

#endif // DATA_SETS_RATIO_HPP