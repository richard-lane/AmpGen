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

  private:
    void verifyInputs(std::vector<double> &myNumeratorBinLimits,
                      std::vector<size_t> &myNumeratorData,
                      std::vector<double> &myDenominatorBinLimits,
                      std::vector<size_t> &myDenominatorData);

    size_t numBins{0};

    std::vector<size_t> numeratorData{};
    std::vector<size_t> denominatorData{};
    std::vector<double> binLimits{};
    std::vector<double> binRatios{};
};

#endif // DATA_SETS_RATIO_HPP