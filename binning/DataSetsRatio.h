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
    explicit DataSetsRatio(std::vector<double> numeratorBinLimits,
                           std::vector<size_t> numeratorData,
                           std::vector<double> denominatorBinLimits,
                           std::vector<size_t> denominatorData);

    void _setBinRatios(std::vector<size_t> numeratorData, std::vector<size_t> denominatorData);

  private:
    void verifyInputs(std::vector<double> numeratorBinLimits,
                      std::vector<size_t> numeratorData,
                      std::vector<double> denominatorBinLimits,
                      std::vector<size_t> denominatorData);
    void _setBinLimits(std::vector<double> myBinLimits);
    void _setData();

    std::vector<size_t> numeratorData{};
    std::vector<size_t> denominatorData{};
    std::vector<double> binLimits{};
    std::vector<double> binRatios{};
};

#endif // DATA_SETS_RATIO_HPP