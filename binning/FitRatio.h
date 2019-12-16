#ifndef FITRATIO_H
#define FITRATIO_H

#include <vector>

/*
 * Class for fitting the ratio of CF to DCS decays
 *
 */
class FitRatio
{
  public:
    explicit FitRatio(std::vector<double> &binLimits, std::vector<double> &ratioData, std::vector<double> &ratioErrors);

    void fitRatios();

  private:
    void verifyInputs(std::vector<double> &binLimits, std::vector<double> &ratioData, std::vector<double> &ratioErrors);

    size_t _numBins = 0;

    std::vector<double> ratioData{};
    std::vector<double> ratioErrors{};
    std::vector<double> timeData{};
    std::vector<double> timeErrors{};
};

#endif // FITRATIO_H