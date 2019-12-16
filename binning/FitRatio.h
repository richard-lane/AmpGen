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

    void _setBinRatios();
    void plotBinRatios();
    void fitRatios();

  private:
    void verifyInputs(std::vector<double> &binLimits, std::vector<double> &ratioData, std::vector<double> &ratioErrors);

    std::vector<double> ratioData{};
    std::vector<double> ratioErrors{};
    std::vector<double> binLimits{};
};

#endif // FITRATIO_H