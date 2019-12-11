#ifndef BIN_GENERATED_DECAYS_HPP
#define BIN_GENERATED_DECAYS_HPP

#include <string>

#include "TLorentzVector.h"
#include "TTree.h"

/*
 * Class representing the data stored in a ROOT file for a series of decays
 *
 */
class DecaysData
{
  public:
    explicit DecaysData(TFile *myTFile, std::string treeName);

    const std::vector<TLorentzVector> particleData(std::string particleName);

    size_t numEvents;

  protected:
    TTree *myTree;

  private:
    void getNumEvents();
};

/*
 * Class for D or Dbar decays to K3pi
 */
class D2K3PiData : public DecaysData
{
    using DecaysData::DecaysData;

  public:
    void setDecayTimes(std::string timesBranchName);

    std::vector<double>         decayTimes{};
    std::vector<TLorentzVector> kVectors{};
    std::vector<TLorentzVector> pi1Vectors{};
    std::vector<TLorentzVector> pi2Vectors{};
    std::vector<TLorentzVector> pi3Vectors{};
};

#endif // BIN_GENERATED_DECAYS_HPP