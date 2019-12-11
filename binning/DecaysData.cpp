#ifndef DECAYSDATA_CPP
#define DECAYSDATA_CPP

#include "DecaysData.h"
#include "binning_helpers.cpp"
#include "rootfile_helpers.cpp"

#include "TFile.h"

DecaysData::DecaysData(TFile *myTFile, std::string treeName)
{
    myTFile->GetObject(treeName.c_str(), myTree);
    getNumEvents();
}

/*
 * Write the data for a given particle to a vector of TLorentzVectors
 */
const std::vector<TLorentzVector> DecaysData::particleData(std::string particleName)
{
    return writeVector(*myTree, particleName);
}

/*
 * Find how many events in our tree and set it to numEvents
 */
void DecaysData::getNumEvents()
{
    numEvents = myTree->GetEntries();
}

/*
 * Set decay times using the right branch name
 */
void D2K3PiData::setDecayTimes(std::string timesBranchName)
{
    decayTimes = std::vector<double>(numEvents, -1);
    saveBranchToVector(*myTree, timesBranchName, decayTimes);
}

/*
 * Populate a D2K3PiData class with particle data
 */
void D2K3PiData::populate(std::string timesBranchName)
{
    kVectors   = particleData("_1_K~");
    pi1Vectors = particleData("_2_pi#");
    pi2Vectors = particleData("_3_pi#");
    pi3Vectors = particleData("_4_pi~");
    setDecayTimes(timesBranchName);
}

#endif // DECAYSDATA_CPP