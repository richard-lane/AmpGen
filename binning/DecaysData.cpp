#ifndef DECAYSDATA_CPP
#define DECAYSDATA_CPP

#include "DecaysData.h"
#include "binning_helpers.cpp"

#include "TFile.h"

/*
 * All the information we need to extract is in the ROOT file on a given Tree, so this is all the constructor needs.
 *
 */
DecaysData::DecaysData(TFile *myTFile, std::string treeName)
{
    myTFile->GetObject(treeName.c_str(), myTree);
    numEvents = myTree->GetEntries();
}


/*
 * Write the data on branchName to the index'th position of each TLorentzVector in myVector.
 * e.g. to write x-momenta of a particle described by ROOT branch foo_Px, call saveBranchToVector("foo_Px", myVector, 0)
 *
 * The TLorentzVector should be of the form (Px, Py, Pz, E).
 */
void DecaysData::writeBranchToLorentzVectors(const std::string &          branchName,
                                             std::vector<TLorentzVector> &myVector,
                                             const size_t &               index)
{
    double myData{0.0};

    myTree->SetBranchAddress(branchName.c_str(), &myData);
    for (Long64_t i = 0; i < myTree->GetEntries(); ++i) {
        myTree->GetEntry(i);
        myVector[i][index] = myData;
    }

    // Reset all branch addresses to avoid a bug where repeatedly calling this function would set an array to the wrong
    // values
    myTree->ResetBranchAddresses();
}

/*
 * Write the data for a given particle to a vector of TLorentzVectors
 *
 */
const std::vector<TLorentzVector> DecaysData::particleData(std::string particleName)
{
    std::vector<TLorentzVector> myVector = std::vector<TLorentzVector>(numEvents);

    writeBranchToLorentzVectors(particleName + "_Px", myVector, 0);
    writeBranchToLorentzVectors(particleName + "_Py", myVector, 1);
    writeBranchToLorentzVectors(particleName + "_Pz", myVector, 2);
    writeBranchToLorentzVectors(particleName + "_E", myVector, 3);

    return myVector;
}

/*
 * Set decay times on a branch
 * timesBranchName might be set to e.g. "D_decayTime"
 *
 */
void D2K3PiData::setDecayTimes(std::string timesBranchName)
{

    double myData{0.0};
    decayTimes = std::vector<double>(numEvents, -1);

    myTree->SetBranchAddress(timesBranchName.c_str(), &myData);
    for (Long64_t i = 0; i < myTree->GetEntries(); ++i) {
        myTree->GetEntry(i);
        decayTimes[i] = myData;
    }

    // Reset all branch addresses to avoid a bug where repeatedly calling this function would set an array to the wrong
    // values
    myTree->ResetBranchAddresses();
}

/*
 * Populate a D2K3PiData class with particle data
 *
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