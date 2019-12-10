#include "TLorentzVector.h"
#include "TTree.h"

/*
 * Write the data on branchName to myVector.
 */
void saveBranchToVector(TTree &myTree, const std::string &branchName, std::vector<double> &myVector)
{
    double myData{0.0};

    myTree.SetBranchAddress(branchName.c_str(), &myData);
    for (Long64_t i = 0; i < myTree.GetEntries(); ++i) {
        myTree.GetEntry(i);
        myVector[i] = myData;
    }
    // Reset all branch addresses to avoid a bug where repeatedly calling this function would set an array to the wrong
    // values
    myTree.ResetBranchAddresses();
}

/*
 * Write the data on branchName to the index'th position of each TLorentzVector in myVector.
 * e.g. to write x-momenta of a particle described by ROOT branch foo_Px, call saveBranchToVector("foo_Px", myVector, 0)
 *
 * The TLorentzVector should be of the form (Px, Py, Pz, E).
 */
void writeBranchToLorentzVectors(TTree &                      myTree,
                                 const std::string &          branchName,
                                 std::vector<TLorentzVector> &myVector,
                                 const size_t &               index)
{
    double myData{0.0};

    myTree.SetBranchAddress(branchName.c_str(), &myData);
    for (Long64_t i = 0; i < myTree.GetEntries(); ++i) {
        myTree.GetEntry(i);
        myVector[i][index] = myData;
    }
    // Reset all branch addresses to avoid a bug where repeatedly calling this function would set an array to the wrong
    // values
    myTree.ResetBranchAddresses();
}

/*
 * Write a vector of TLorentzVectors containing the data for the particle branchName
 */
std::vector<TLorentzVector> writeVector(TTree &myTree, const std::string &particleName)
{
    std::vector<TLorentzVector> myVector = std::vector<TLorentzVector>(myTree.GetEntries());

    writeBranchToLorentzVectors(myTree, particleName + "_Px", myVector, 0);
    writeBranchToLorentzVectors(myTree, particleName + "_Py", myVector, 1);
    writeBranchToLorentzVectors(myTree, particleName + "_Pz", myVector, 2);
    writeBranchToLorentzVectors(myTree, particleName + "_E", myVector, 3);

    return myVector;
}
