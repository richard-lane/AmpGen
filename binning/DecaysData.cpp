#ifndef DECAYSDATA_CPP
#define DECAYSDATA_CPP

#include <iostream>

#include "DecaysData.h"

#include "TFile.h"
#include "k3pi_binning.h"

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

/*
 * Sort the events into bins based on $BIN_LIMITS
 *
 * Sets the binnedTimes vector of vectors
 */
void D2K3PiData::binTimes(void)
{
    // Define our bins
    k3pi_binning::binning bins(dcsFile, cfFile, dcs_offset, {BIN_LIMITS});

    // Iterate over all events, sorting them into bins
    for (size_t i = 0; i < numEvents; ++i) {
        // Create a vector of TLorentzVectors for this event (K+, pi-, pi-, pi+)
        std::vector<TLorentzVector> eventVector{kVectors[i], pi1Vectors[i], pi2Vectors[i], pi3Vectors[i]};
        auto                        event = k3pi_binning::eventFromVectors(eventVector);

        // Find which bin the event belongs in
        // The 1 tags the sign of the K meson in the D->K3pi decay
        int bin = bins.bin(eventVector, 1);

        // Log the time of the event in this bin
        // This might be slow because it's dynamically resizing the vector, but should be ok for our purposes.
        binnedTimes[bin].push_back(decayTimes[i]);
    }
}

/*
 * Sort our vectors of binned times into increasing order
 */
void D2K3PiData::sortBinnedTimes()
{
    for (size_t i = 0; i < NUM_BINS; ++i) {
        std::sort(binnedTimes[i].begin(), binnedTimes[i].end());
    }
}

/*
 * Set time bins
 */
void D2K3PiData::setTimeBins()
{

    for (double i = 0; i < 341; ++i) {
        timeBinLimits.push_back(i / 100000);
    }
}

/*
 * Split a (sorted) vector based on timeBinLimits (it wouldn't be hard to make this work with unsorted vectors, though
 * it's slightly more efficient if we demand that vectors are sorted)
 *
 * Bin limits define lowest, highest and intermediate points
 * Bin limits should be sorted
 */
std::vector<std::vector<double>> D2K3PiData::splitVectorWithLimits(std::vector<double> &myVector)
{
    // Bin limits should be sorted
    if (!std::is_sorted(timeBinLimits.begin(), timeBinLimits.end())) {
        std::cerr << "Bad time bin limits; should be sorted" << std::endl;
        throw;
    }

    // Vector should be sorted
    if (!std::is_sorted(myVector.begin(), myVector.end())) {
        std::cerr << "Vector must be sorted in order to split" << std::endl;
        throw;
    }

    // First point in the vector should be in a bin; if it is smaller than the lowest bin emit a warning
    if (myVector[0] < timeBinLimits[0]) {
        std::cerr << "[warning] Vector contains values smaller than lowest bin limit" << std::endl;
    }

    // Last point in the vector should be in a bin; throw if it is larger than the highest bin
    size_t vectorSize = myVector.size();
    size_t numBins    = timeBinLimits.size() - 1;
    if (myVector[vectorSize - 1] > timeBinLimits[numBins]) {
        std::cerr << "Vector extends over range larger than largest bin limit" << std::endl;
        throw;
    }

    // Split our vector into smaller vectors along the bin limits
    std::vector<std::vector<double>> splitVector(numBins);

    // Left edge of bin we are currently inserting points into
    size_t currentTimeBin = 0;
    for (size_t i = 0; i < vectorSize; ++i) {
        // If this value is more than the bin limit, we want to put subsequent points in a higher bin
        while (myVector[i] > timeBinLimits[currentTimeBin + 1]) {
            currentTimeBin += 1;

            // I think this code never gets called...?
            if (currentTimeBin > numBins) {
                std::cerr << "Attempted to insert points into a bin that does not exist" << std::endl;
                throw;
            }
        }
        splitVector[currentTimeBin].push_back(myVector[i]);
    }

    return splitVector;
}

/*
 * Populate timeBinnedTimes with the times binned along timeBinLimits
 */
void D2K3PiData::splitTimes()
{
    for (size_t bin = 0; bin < NUM_BINS; ++bin) {
        timeBinnedTimes[bin] = splitVectorWithLimits(binnedTimes[bin]);
    }
}

#endif // DECAYSDATA_CPP