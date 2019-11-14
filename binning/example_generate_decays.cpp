#include <cmath>
#include <complex>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TRandom.h"
#include "TTree.h"

#include "k3pi_binning.h"

// ---- Magic Numbers
#define NUM_EVENTS 1000

#define NUMBER_PRODUCTS 4
#define K_MASS 0.493677
#define PION_MASS 0.13957018

#define PARTICLE_MOMENTUM 0, 0, 0, 1.86962 // Decaying particle momentum in GeV

/*
 * From a vector of TLorentzVectors, return a vector of a particle's properties
 *
 * e.g. findParticleData(i,j,myVector) will return a vector of the ith particle's j'th momentum component
 *
 * Uses the convention that 4-momentum is (px, py, pz, E)
 */
std::vector<double> findParticleData(const size_t                                    particleIndex,
                                     const size_t                                    momentumIndex,
                                     const std::vector<std::vector<TLorentzVector>> &eventsVector)
{
    size_t vectorLength = eventsVector.size();

    std::vector<double> dataVector = std::vector<double>(vectorLength);

    for (size_t i = 0; i < vectorLength; ++i) {
        dataVector[i] = eventsVector[i][particleIndex][momentumIndex];
    }

    return dataVector;
}

void example_generate_decays()
{
    // Create a phase space object to generate decays from
    TLorentzVector      dMomentum(PARTICLE_MOMENTUM);
    std::vector<double> daughterMasses = {K_MASS, PION_MASS, PION_MASS, PION_MASS};

    TGenPhaseSpace phsp;
    phsp.SetDecay(dMomentum, NUMBER_PRODUCTS, daughterMasses.data());

    // ---- Create decays using the k3pi scheme
    std::vector<std::vector<TLorentzVector>> eventVectors = std::vector<std::vector<TLorentzVector>>(NUM_EVENTS);
    for (size_t i = 0; i < NUM_EVENTS; ++i) {
        eventVectors[i] = k3pi_binning::makeUnweighted(phsp);
    }

    // Create vectors of K and pi data
    std::vector<double> kE  = findParticleData(0, 3, eventVectors);
    std::vector<double> kPx = findParticleData(0, 0, eventVectors);
    std::vector<double> kPy = findParticleData(0, 1, eventVectors);
    std::vector<double> kPz = findParticleData(0, 2, eventVectors);

    std::vector<double> pi1E  = findParticleData(1, 3, eventVectors);
    std::vector<double> pi1Px = findParticleData(1, 0, eventVectors);
    std::vector<double> pi1Py = findParticleData(1, 1, eventVectors);
    std::vector<double> pi1Pz = findParticleData(1, 2, eventVectors);

    std::vector<double> pi2E  = findParticleData(2, 3, eventVectors);
    std::vector<double> pi2Px = findParticleData(2, 0, eventVectors);
    std::vector<double> pi2Py = findParticleData(2, 1, eventVectors);
    std::vector<double> pi2Pz = findParticleData(2, 2, eventVectors);

    std::vector<double> pi3E  = findParticleData(3, 3, eventVectors);
    std::vector<double> pi3Px = findParticleData(3, 0, eventVectors);
    std::vector<double> pi3Py = findParticleData(3, 1, eventVectors);
    std::vector<double> pi3Pz = findParticleData(3, 2, eventVectors);

    // ---- Write them to a root file
    // This file must be initialised before the tree is created so that ROOT knows which file to write the tree to
    TFile outFile("GeneratedDecays.root", "RECREATE");

    //     Create a TTree called DalitzEventList
    TTree *myTree = new TTree("DalitzEventList", "Generated D->K3pi decays");

    //     Write the data to branches on the tree called the right things
    unsigned long long bufsize = sizeof(double) * eventVectors.size();
    myTree->Branch("_1_K~_E", &kE);
    myTree->Branch("_1_K~_Px", &kPx);
    myTree->Branch("_1_K~_Py", &kPy);
    myTree->Branch("_1_K~_Pz", &kPz);

    myTree->Branch("_2_pi#_Px", &pi1Px);
    myTree->Branch("_2_pi#_Py", &pi1Py);
    myTree->Branch("_2_pi#_Pz", &pi1Pz);
    myTree->Branch("_2_pi#_E", &pi1E);

    myTree->Branch("_3_pi#_Px", &pi2Px);
    myTree->Branch("_3_pi#_Py", &pi2Py);
    myTree->Branch("_3_pi#_Pz", &pi2Pz);
    myTree->Branch("_3_pi#_E", &pi2E);

    myTree->Branch("_4_pi~_Px", &pi3Px);
    myTree->Branch("_4_pi~_Py", &pi3Py);
    myTree->Branch("_4_pi~_Pz", &pi3Pz);
    myTree->Branch("_4_pi~_E", &pi3E);

    myTree->Fill();

    outFile.Write();
    outFile.Close();
}
