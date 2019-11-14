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
#define NUM_EVENTS 1000000

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
    std::vector<double> kE = findParticleData(1, 3, eventVectors);

    auto  kCanvas = new TCanvas("kE", "kE", 600, 600);
    TH1D *hist    = new TH1D("kE", "kE", 100, 0, 1);
    hist->FillN(kE.size(), kE.data(), 0);
    hist->Draw();

    // ---- Write them to a root file
    //     Create a new ROOT file
    // TFile outFile("GeneratedDecays.root", "CREATE");

    //     Create a TTree called DalitzEventList
    // TTree *myTree = new TTree("DalitzEventList", "Generated D->K3pi decays");

    //     Write the data to branches on the tree called the right things
}
