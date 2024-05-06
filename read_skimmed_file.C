/***************************************************************************
 *  Read the skimmed version of delphes and print data from the TLorentzvectors 
****************************************************************************/

#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

void readRootFile(const char* filename) {
    // Open the root file
    TFile inputFile(filename, "READ");
    if (!inputFile.IsOpen()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // Get the TTree from the root file
    TTree* ntuple = dynamic_cast<TTree*>(inputFile.Get("Ntuple"));
    if (!ntuple) {
        std::cerr << "Error: Could not find TTree 'Ntuple' in file " << filename << std::endl;
        return;
    }

    // Create TClonesArrays to store TLorentzVectors
    TClonesArray* jetArray = nullptr;
    TClonesArray* bjetArray = nullptr;
    TClonesArray* electronArray = nullptr;
    TClonesArray* muonArray = nullptr;
    TClonesArray* photonArray = nullptr;
    TClonesArray* METArray = nullptr;
    Double32_t cmsEnergy = 0;

    // Set branch addresses for each TClonesArray
    ntuple->SetBranchAddress("Jets", &jetArray);
    ntuple->SetBranchAddress("BJets", &bjetArray);
    ntuple->SetBranchAddress("Electrons", &electronArray);
    ntuple->SetBranchAddress("Muons", &muonArray);
    ntuple->SetBranchAddress("Photons", &photonArray);
    ntuple->SetBranchAddress("MET", &METArray);
    ntuple->SetBranchAddress("CMS_Energy", &cmsEnergy);

    std::vector<Double_t>* weightVector = nullptr;
    ntuple->SetBranchAddress("Weight", &weightVector);

    // Loop over all entries in the TTree
    Long64_t numEntries = ntuple->GetEntries();
    for (Long64_t i = 0; i < numEntries; ++i) {
        ntuple->GetEntry(i);

        std::cout << "Event " << i << ":" << std::endl;

        std::cout << " CMS_Energy: " << cmsEnergy << std::endl;

        // Print TLorentzVectors for jets
        std::cout << " Jets:" << std::endl;
        for (int j = 0; j < jetArray->GetEntries(); ++j) {
            TLorentzVector* jet = dynamic_cast<TLorentzVector*>(jetArray->At(j));
            std::cout << "   Jet " << j << ": (Pt, Eta, Phi, M) = (" << jet->Pt() << ", " << jet->Eta() << ", " << jet->Phi() << ", " << jet->M() << ")" << std::endl;
        }

        // Print TLorentzVectors for b-jets
        std::cout << "  B-Jets:" << std::endl;
        for (int j = 0; j < bjetArray->GetEntries(); ++j) {
            TLorentzVector* bjet = dynamic_cast<TLorentzVector*>(bjetArray->At(j));
            std::cout << "   B-Jet " << j << ": (Pt, Eta, Phi, M) = (" << bjet->Pt() << ", " << bjet->Eta() << ", " << bjet->Phi() << ", " << bjet->M() << ")" << std::endl;
        }

        // Print TLorentzVectors for electrons
        std::cout << " Electrons:" << std::endl;
        for (int j = 0; j < electronArray->GetEntries(); ++j) {
            TLorentzVector* electron = dynamic_cast<TLorentzVector*>(electronArray->At(j));
            std::cout << "   Electron " << j << ": (Pt, Eta, Phi, M) = (" << electron->Pt() << ", " << electron->Eta() << ", " << electron->Phi() << ", " << electron->M() << ")" << std::endl;
        }

        // Print TLorentzVectors for muons
        std::cout << " Muons:" << std::endl;
        for (int j = 0; j < muonArray->GetEntries(); ++j) {
            TLorentzVector* muon = dynamic_cast<TLorentzVector*>(muonArray->At(j));
            std::cout << "   Muon " << j << ": (Pt, Eta, Phi, M) = (" << muon->Pt() << ", " << muon->Eta() << ", " << muon->Phi() << ", " << muon->M() << ")" << std::endl;
        }

        // Print TLorentzVectors for photons
        std::cout << " Photons:" << std::endl;
        for (int j = 0; j < photonArray->GetEntries(); ++j) {
            TLorentzVector* photon = dynamic_cast<TLorentzVector*>(photonArray->At(j));
            std::cout << "   Photon " << j << ": (Pt, Eta, Phi, M) = (" << photon->Pt() << ", " << photon->Eta() << ", " << photon->Phi() << ", " << photon->M() << ")" << std::endl;
        }

        // Print TLorentzVectors for photons
        std::cout << " METs:" << std::endl;
        for (int j = 0; j < METArray->GetEntries(); ++j) {
            TLorentzVector* MET = dynamic_cast<TLorentzVector*>(METArray->At(j));
            std::cout << "   MET " << j << ": (Pt, Eta, Phi, M) = (" << MET->Pt() << ", " << MET->Eta() << ", " << MET->Phi() << ", " << MET->M() << ")" << std::endl;
        }

        // Print Values of Weights
        std::cout << " Weights:" << std::endl;
        if (weightVector) {
            for (size_t j = 0; j < weightVector->size(); ++j) {
                std::cout << "  Weight " << j << ": " << (*weightVector)[j] << endl;
            }
        } else {
            std::cerr << "Error: Weight vector is null." << std::endl;
        }

        std::cout << std::endl;
    }

    // Close the input file
    inputFile.Close();
}

void read_skimmed_file(const char* inputFileName) {
    // Call the function to read the root file and print TLorentzVectors
    readRootFile(inputFileName);

    return 0;
}


