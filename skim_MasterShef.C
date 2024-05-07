/***************************************************************************
 *  Read MasterShef file and save a skimmed version of MasterShef (root) file with TLorentzvectors
****************************************************************************/

#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <iostream>
#include <TLeaf.h>

void readEventData(const char* inputFileName, const char* outputFileName, Double32_t cmsEnergy) {
    // Open the ROOT file
    TFile inputFile(inputFileName);
    if (inputFile.IsZombie()) {
        std::cerr << "Error: Unable to open file " << inputFileName << std::endl;
        return;
    }
    
    // Access the Delphes TTree
    TTree* tree = dynamic_cast<TTree*>(inputFile.Get("Nominal"));
    if (!tree) {
        std::cerr << "Error: Unable to access tree Nominal in file " << inputFileName << std::endl;
        return;
    }
    
    
    vector<float>*  jet_pt = 0;
    vector<float>   *jet_eta = 0;
    vector<float>   *jet_phi = 0;
    vector<float>   *jet_e = 0;
    vector<int>     *jet_IsBJet = 0;
    
    vector<float>   *el_pt = 0;
    vector<float>   *el_eta = 0;
    vector<float>   *el_phi = 0;
    vector<float>   *el_e = 0;
    
    vector<float>   *mu_pt = 0;
    vector<float>   *mu_eta = 0;
    vector<float>   *mu_phi = 0;
    vector<float>   *mu_e = 0;
    
    vector<float>   *ph_pt = 0;
    vector<float>   *ph_eta = 0;
    vector<float>   *ph_phi = 0;
    vector<float>   *ph_e = 0;
    
    Float_t         MET_pt = 0;
    Float_t         MET_phi = 0;
    Float_t         MET_pt_prime = 0;
    Float_t         MET_phi_prime = 0;
    Float_t         nominalWeight;
    
    Int_t           nJets;
    Int_t           nEl;
    Int_t           nMu;
    Int_t           nPh;
    
    TBranch        *b_jet_pt;   //!
    TBranch        *b_jet_eta;   //!
    TBranch        *b_jet_phi;   //!
    TBranch        *b_jet_e;   //!
    TBranch        *b_jet_IsBJet;   //!
    
    TBranch        *b_el_pt;   //!
    TBranch        *b_el_eta;   //!
    TBranch        *b_el_phi;   //!
    TBranch        *b_el_e;   //!
    
    TBranch        *b_mu_pt;   //!
    TBranch        *b_mu_eta;   //!
    TBranch        *b_mu_phi;   //!
    TBranch        *b_mu_e;   //!
    //!
    TBranch        *b_ph_pt;   //!
    TBranch        *b_ph_eta;   //!
    TBranch        *b_ph_phi;   //!
    TBranch        *b_ph_e;   //!
    
    TBranch        *b_MET_pt;   //!
    TBranch        *b_MET_phi;   //!
    TBranch        *b_MET_pt_prime;   //!
    TBranch        *b_MET_phi_prime;   //!
    
    TBranch        *b_nJets;   //!
    TBranch        *b_nEl;   //!
    TBranch        *b_nMu;   //!
    TBranch        *b_nPh;   //!
    
    TBranch        *b_nominalWeight;   //!

    
    tree->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
    tree->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
    tree->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
    tree->SetBranchAddress("jet_e", &jet_e, &b_jet_e);
    tree->SetBranchAddress("jet_IsBJet", &jet_IsBJet, &b_jet_IsBJet);
    
    tree->SetBranchAddress("el_pt", &el_pt, &b_el_pt);
    tree->SetBranchAddress("el_eta", &el_eta, &b_el_eta);
    tree->SetBranchAddress("el_phi", &el_phi, &b_el_phi);
    tree->SetBranchAddress("el_e", &el_e, &b_el_e);
    
    tree->SetBranchAddress("mu_pt", &mu_pt, &b_mu_pt);
    tree->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
    tree->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
    tree->SetBranchAddress("mu_e", &mu_e, &b_mu_e);
    
    tree->SetBranchAddress("ph_pt", &ph_pt, &b_ph_pt);
    tree->SetBranchAddress("ph_eta", &ph_eta, &b_ph_eta);
    tree->SetBranchAddress("ph_phi", &ph_phi, &b_ph_phi);
    tree->SetBranchAddress("ph_e", &ph_e, &b_ph_e);
    
    tree->SetBranchAddress("MET_pt", &MET_pt, &b_MET_pt);
    tree->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
    tree->SetBranchAddress("MET_pt_prime", &MET_pt_prime, &b_MET_pt_prime);
    tree->SetBranchAddress("MET_phi_prime", &MET_phi_prime, &b_MET_phi_prime);
    
    tree->SetBranchAddress("nJets", &nJets, &b_nJets);
    tree->SetBranchAddress("nEl", &nEl, &b_nEl);
    tree->SetBranchAddress("nMu", &nMu, &b_nMu);
    tree->SetBranchAddress("nPh", &nPh, &b_nPh);
    
    tree->SetBranchAddress("nominalWeight", &nominalWeight, &b_nominalWeight);

    // Loop over entries in the TTree and create Vectors
    std::vector<std::vector<TLorentzVector>> jet_Vectors;
    std::vector<std::vector<TLorentzVector>> bjet_Vectors;
    std::vector<std::vector<TLorentzVector>> el_Vectors;
    std::vector<std::vector<TLorentzVector>> mu_Vectors;
    std::vector<std::vector<TLorentzVector>> photon_Vectors;
    std::vector<std::vector<TLorentzVector>> MET_Vectors;
    std::vector<std::vector<Double_t>> Weight_Vectors;

    Long64_t numEntries = tree->GetEntries();
    cout << "no. of events :" << numEntries << endl;
    for (Long64_t i = 0; i < numEntries; ++i) {
        // Clear the vector for the current event
        //jet_Vectors.clear();
        std::vector<TLorentzVector> event_JetVectors;
        std::vector<TLorentzVector> event_bJetVectors;
        std::vector<TLorentzVector> event_elVectors;
        std::vector<TLorentzVector> event_muVectors;
        std::vector<TLorentzVector> event_phVectors;
        std::vector<TLorentzVector> event_METVectors;
        std::vector<Double_t> event_WeightVectors;

        // Get entry i
        tree->GetEntry(i);
        //cout << "Event " << i << endl;
        if(i%1000==0) cout << "Event #" <<  i << endl;
        
        for (Int_t j = 0; j < nJets; ++j) {
            TLorentzVector jet_Vector;
            //cout << "jet_pt :" << jet_pt->at(j) << endl;
            jet_Vector.SetPtEtaPhiE(jet_pt->at(j)/1000.0, jet_eta->at(j), jet_phi->at(j), jet_e->at(j)/1000.0);
	    event_JetVectors.push_back(jet_Vector);
        }
        jet_Vectors.push_back(event_JetVectors);
 
        
        for (Int_t j = 0; j < nJets; ++j) {
            TLorentzVector bjet_Vector;
            //cout << j << " isbTag :" << Jet_BTag[j] << endl;
            if ((jet_IsBJet->at(j))==1){
                bjet_Vector.SetPtEtaPhiE(jet_pt->at(j)/1000.0, jet_eta->at(j), jet_phi->at(j), jet_e->at(j)/1000.0);

                event_bJetVectors.push_back(bjet_Vector);
            }
        }
        bjet_Vectors.push_back(event_bJetVectors);
        

        for (Int_t j = 0; j < nEl; ++j) {
            TLorentzVector electronVector;
            electronVector.SetPtEtaPhiE(el_pt->at(j)/1000.0, el_eta->at(j), el_phi->at(j),  el_e->at(j)/1000.0);
            event_elVectors.push_back(electronVector);
        }
        el_Vectors.push_back(event_elVectors);

        for (Int_t j = 0; j < nMu; ++j) {
            TLorentzVector muonVector;
            muonVector.SetPtEtaPhiE(mu_pt->at(j)/1000.0, mu_eta->at(j), mu_phi->at(j), mu_e->at(j)/1000.0);
            event_muVectors.push_back(muonVector);
        }
        mu_Vectors.push_back(event_muVectors);

          // Create TLorentzVector objects for photons and store them in the vector
        for (Int_t j = 0; j < nPh; ++j) {
            TLorentzVector photonVector;
            photonVector.SetPtEtaPhiE(ph_pt->at(j)/1000.0, ph_eta->at(j), ph_phi->at(j), ph_e->at(j)/1000.0);
            event_phVectors.push_back(photonVector);
        }
        photon_Vectors.push_back(event_phVectors);
        
        // Create TLorentzVector objects for MET and store them in the vector
      for (Int_t j = 0; j < 1; ++j) {
          TLorentzVector METVector;
          METVector.SetPtEtaPhiM(MET_pt/1000.0, 0, MET_phi, 0);
          event_METVectors.push_back(METVector);
      }
      MET_Vectors.push_back(event_METVectors);

      if (nominalWeight) {
         for (Int_t j = 0; j < 1; ++j) {
             event_WeightVectors.push_back(nominalWeight);
         }
      }
      if (!nominalWeight){
         for (Int_t j = 0; j < 1; ++j) {
             event_WeightVectors.push_back(1.0);
         }
      }
      Weight_Vectors.push_back(event_WeightVectors); 

    }
    // Close the file
    inputFile.Close();
    
    TFile outputFile(outputFileName, "RECREATE");
    TTree ntuple("Ntuple", "Ntuple containing TLorentzVector objects");

    // Create TClonesArrays to store TLorentzVectors for each object type
    TClonesArray jetArray("TLorentzVector", 20);
    TClonesArray bjetArray("TLorentzVector", 20);
    TClonesArray electronArray("TLorentzVector", 20);
    TClonesArray muonArray("TLorentzVector", 20);
    TClonesArray photonArray("TLorentzVector", 20);
    TClonesArray METArray("TLorentzVector", 20);
    std::vector<Double_t> *WeightArray = new std::vector<Double_t>();

    // Add branches to the TTree for each TClonesArray
    ntuple.Branch("Jets", &jetArray, 256000, 0);
    ntuple.Branch("BJets", &bjetArray, 256000, 0);
    ntuple.Branch("Electrons", &electronArray, 256000, 0);
    ntuple.Branch("Muons", &muonArray, 256000, 0);
    ntuple.Branch("Photons", &photonArray, 256000, 0);
    ntuple.Branch("MET", &METArray, 256000, 0);
    ntuple.Branch("Weight", "std::vector<Double_t>", &WeightArray, 256000, 0);

    ntuple.Branch("CMS_Energy", &cmsEnergy, "CMS_Energy/D");

    // Loop over all events and fill TLorentzVectors into TClonesArrays
    for (size_t i = 0; i < jet_Vectors.size(); ++i) {
        // Clear the TClonesArrays before filling them for the current event
        jetArray.Clear();
        bjetArray.Clear();
        electronArray.Clear();
        muonArray.Clear();
        photonArray.Clear();
        METArray.Clear();
        WeightArray->clear();

        // Fill TClonesArrays with TLorentzVectors for jets
        for (const auto& jet : jet_Vectors[i]) {
            new (jetArray[jetArray.GetEntries()]) TLorentzVector(jet);
        }

        // Fill TClonesArrays with TLorentzVectors for b-jets
        for (const auto& bjet : bjet_Vectors[i]) {
            new (bjetArray[bjetArray.GetEntries()]) TLorentzVector(bjet);
        }

        // Fill TClonesArrays with TLorentzVectors for electrons
        for (const auto& electron : el_Vectors[i]) {
            new (electronArray[electronArray.GetEntries()]) TLorentzVector(electron);
        }

        // Fill TClonesArrays with TLorentzVectors for muons
        for (const auto& muon : mu_Vectors[i]) {
            new (muonArray[muonArray.GetEntries()]) TLorentzVector(muon);
        }

        // Fill TClonesArrays with TLorentzVectors for photons
        for (const auto& photon : photon_Vectors[i]) {
            new (photonArray[photonArray.GetEntries()]) TLorentzVector(photon);
        }
        
        // Fill TClonesArrays with TLorentzVectors for METs
        for (const auto& MET : MET_Vectors[i]) {
            new (METArray[METArray.GetEntries()]) TLorentzVector(MET);
        }

        // Fill WeightArray with values from Weight_Vectors[i]
        const std::vector<Double_t>& event_W = Weight_Vectors[i];
        for (size_t j = 0; j < event_W.size(); ++j) {
            WeightArray->push_back(event_W[j]);
        }
        
        // Fill the TTree for the current event
        ntuple.Fill();
    }

    // Write the TTree to the output file
    ntuple.Write();

    // Close the output file
    outputFile.Close();
    
    // Read the vectors
    cout << "jet_Vectors.size :" << jet_Vectors.size() << endl;
    cout << "bjet_Vectors.size :" << bjet_Vectors.size() << endl;
    cout << "el_Vectors.size :" << el_Vectors.size() << endl;
    cout << "mu_Vectors.size :" << mu_Vectors.size() << endl;
    cout << "photon_Vectors.size :" << photon_Vectors.size() << endl;
    cout << "MET_Vectors.size :" << MET_Vectors.size() << endl;
    cout << "Weight_Vectors.size :" << Weight_Vectors.size() << endl;
    
/*
    //Read vectors
    // Loop over all events in jet_Vectors
    for (size_t i = 0; i < jet_Vectors.size(); ++i) {
        // Access the vector of TLorentzVectors for the current event
        const std::vector<TLorentzVector>& event_Jets = jet_Vectors[i];
        const std::vector<TLorentzVector>& event_bJets = bjet_Vectors[i];
        const std::vector<TLorentzVector>& event_el = el_Vectors[i];
        const std::vector<TLorentzVector>& event_mu = mu_Vectors[i];
        const std::vector<TLorentzVector>& event_ph = photon_Vectors[i];
        const std::vector<TLorentzVector>& event_MET = MET_Vectors[i];

        // Print the number of jets in the current event
        //cout << "For event# " << i << " #jets: " << event_Jets.size() << " #bjets: " << event_bJets.size() << " #electrons: " << event_el.size() << " #muons: " << event_mu.size() << " #photons: " << event_ph.size() << endl;
     
        // Loop over all jets in the current event
        for (size_t j = 0; j < event_MET.size(); ++j) {
            // Access the j-th TLorentzVector in the current event
            const TLorentzVector& MET = event_MET[j];
         
            // Print the values of the j-th jet in the current event
            Float_t pt = MET.Pt();
            Float_t eta = MET.Eta();
            Float_t phi = MET.Phi();
            Float_t Energy = MET.M();
         
            cout << "Event " << i << ", MET " << j << ": PT=" << pt << ", Eta=" << eta << ", Phi=" << phi << ", Energy=" << Energy << endl;
     }
        // Check if there are any jets in the current event
        if (!event_Jets.empty()) {
            // Loop over all jets in the current event
            for (size_t j = 0; j < event_Jets.size(); ++j) {
                // Access the j-th TLorentzVector in the current event
                const TLorentzVector& jet = event_Jets[j];
                
                // Print the values of the j-th jet in the current event
                Float_t pt = jet.Pt();
                Float_t eta = jet.Eta();
                Float_t phi = jet.Phi();
                Float_t mass = jet.M();
                
                //cout << "Event " << i << ", Jet " << j << ": PT=" << pt << ", Eta=" << eta << ", Phi=" << phi << ", Mass=" << mass << endl;
            }
        } else {
            cout << "Event " << i << ": No jets." << endl;
        }
    }
    */
}


void skim_MasterShef(const char* inputFileName = nullptr, const char* outputFileName = nullptr, Double32_t cmsEnergy = 0.0) {
    // Check if both input and output file names are provided
    if (!inputFileName || !outputFileName) {
        std::cerr << "Usage: root -b -q skim_MasterShef.C(\"inputFileName\", \"outputFileName\", cmsEnergy)" << std::endl;
        return;
    }
    // Call the skim_MasterShef function with the provided file names
    readEventData(inputFileName, outputFileName, cmsEnergy);
}



