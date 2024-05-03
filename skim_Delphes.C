/***************************************************************************
 *  Read Delphes file and save a skimmed version of root file with TLorentzvectors 
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
    TTree* tree = dynamic_cast<TTree*>(inputFile.Get("Delphes"));
    if (!tree) {
        std::cerr << "Error: Unable to access tree Delphes in file " << inputFileName << std::endl;
        return;
    }
    
    // Find values for kMaxJet,kMaxElectron, kMaxMuon, kMaxPhoton
    Int_t kMaxJet = 0;
    Int_t kMaxElectron = 0;
    Int_t kMaxPhoton = 0;
    Int_t kMaxMuon = 0;
    Int_t kMaxMissingET = 0;

    TBranch* jet_Pt_Branch = tree->GetBranch("Jet.PT");
    TBranch* el_Pt_Branch = tree->GetBranch("Electron.PT");
    TBranch* mu_Pt_Branch = tree->GetBranch("Muon.PT");
    TBranch* ph_Pt_Branch = tree->GetBranch("Photon.PT");
    TBranch* MET_Branch = tree->GetBranch("MissingET.MET");

    TLeaf* jet_leaf = jet_Pt_Branch->GetLeaf("Jet.PT");
    TLeaf* el_leaf = el_Pt_Branch->GetLeaf("Electron.PT");
    TLeaf* mu_leaf = mu_Pt_Branch->GetLeaf("Muon.PT");
    TLeaf* ph_leaf = ph_Pt_Branch->GetLeaf("Photon.PT");
    TLeaf* MET_leaf = MET_Branch->GetLeaf("MissingET.MET");

    Long64_t numEntries1 = tree->GetEntries();
    for (Long64_t i = 0; i < numEntries1; ++i) {
        jet_Pt_Branch->GetEntry(i);
        Int_t len_jet = jet_leaf->GetLen();
        if (len_jet > kMaxJet) {
            kMaxJet = len_jet;
        }
    }
    
    for (Long64_t i = 0; i < numEntries1; ++i) {
        el_Pt_Branch->GetEntry(i);
        Int_t len_el = el_leaf->GetLen();
        if (len_el > kMaxElectron) {
            kMaxElectron = len_el;
        }
    }
      
    for (Long64_t i = 0; i < numEntries1; ++i) {
        mu_Pt_Branch->GetEntry(i);
        Int_t len_mu = mu_leaf->GetLen();
        if (len_mu > kMaxMuon) {
            kMaxMuon = len_mu;
        }
    }
    
    for (Long64_t i = 0; i < numEntries1; ++i) {
        ph_Pt_Branch->GetEntry(i);
        Int_t len_ph = ph_leaf->GetLen();
        if (len_ph > kMaxPhoton) {
            kMaxPhoton = len_ph;
        }
    }

    for (Long64_t i = 0; i < numEntries1; ++i) {
        MET_Branch->GetEntry(i);
        Int_t len_MET = MET_leaf->GetLen();
        if (len_MET > kMaxMissingET) {
            kMaxMissingET = len_MET;
        }
    }

    std::cout << "kMaxJet: " << kMaxJet << std::endl;
    std::cout << "kMaxElectron: " << kMaxElectron << std::endl;
    std::cout << "kMaxMuon: " << kMaxMuon << std::endl;
    std::cout << "kMaxPhoton: " << kMaxPhoton << std::endl;
    std::cout << "kMaxMissingET: " << kMaxMissingET << std::endl;

    
    Float_t         Jet_PT[kMaxJet];   //[Jet_]
    Float_t         Jet_Eta[kMaxJet];   //[Jet_]
    Float_t         Jet_Phi[kMaxJet];   //[Jet_]
    Float_t         Jet_Mass[kMaxJet];   //[Jet_]
    UInt_t          Jet_BTag[kMaxJet];   //[Jet_]

    Float_t         Electron_PT[kMaxElectron];   //[Electron_]
    Float_t         Electron_Eta[kMaxElectron];   //[Electron_]
    Float_t         Electron_Phi[kMaxElectron];   //[Electron_]
    //Float_t         Electron_Mass[kMaxElectron];   //[Electron_]

    Float_t         Muon_PT[kMaxMuon];   //[Muon_]
    Float_t         Muon_Eta[kMaxMuon];   //[Muon_]
    Float_t         Muon_Phi[kMaxMuon];   //[Muon_]
    //Float_t         Muon_Mass[kMaxMuon];   //[Muon_]

    Float_t         Photon_PT[kMaxPhoton];   //[Photon_]
    Float_t         Photon_Eta[kMaxPhoton];   //[Photon_]
    Float_t         Photon_Phi[kMaxPhoton];   //[Photon_]
    Float_t         Photon_E[kMaxPhoton];   //[Photon_]

    Float_t         MissingET_MET[kMaxMissingET];   //[MissingET_]
    Float_t         MissingET_Eta[kMaxMissingET];   //[MissingET_]
    Float_t         MissingET_Phi[kMaxMissingET];   //[MissingET_]    

    TBranch        *b_Jet_PT;   //!
    TBranch        *b_Jet_Eta;   //!
    TBranch        *b_Jet_Phi;   //!
    TBranch        *b_Jet_T;   //!
    TBranch        *b_Jet_Mass;   //!
    TBranch        *b_Jet_BTag;   //!
    TBranch        *b_Electron_PT;   //!
    TBranch        *b_Electron_Eta;   //!
    TBranch        *b_Electron_Phi;   //!
    //TBranch        *b_Electron_Mass;   //!
    TBranch        *b_Muon_PT;   //!
    TBranch        *b_Muon_Eta;   //!
    TBranch        *b_Muon_Phi;   //!
    //TBranch        *b_Muon_Mass;   //!
    TBranch        *b_Photon_PT;   //!
    TBranch        *b_Photon_Eta;   //!
    TBranch        *b_Photon_Phi;   //!
    TBranch        *b_Photon_E;   //!

    TBranch        *b_MissingET_MET;   //!
    TBranch        *b_MissingET_Eta;   //!
    TBranch        *b_MissingET_Phi;   //!
    
    tree->SetBranchAddress("Jet.PT", Jet_PT, &b_Jet_PT);
    tree->SetBranchAddress("Jet.Eta", Jet_Eta, &b_Jet_Eta);
    tree->SetBranchAddress("Jet.Phi", Jet_Phi, &b_Jet_Phi);
    tree->SetBranchAddress("Jet.Mass", Jet_Mass, &b_Jet_Mass);
    tree->SetBranchAddress("Jet.BTag", Jet_BTag, &b_Jet_BTag);
    
    tree->SetBranchAddress("Electron.PT", Electron_PT, &b_Electron_PT);
    tree->SetBranchAddress("Electron.Eta", Electron_Eta, &b_Electron_Eta);
    tree->SetBranchAddress("Electron.Phi", Electron_Phi, &b_Electron_Phi);
    //tree->SetBranchAddress("Electron.Mass", Electron_Mass, &b_Electron_Mass);

    tree->SetBranchAddress("Muon.PT", Muon_PT, &b_Muon_PT);
    tree->SetBranchAddress("Muon.Eta", Muon_Eta, &b_Muon_Eta);
    tree->SetBranchAddress("Muon.Phi", Muon_Phi, &b_Muon_Phi);
    //tree->SetBranchAddress("Muon.Mass", Muon_Mass, &b_Muon_Mass);

    tree->SetBranchAddress("Photon.PT", Photon_PT, &b_Photon_PT);
    tree->SetBranchAddress("Photon.Eta", Photon_Eta, &b_Photon_Eta);
    tree->SetBranchAddress("Photon.Phi", Photon_Phi, &b_Photon_Phi);
    tree->SetBranchAddress("Photon.E", Photon_E, &b_Photon_E);

    tree->SetBranchAddress("MissingET.Eta", MissingET_Eta, &b_MissingET_Eta);
    tree->SetBranchAddress("MissingET.Phi", MissingET_Phi, &b_MissingET_Phi);
    tree->SetBranchAddress("MissingET.MET", MissingET_MET, &b_MissingET_MET);
    
    // Loop over entries in the TTree and create Vectors
    std::vector<std::vector<TLorentzVector>> jet_Vectors;
    std::vector<std::vector<TLorentzVector>> bjet_Vectors;
    std::vector<std::vector<TLorentzVector>> el_Vectors;
    std::vector<std::vector<TLorentzVector>> mu_Vectors;
    std::vector<std::vector<TLorentzVector>> photon_Vectors;
    std::vector<std::vector<TLorentzVector>> MET_Vectors;

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

        // Get entry i
        tree->GetEntry(i);
        //cout << "Event " << i << endl; 
        if(i%1000==0) cout << "Event #" <<  i << endl;

        for (Int_t j = 0; j < kMaxJet; ++j) {
            TLorentzVector jet_Vector;
            jet_Vector.SetPtEtaPhiM(Jet_PT[j], Jet_Eta[j], Jet_Phi[j], Jet_Mass[j]);
	    event_JetVectors.push_back(jet_Vector);
        }
        jet_Vectors.push_back(event_JetVectors);

        
        for (Int_t j = 0; j < kMaxJet; ++j) {
            TLorentzVector bjet_Vector;
            //cout << j << " isbTag :" << Jet_BTag[j] << endl;
            if (Jet_BTag[j]==1){
                bjet_Vector.SetPtEtaPhiM(Jet_PT[j], Jet_Eta[j], Jet_Phi[j], Jet_Mass[j]);
                event_bJetVectors.push_back(bjet_Vector);
            }
        }
        bjet_Vectors.push_back(event_bJetVectors);

        for (Int_t j = 0; j < kMaxElectron; ++j) {
            TLorentzVector electronVector;
            electronVector.SetPtEtaPhiM(Electron_PT[j], Electron_Eta[j], Electron_Phi[j],  0.000511 );
            event_elVectors.push_back(electronVector);
        }
        el_Vectors.push_back(event_elVectors);

        for (Int_t j = 0; j < kMaxMuon; ++j) {
            TLorentzVector muonVector;
            muonVector.SetPtEtaPhiM(Muon_PT[j], Muon_Eta[j], Muon_Phi[j], 0.1057 );
            event_muVectors.push_back(muonVector);
        }
        mu_Vectors.push_back(event_muVectors);


          // Create TLorentzVector objects for photons and store them in the vector
        for (Int_t j = 0; j < kMaxPhoton; ++j) {
            TLorentzVector photonVector;
            photonVector.SetPtEtaPhiM(Photon_PT[j], Photon_Eta[j], Photon_Phi[j], 0.0);
            event_phVectors.push_back(photonVector);
        }
        photon_Vectors.push_back(event_phVectors);

        // Create TLorentzVector objects for MET and store them in the vector
        for (Int_t j = 0; j < kMaxMissingET; ++j) {
            TLorentzVector METVector;
            METVector.SetPtEtaPhiM(MissingET_MET[j], 0.0, MissingET_Phi[j], 0.0);
            event_METVectors.push_back(METVector);
        }
        MET_Vectors.push_back(event_METVectors);

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

    // Add branches to the TTree for each TClonesArray
    ntuple.Branch("Jets", &jetArray, 256000, 0);
    ntuple.Branch("BJets", &bjetArray, 256000, 0);
    ntuple.Branch("Electrons", &electronArray, 256000, 0);
    ntuple.Branch("Muons", &muonArray, 256000, 0);
    ntuple.Branch("Photons", &photonArray, 256000, 0);
    ntuple.Branch("MET", &METArray, 256000, 0);
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

        // Print the number of jets in the current event
        cout << "For event# " << i << " #jets: " << event_Jets.size() << " #bjets: " << event_bJets.size() << " #electrons: " << event_el.size() << " #muons: " << event_mu.size() << " #photons: " << event_ph.size() << endl;
        
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

void skim_Delphes(const char* inputFileName = nullptr, const char* outputFileName = nullptr, Double32_t cmsEnergy = 0.0) {
    // Check if both input and output file names are provided
    if (!inputFileName || !outputFileName) {
        std::cerr << "Usage: root -b -q skim_Delphes.C(\"inputFileName\", \"outputFileName\", cmsEnergy)" << std::endl;
        return;
    }
    // Call the skim_Delphes_file function with the provided file names
    readEventData(inputFileName, outputFileName, cmsEnergy);
}


