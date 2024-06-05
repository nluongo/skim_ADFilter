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
    TTree* tree = dynamic_cast<TTree*>(inputFile.Get("Nominal"));
    if (!tree) {
        std::cerr << "Error: Unable to access tree Delphes in file " << inputFileName << std::endl;
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

    Long64_t numEntries = tree->GetEntries();
    cout << "no. of events :" << numEntries << endl;

    std::vector<Double32_t> bJET_pt, bJET_eta, bJET_phi, bJET_mass;
    Int_t  N_bJET;

    std::vector<Double32_t> JET_pt, JET_eta, JET_phi, JET_mass;
    Int_t N_JET;

    std::vector<Double32_t> EL_pt, EL_eta, EL_phi;
    Int_t N_EL;

    std::vector<Double32_t> MU_pt, MU_eta, MU_phi;
    Int_t N_MU;

    std::vector<Double32_t> PH_pt, PH_eta, PH_phi, PH_e;
    Int_t N_PH;

    std::vector<Double32_t> MET_eta, MET_met, met_phi;

    std::vector<Double32_t> Evt_Weight;

    TFile outputFile(outputFileName, "RECREATE");
    TTree ntuple("Ntuple", "Ntuple containing TLorentzVector objects");

    // jet branches
    ntuple.Branch("JET_n",  &N_JET);
    ntuple.Branch("JET_pt", &JET_pt, 256000, 0);
    ntuple.Branch("JET_eta", &JET_eta, 256000, 0);
    ntuple.Branch("JET_phi", &JET_phi, 256000, 0);
    ntuple.Branch("JET_mass", &JET_mass, 256000, 0);

    // b jet branches
    ntuple.Branch("bJET_n",  &N_bJET);
    ntuple.Branch("bJET_pt", &bJET_pt, 256000, 0);
    ntuple.Branch("bJET_eta", &bJET_eta, 256000, 0);
    ntuple.Branch("bJET_phi", &bJET_phi, 256000, 0);
    ntuple.Branch("bJET_mass", &bJET_mass, 256000, 0);
    
    // Electron branches
    ntuple.Branch("EL_n", &N_EL);
    ntuple.Branch("EL_pt", &EL_pt, 256000, 0);
    ntuple.Branch("EL_eta", &EL_eta, 256000, 0);
    ntuple.Branch("EL_phi", &EL_phi, 256000, 0);

    // Muon branches
    ntuple.Branch("MU_n", &N_MU);
    ntuple.Branch("MU_pt", &MU_pt, 256000, 0);
    ntuple.Branch("MU_eta", &MU_eta, 256000, 0);
    ntuple.Branch("MU_phi", &MU_phi, 256000, 0);

    // Photon branches
    ntuple.Branch("PH_n", &N_PH);
    ntuple.Branch("PH_pt", &PH_pt, 256000, 0);
    ntuple.Branch("PH_eta", &PH_eta, 256000, 0);
    ntuple.Branch("PH_phi", &PH_phi, 256000, 0);
    ntuple.Branch("PH_e", &PH_e, 256000, 0);

    // MissingET branches
    ntuple.Branch("MET_eta", &MET_eta, 256000, 0);
    ntuple.Branch("MET_phi", &met_phi, 256000, 0);
    ntuple.Branch("MET_met", &MET_met, 256000, 0);
    
    ntuple.Branch("Evt_Weight", &Evt_Weight, 256000, 0);
    
    TH1F* meta = new TH1F("meta","meta",10,0,10);
    meta->Fill("CMS energy [GeV]",(Double_t)cmsEnergy);
    meta->Write();

    for (Long64_t i = 0; i < numEntries; ++i) {
        // Clear the vector for the current event
        //jet_Vectors.clear();

        JET_pt.clear(); 
        JET_eta.clear();
        JET_phi.clear();
        JET_mass.clear();
        N_JET=0;
        
        bJET_pt.clear();
        bJET_eta.clear();
        bJET_phi.clear();
        bJET_mass.clear();
        N_bJET=0;
        
        EL_pt.clear();
        EL_eta.clear();
        EL_phi.clear();
        N_EL = 0;

        MU_pt.clear();
        MU_eta.clear();
        MU_phi.clear();
        N_MU = 0;

        PH_pt.clear();
        PH_eta.clear();
        PH_phi.clear();
        PH_e.clear();
        N_PH = 0;

        MET_eta.clear();
        met_phi.clear();
        MET_met.clear();
        
        Evt_Weight.clear();

        // Get entry i
        tree->GetEntry(i);
        //cout << "Event " << i << endl; 
        if(i%1000==0) cout << "Event #" <<  i << endl;

        cout << " " << endl;

        N_JET=nJets;
        for (Int_t j = 0; j < nJets; ++j) {
            cout << "jet : " << j << " " << jet_pt->at(j)/1000.0 << endl;
            JET_pt.push_back(jet_pt->at(j)/1000.0);
            JET_eta.push_back(jet_eta->at(j));
            JET_phi.push_back(jet_phi->at(j));
            JET_mass.push_back(jet_e->at(j)/1000.0);
        }
      
        for (Int_t j = 0; j < nJets; ++j) {
            if ((jet_IsBJet->at(j)) != 1)
                continue;
            cout << "bjet : " << j << " " << jet_pt->at(j)/1000.0 << endl;
            bJET_pt.push_back(jet_pt->at(j)/1000.0);
            bJET_eta.push_back(jet_eta->at(j));
            bJET_phi.push_back(jet_phi->at(j));
            bJET_mass.push_back(jet_e->at(j)/1000.0);
        }
        
        // Electron data
        for (Int_t j = 0; j < nEl; ++j) {
            cout << "electron: " << j << " " << el_pt->at(j)/1000.0 << endl;
            EL_pt.push_back(el_pt->at(j)/1000.0);
            EL_eta.push_back(el_eta->at(j));
            EL_phi.push_back(el_phi->at(j));
        }

        // Muon data
        for (Int_t j = 0; j < nMu; ++j) {
            cout << "muon: " << j << " " << mu_pt->at(j)/1000.0 << endl;
            MU_pt.push_back(mu_pt->at(j)/1000.0);
            MU_eta.push_back(mu_eta->at(j));
            MU_phi.push_back(mu_phi->at(j));
        }

        // Photon data
        for (Int_t j = 0; j < nPh; ++j) {
            cout << "photon: " << j << " " << ph_pt->at(j)/1000.0 << endl;
            PH_pt.push_back(ph_pt->at(j)/1000.0);
            PH_eta.push_back(ph_eta->at(j));
            PH_phi.push_back(ph_phi->at(j));
            PH_e.push_back(ph_e->at(j)/1000.0);
        }

        // MissingET data
        for (Int_t j = 0; j < 1; ++j) {
            cout << "missingET pt: " << j << " " << MET_pt/1000.0 << endl;
            cout << "missingET phi: " << j << " " << MET_phi << endl;
            MET_eta.push_back(0);
            met_phi.push_back(MET_phi);
            MET_met.push_back(MET_pt/1000.0);
        }
        
        // Event.Weight data
        for (Int_t j = 0; j < 1; ++j) {
            cout << "Event.Weight: " << j << " " <<  nominalWeight << endl;
            Evt_Weight.push_back(nominalWeight);
        }
      
       // fill ntuples
       ntuple.Fill();

    }

     // Close the file
    inputFile.Close();

   // Write the TTree to the output file
    ntuple.Write();
 
    // Close the output file
    outputFile.Close();
    
    
}

void skim_MasterShef(const char* inputFileName = nullptr, const char* outputFileName = nullptr, Double32_t cmsEnergy = 0.0) {
    // Check if both input and output file names are provided
    if (!inputFileName || !outputFileName) {
        std::cerr << "Usage: root -b -q skim_MasterShef.C(\"inputFileName\", \"outputFileName\", cmsEnergy)" << std::endl;
        return;
    }
    // Call the skim_Delphes_file function with the provided file names
    readEventData(inputFileName, outputFileName, cmsEnergy);
}


