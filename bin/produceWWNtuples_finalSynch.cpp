#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <iomanip>
#include <ctime>
#include <map>
#include <math.h>

#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TChain.h"
#include "TProfile.h"
#include "TMath.h"
#include "TString.h"
#include "TClass.h"
#include "TApplication.h"
#include "TLorentzVector.h"

#include "../interface/setInputTree.h"
#include "../interface/setOutputTreeSynch.h"
#include "../interface/METzCalculator.h"
#include "../interface/METzCalculator_Run2.h"
#include "../interface/analysisUtils.h"

using namespace std;

//*******MAIN*******************************************************************

int main (int argc, char** argv)
{ 
  std::string inputFolder = argv[1];
  std::string outputFile = argv[2];
  bool isMC = argv[3];
  std::string leptonName = argv[4];
  std::string inputTreeName = argv[5];
  std::string inputFile = argv[6];
  std::string xSecWeight = argv[7];
  std::string numberOfEntries = argv[8];
  float weight = std::atof(xSecWeight.c_str())/std::atof(numberOfEntries.c_str());
  if (strcmp(leptonName.c_str(),"el")!=0 && strcmp(leptonName.c_str(),"mu")!=0) {
    std::cout<<"Error: wrong lepton category"<<std::endl;
    return(-1);
  }

  int evento=-17;

  TLorentzVector W,MET,LEP;
  TLorentzVector NU0,NU1,NU2;
  TLorentzVector JET, HADW, AK4;
  TLorentzVector VBF1,VBF2,TOT;
  TLorentzVector ELE,MU;

  std::vector<TLorentzVector> tightMuon;
  std::vector<TLorentzVector> looseMuon;
  std::vector<TLorentzVector> tightEle;
  std::vector<TLorentzVector> looseEle;

  setInputTree *ReducedTree = new setInputTree (inputTreeName.c_str());
  ReducedTree->Init();


  char command1[3000];
  sprintf(command1, "xrd eoscms dirlist %s/%s/  | awk '{print \"root://eoscms.cern.ch/\"$5}' > listTemp_%s.txt", (inputFolder).c_str(), (inputFile).c_str(), outputFile.c_str());
  std::cout<<command1<<std::endl;
  system(command1);
  char list1[2000];
  sprintf (list1, "listTemp_%s.txt", outputFile.c_str());
  ifstream rootList (list1);

  while (!rootList.eof())
    {
      char iRun_tW[700];
      rootList >> iRun_tW;
      ReducedTree->fChain->Add(iRun_tW);
    }

  char command3[300];
  sprintf(command3, "rm listTemp_%s.txt", outputFile.c_str());
  system(command3);

  int cutEff[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  //---------output tree----------------
  TFile* outROOT = TFile::Open((std::string("output/output_synch_")+leptonName+std::string("/")+outputFile+(".root")).c_str(),"recreate");
  outROOT->cd();
  TTree* outTree = new TTree("otree", "otree");
  outTree->SetDirectory(0);
  setOutputTreeSynch *WWTree = new setOutputTreeSynch(outTree);

  //---------start loop on events------------
  for (Long64_t jentry=0; jentry<ReducedTree->fChain->GetEntries();jentry++) {

    Long64_t iEntry = ReducedTree->LoadTree(jentry);
    if (iEntry < 0) break;
    int nb = ReducedTree->fChain->GetEntry(jentry);   
    // if (Cut(ientry) < 0) continue;                                                                                                                           

    tightMuon.clear();
    tightEle.clear();
    looseMuon.clear();
    looseEle.clear();

    if(iEntry % 1000 == 0)    
      cout << "read entry: " << iEntry << endl;

    WWTree->initializeVariables(); //initialize all variables
    
    WWTree->issignal = 0;
    WWTree->wSampleWeight = weight; //xsec/numberOfEntries
    WWTree->totalEventWeight = 1.; //temporary value
    WWTree->eff_and_pu_Weight = 1.; //temporary value
    
    //require at least one lepton and one jet
    //    if ( strcmp(leptonName.c_str(),"el")==0 && ReducedTree->ElectronsNum==0) continue; 
    //    if ( strcmp(leptonName.c_str(),"mu")==0 && ReducedTree->MuonsNum==0) continue;      
        
    //save event variables
    WWTree->run   = ReducedTree->RunNum;
    WWTree->event = ReducedTree->EvtNum;
    WWTree->lumi = ReducedTree->LumiBlockNum;
   // WWTree->njets = ReducedTree->NJets;
    WWTree->nPV  = ReducedTree->NVtx;

        if(WWTree->event==evento)  std::cout<<"debug: "<<std::endl;
    /////////////////THE SELECTED LEPTON
    int nTightLepton=0;
    if (strcmp(leptonName.c_str(),"el")==0) {
      float tempPt=0.;
      for (int i=0; i<ReducedTree->ElectronsNum; i++) {
	if (ReducedTree->Electrons_isHEEP[i]==false) continue;       
        if (ReducedTree->ElectronsPt[i]<=90) continue;
	if (ReducedTree->ElectronsPt[i]<tempPt) continue;
	ELE.SetPtEtaPhiE(ReducedTree->ElectronsPt[i],ReducedTree->ElectronsEta[i],ReducedTree->ElectronsPhi[i],ReducedTree->ElectronsE[i]);
	tightEle.push_back(ELE);
	WWTree->l_pt  = ReducedTree->ElectronsPt[i];
	WWTree->l_eta = ReducedTree->ElectronsEta[i];
	WWTree->l_phi = ReducedTree->ElectronsPhi[i];	
	WWTree->l_e= ReducedTree->ElectronsE[i];	
	tempPt = WWTree->l_pt;
	nTightLepton++;
      }
    }
    else if (strcmp(leptonName.c_str(),"mu")==0) {
      float tempPt=0.;
      for (int i=0; i<ReducedTree->MuonsNum; i++) {
	if (ReducedTree->Muons_isHighPt[i]==false) continue;
        if(WWTree->event==evento)  std::cout<<"debug: "<<i<<std::endl;
        if (ReducedTree->Muons_isPFMuon[i]==false) continue; 
        if(WWTree->event==evento)  std::cout<<"debug: "<<i<<std::endl;
        if ((ReducedTree->Muons_trackIso[i]/ReducedTree->MuonsPt[i])>=0.1) continue;
        if(WWTree->event==evento)  std::cout<<"debug: "<<i<<std::endl;
        if (ReducedTree->MuonsPt[i]<50) continue;
        if(WWTree->event==evento)  std::cout<<"debug: "<<i<<std::endl;
        if (fabs(ReducedTree->MuonsEta[i])>=2.1) continue;
        if(WWTree->event==evento)  std::cout<<"debug: "<<i<<std::endl;
	MU.SetPtEtaPhiE(ReducedTree->MuonsPt[i],ReducedTree->MuonsEta[i],ReducedTree->MuonsPhi[i],ReducedTree->MuonsE[i]);
	tightMuon.push_back(MU);
	if (ReducedTree->MuonsPt[i]<tempPt) continue;
	WWTree->l_pt  = ReducedTree->MuonsPt[i];
	WWTree->l_eta = ReducedTree->MuonsEta[i];
	WWTree->l_phi = ReducedTree->MuonsPhi[i];
	WWTree->l_e = ReducedTree->MuonsE[i];
	tempPt = WWTree->l_pt;
	nTightLepton++;
      }
    }
    if (nTightLepton==0) continue; //no leptons with required ID
        if(WWTree->event==evento)  std::cout<<"debug: "<<std::endl;

    //VETO ADDITIONAL LEPTONS
    int nLooseLepton=0;
    for (int i=0; i<ReducedTree->ElectronsNum; i++) {
      if (ReducedTree->Electrons_isHEEP[i]==false) continue;       
      if (ReducedTree->ElectronsPt[i]<35) continue;       
      ELE.SetPtEtaPhiE(ReducedTree->ElectronsPt[i],ReducedTree->ElectronsEta[i],ReducedTree->ElectronsPhi[i],ReducedTree->ElectronsE[i]);
      looseEle.push_back(ELE);      
      nLooseLepton++;
    }
    WWTree->nLooseEle = nLooseLepton;

    for (int i=0; i<ReducedTree->MuonsNum; i++) {
      if (ReducedTree->Muons_isHighPt[i]==false) continue;
      if ((ReducedTree->Muons_trackIso[i]/ReducedTree->MuonsPt[i])>=0.1) continue;
      if (fabs(ReducedTree->MuonsEta[i])>=2.4) continue;
      if (ReducedTree->MuonsPt[i]<20) continue;
      MU.SetPtEtaPhiE(ReducedTree->MuonsPt[i],ReducedTree->MuonsEta[i],ReducedTree->MuonsPhi[i],ReducedTree->MuonsE[i]);
      looseMuon.push_back(MU);
      nLooseLepton++;
    }
    WWTree->nLooseMu = nLooseLepton-WWTree->nLooseEle;
    
    if (nLooseLepton!=1) continue; //no additional leptons
    cutEff[0]++;
        if(WWTree->event==evento)  std::cout<<"debug: "<<std::endl;

    //preselection on jet pt and met
    if (ReducedTree->METPt < 40) continue; 
    cutEff[1]++;
        if(WWTree->event==evento)  std::cout<<"debug: "<<std::endl;

    MET.SetPtEtaPhiE(ReducedTree->METPt,0.,ReducedTree->METPhi,0.);
    LEP.SetPtEtaPhiE(WWTree->l_pt,WWTree->l_eta,WWTree->l_phi,WWTree->l_e);

    //    if (ReducedTree->AK8Jets_PtCorr[0] < 150) continue;     
    
    //lepton Pt preselection
    //    if ( strcmp(leptonName.c_str(),"el")==0 && ReducedTree->ElectronsPt[0]<30) continue; 
    //    if ( strcmp(leptonName.c_str(),"mu")==0 && ReducedTree->MuonsPt[0]<30) continue; 

    //////////////THE MET

    // Calculate Neutrino Pz using all the possible choices : 
    // type0 -> if real roots, pick the one nearest to the lepton Pz except when the Pz so chosen
    //               is greater than 300 GeV in which case pick the most central root.
    // type1 -> type = 1: if real roots, choose the one closest to the lepton Pz if complex roots, use only the real part.
    //          type = 2: if real roots, choose the most central solution. if complex roots, use only the real part. 
    //          type = 3: if real roots, pick the largest value of the cosine*

    float Wmass = 80.385;

    TLorentzVector W_mu, W_Met;

    W_mu.SetPtEtaPhiE(WWTree->l_pt,WWTree->l_eta,WWTree->l_phi,WWTree->l_e);
    W_Met.SetPxPyPzE(ReducedTree->METPt * TMath::Cos(ReducedTree->METPhi), ReducedTree->METPt * TMath::Sin(ReducedTree->METPhi), 0., sqrt(ReducedTree->METPt*ReducedTree->METPt));

    if(W_mu.Pt()<=0 || W_Met.Pt() <= 0 ){ std::cerr<<" Negative Lepton - Neutrino Pt "<<std::endl; continue ; }

    // type0 calculation of neutrino pZ
    METzCalculator NeutrinoPz_type0;
    METzCalculator_Run2 NeutrinoPz_run2;
    NeutrinoPz_type0.SetMET(W_Met);
    NeutrinoPz_type0.SetLepton(W_mu);
    NeutrinoPz_type0.SetLeptonType(leptonName.c_str());

    NeutrinoPz_run2.SetMET(W_Met);
    NeutrinoPz_run2.SetLepton(W_mu);
    NeutrinoPz_run2.SetLeptonType(leptonName.c_str());

    double pz1_type0 = NeutrinoPz_type0.Calculate(); // Default one -> according to type0
    double pz2_type0 = NeutrinoPz_type0.getOther(); // Default one

    double pz1_run2 = NeutrinoPz_run2.Calculate(); 

    // don't touch the neutrino pT
    TLorentzVector W_neutrino_type0_met; 
    W_neutrino_type0_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));

    //    W_mass_type0_met = (W_neutrino_type0_met+W_mu).M();
    //    W_pz_type0_met = (W_neutrino_type0_met+W_mu).Pz();
    //    W_nu1_pz_type0_met = pz1_type0;
    //    W_nu2_pz_type0_met = pz2_type0;

    // chenge the neutrino pT in case of complex solution in order to make it real
    TLorentzVector W_neutrino_type0; 
    W_neutrino_type0.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));

    if (NeutrinoPz_type0.IsComplex()) {// if this is a complix, change MET
      double nu_pt1 = NeutrinoPz_type0.getPtneutrino(1);
      double nu_pt2 = NeutrinoPz_type0.getPtneutrino(2);
      TLorentzVector W_neutrino_1;
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(ReducedTree->METPhi),
			      nu_pt1 * TMath::Sin(ReducedTree->METPhi), pz1_type0, sqrt(nu_pt1*nu_pt1 + pz1_type0*pz1_type0) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(ReducedTree->METPhi),
			      nu_pt2 * TMath::Sin(ReducedTree->METPhi), pz1_type0, sqrt(nu_pt2*nu_pt2 + pz1_type0*pz1_type0) );

      if ( fabs((W_mu+W_neutrino_1).M()-Wmass) < fabs((W_mu+W_neutrino_2).M()-Wmass) ) W_neutrino_type0 = W_neutrino_1;
      else W_neutrino_type0 = W_neutrino_2;
    }

    //    W_mass_type0 = (W_mu+W_neutrino_type0).M();
    //    W_pz_type0 = (W_mu+W_neutrino_type0).Pz();
    //    W_nu1_pz_type0 = pz1_type0;
    //    W_nu2_pz_type0 = pz2_type0;

    // type2 calculation of neutrino pZ
    METzCalculator NeutrinoPz_type2;
    NeutrinoPz_type2.SetMET(W_Met);
    NeutrinoPz_type2.SetLepton(W_mu);
    NeutrinoPz_type2.SetLeptonType(leptonName.c_str());

    double pz1_type2 = NeutrinoPz_type2.Calculate(2); // Default one -> according to type2
    double pz2_type2 = NeutrinoPz_type2.getOther(); // Default one

    // don't touch the neutrino pT
    TLorentzVector W_neutrino_type2_met; 
    W_neutrino_type2_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));
    //    W_mass_type2_met = (W_neutrino_type2_met+W_mu).M();
    //    W_pz_type2_met = (W_neutrino_type2_met+W_mu).Pz();
    //    W_nu1_pz_type2_met = pz1_type2;
    //    W_nu2_pz_type2_met = pz2_type2;

    // chenge the neutrino pT in case of complex solution in order to make it real
    TLorentzVector W_neutrino_type2; 
    W_neutrino_type2.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));

    if (NeutrinoPz_type2.IsComplex()) {// if this is a complix, change MET
      double nu_pt1 = NeutrinoPz_type2.getPtneutrino(1);
      double nu_pt2 = NeutrinoPz_type2.getPtneutrino(2);
      TLorentzVector W_neutrino_1;
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(ReducedTree->METPhi),
			      nu_pt1 * TMath::Sin(ReducedTree->METPhi), pz1_type2, sqrt(nu_pt1*nu_pt1 + pz1_type2*pz1_type2) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(ReducedTree->METPhi),
			      nu_pt2 * TMath::Sin(ReducedTree->METPhi), pz1_type2, sqrt(nu_pt2*nu_pt2 + pz1_type2*pz1_type2) );

      if ( fabs((W_mu+W_neutrino_1).M()-Wmass) < fabs((W_mu+W_neutrino_2).M()-Wmass) ) W_neutrino_type2 = W_neutrino_1;
      else W_neutrino_type2 = W_neutrino_2;
    }

    //    W_mass_type2 = (W_mu+W_neutrino_type2).M();
    //    W_pz_type2 = (W_mu+W_neutrino_type2).Pz();
    //    W_nu1_pz_type2 = pz1_type2;
    //    W_nu2_pz_type2 = pz2_type2;

    WWTree->pfMET   = sqrt(ReducedTree->METPt*ReducedTree->METPt);
    WWTree->pfMETPhi = ReducedTree->METPhi;
    WWTree->nu_pz_type0 = pz1_type0;
    WWTree->nu_pz_type2 = pz1_type2;
    WWTree->nu_pz_run2 = pz1_run2;
    WWTree->nu_pz_isre = 1-NeutrinoPz_run2.IsComplex();
    WWTree->nu_pz_run2_oth = NeutrinoPz_run2.getOther();
    WWTree->nu_pz_run2_type = NeutrinoPz_run2.getType();


    /////////////////THE LEPTONIC W
    
    LEP.SetPtEtaPhiE(WWTree->l_pt,WWTree->l_eta,WWTree->l_phi,WWTree->l_e);
    NU0.SetPxPyPzE(ReducedTree->METPt*TMath::Cos(ReducedTree->METPhi),ReducedTree->METPt*TMath::Sin(ReducedTree->METPhi),WWTree->nu_pz_type0,TMath::Sqrt(WWTree->pfMET*WWTree->pfMET+WWTree->nu_pz_type0*WWTree->nu_pz_type0));
    NU2.SetPxPyPzE(ReducedTree->METPt*TMath::Cos(ReducedTree->METPhi),ReducedTree->METPt*TMath::Sin(ReducedTree->METPhi),WWTree->nu_pz_type2,TMath::Sqrt(WWTree->pfMET*WWTree->pfMET+WWTree->nu_pz_type2*WWTree->nu_pz_type2));
    NU1.SetPxPyPzE(ReducedTree->METPt*TMath::Cos(ReducedTree->METPhi),ReducedTree->METPt*TMath::Sin(ReducedTree->METPhi),WWTree->nu_pz_run2,TMath::Sqrt(WWTree->pfMET*WWTree->pfMET+WWTree->nu_pz_run2*WWTree->nu_pz_run2));
    W = LEP + NU2;
    
    WWTree->W_pt = W.Pt();
    WWTree->W_eta = W.Eta();
    WWTree->W_phi = W.Phi();
    WWTree->W_mt = TMath::Sqrt(2*LEP.Et()*NU2.Et()*(1-TMath::Cos(LEP.DeltaPhi(NU2))));
    //    W_mt = W.Mt();


    //FOR THE SYNCHORNIZATION!!! REMOVE IT FOR THE REAL ANALYSIS!!!!
    //    NU2.SetPtEtaPhiE(ReducedTree->METPt,0.,ReducedTree->METPhi,0.);
    //    W = NU2+LEP; 
    ////
        if(WWTree->event==evento)  std::cout<<"debug: "<<std::endl;

    if (W.Pt()<200) continue;
    cutEff[2]++;
        if(WWTree->event==evento)  std::cout<<"debug: "<<std::endl;


    //    if (WWTree->W_pt < 150) continue;
//    if (WWTree->deltaR_lak8jet < (TMath::Pi()/2.0))   continue;

    ///////////THE FAT JET
    float tempPt=0.;
    int nGoodAK8jets=0;
    if (ReducedTree->AK8JetsNum < 1 ) continue; 
        if(WWTree->event==evento)  std::cout<<"debug: "<<std::endl;
    for (unsigned int i=0; i<ReducedTree->AK8JetsNum; i++)
      {
	bool isCleanedJet = true;
	if (ReducedTree->AK8Jets_PtCorr[i]<200 || fabs(ReducedTree->AK8JetsEta[i])>2.4)  continue; //be careful: this is not inside the synchntuple code
	if (ReducedTree->AK8Jets_PtCorr[i]<=tempPt) continue; //save the jet with the largest pt
	if (ReducedTree->AK8Jets_AK8isLooseJetId[i]==false) continue; //fat jet must satisfy loose ID

	//CLEANING FROM LEPTONS
	for (int j=0; j<tightEle.size(); j++) {
	  if (deltaR(tightEle.at(j).Eta(), tightEle.at(j).Phi(),
		     ReducedTree->AK8JetsEta[i],   ReducedTree->AK8JetsPhi[i]) <1.0)
	    isCleanedJet = false;
	}
	for (int j=0; j<tightMuon.size(); j++) {
	  if (deltaR(tightMuon.at(j).Eta(), tightMuon.at(j).Phi(),
		     ReducedTree->AK8JetsEta[i],   ReducedTree->AK8JetsPhi[i]) <1.0)
	    isCleanedJet = false;
	}

	/*	for (int j=0; j<ReducedTree->ElectronsNum; j++) {
	  if (ReducedTree->Electrons_isHEEP[j]==false) continue;     
          if (ReducedTree->ElectronsPt[j]<=90) continue;  
	  if (deltaR(ReducedTree->ElectronsEta[j], ReducedTree->ElectronsPhi[j],
		     ReducedTree->AK8JetsEta[i],   ReducedTree->AK8JetsPhi[i]) <1.0)
	    isCleanedJet = false;
	}
	for (int j=0; j<ReducedTree->MuonsNum; j++) {
	  if (ReducedTree->Muons_isHighPt[i]==false) continue;
	  if ((ReducedTree->Muons_trackIso[i]/ReducedTree->MuonsPt[i])>=0.1) continue;
	  if (ReducedTree->MuonsPt[i]<50) continue;
	  if (fabs(ReducedTree->MuonsEta[i])>=2.1) continue;
	  if (deltaR(ReducedTree->MuonsEta[j], ReducedTree->MuonsPhi[j],
		     ReducedTree->AK8JetsEta[i],   ReducedTree->AK8JetsPhi[i]) <1.0)
	    isCleanedJet = false;
	}
	*/

	if (isCleanedJet==false) continue; //jet is overlapped with a lepton

	WWTree->jet_pt  = ReducedTree->AK8Jets_PtCorr[i];
	WWTree->jet_eta = ReducedTree->AK8JetsEta[i];
	WWTree->jet_phi = ReducedTree->AK8JetsPhi[i];
	WWTree->jet_e   = ReducedTree->AK8Jets_ECorr[i];
	WWTree->jet_mass_pruned   = ReducedTree->AK8Jets_prunedMass[i];
        WWTree->jet_mass_softdrop   = ReducedTree->AK8Jets_softDropMass[i];
	WWTree->jet_mass_tr   = ReducedTree->AK8Jets_trimmedMass[i];
	WWTree->jet_mass_fi   = ReducedTree->AK8Jets_filteredMass[i];
	WWTree->jet_tau2tau1   = ReducedTree->AK8Jets_tau2[i]/ReducedTree->AK8Jets_tau1[i];
	tempPt = WWTree->jet_pt;
	nGoodAK8jets++;
      }


    if (nGoodAK8jets==0) continue; //not found a good hadronic W candidate
    cutEff[3]++;
        if(WWTree->event==evento)  std::cout<<"debug: "<<std::endl;

    if (WWTree->jet_pt<200) continue;
    cutEff[4]++;

        if(WWTree->event==evento)  std::cout<<"debug: "<<std::endl;

    //////////////////ANGULAR VARIABLES
    JET.SetPtEtaPhiE(WWTree->jet_pt,WWTree->jet_eta,WWTree->jet_phi,WWTree->jet_e);
    WWTree->deltaR_lak8jet = JET.DeltaR(LEP);
    WWTree->deltaphi_METak8jet = JET.DeltaPhi(NU2);
    WWTree->deltaphi_Vak8jet = JET.DeltaPhi(W);
    if (WWTree->deltaR_lak8jet>(TMath::Pi()/2.0) && fabs(WWTree->deltaphi_METak8jet)>2.0 && fabs(WWTree->deltaphi_Vak8jet)>2.0)
      WWTree->issignal=1;

    //FOUR-BODY INVARIANT MASS
    WWTree->mass_lvj_type0 = (LEP + NU0 + JET).M();
    WWTree->m_lvj = (LEP + NU2 + JET).M();
    WWTree->mass_lvj_run2  = (LEP + NU1 + JET).M();

    /////////VBF and b-tag section
    bool fillVBF = true;

    HADW.SetPtEtaPhiE(WWTree->jet_pt,WWTree->jet_eta,WWTree->jet_phi,WWTree->jet_e); //AK8 fat jet (hadronic W)
    std::vector<int> indexGoodJets;
    indexGoodJets.clear();
    if (indexGoodJets.size()!=0)  fillVBF=false;

    WWTree->njets=0;
    WWTree->nBTagJet_loose=0;
    WWTree->nbtag=0;
    WWTree->nBTagJet_tight=0;

    float oldDeltaR = 1000.;
    float oldDeltaRLep = 1000.;
    int indexCloserJet = -1;
    int indexCloserJetLep = -1;


    //    LEP.SetPtEtaPhiE(WWTree->l_pt,WWTree->l_eta,WWTree->l_phi,WWTree->l_e);
    //    NU2.SetPxPyPzE(ReducedTree->METPt*TMath::Cos(ReducedTree->METPhi),ReducedTree->METPt*TMath::Sin(ReducedTree->METPhi),WWTree->nu_pz_type2,TMath::Sqrt(WWTree->pfMET*WWTree->pfMET+WWTree->nu_pz_type2*WWTree->nu_pz_type2));
    //    W = LEP + NU2;

    for (unsigned int i=0; i<ReducedTree->JetsNum; i++) //loop on AK4 jet
      {
	bool isCleanedJet = true;
	if (ReducedTree->Jets_PtCorr[i]<=30 || ReducedTree->JetsPt[i]<=20 || fabs(ReducedTree->JetsEta[i])>=2.4)  continue;
	if (ReducedTree->Jets_isLooseJetId[i]==false) continue;

	//CLEANING
	if (deltaR(WWTree->jet_eta, WWTree->jet_phi,
		       ReducedTree->JetsEta[i],   ReducedTree->JetsPhi[i]) <0.8)
	  isCleanedJet = false;

	//CLEANING FROM LEPTONS
	for (int j=0; j<tightEle.size(); j++) {
	  if (deltaR(tightEle.at(j).Eta(), tightEle.at(j).Phi(),
		     ReducedTree->JetsEta[i],   ReducedTree->JetsPhi[i]) <0.3) {
	    isCleanedJet = false;
	  }
	}
	for (int j=0; j<tightMuon.size(); j++) {
	  if (deltaR(tightMuon.at(j).Eta(), tightMuon.at(j).Phi(),
		     ReducedTree->JetsEta[i],   ReducedTree->JetsPhi[i]) <0.3) {
	    isCleanedJet = false;
	  }
	}

	/*
	for (int j=0; j<ReducedTree->ElectronsNum; j++) {
	  if (ReducedTree->Electrons_isHEEP[j]==false) continue;       
          if (ReducedTree->ElectronsPt[j]<=90) continue;
	  if (deltaR(ReducedTree->ElectronsEta[j], ReducedTree->ElectronsPhi[j],
		     ReducedTree->JetsEta[i],   ReducedTree->JetsPhi[i]) <0.3)
	    isCleanedJet = false;
	}      
	for (int j=0; j<ReducedTree->MuonsNum; j++) {
	  if (ReducedTree->Muons_isHighPt[i]==false) continue;
	  if ((ReducedTree->Muons_trackIso[i]/ReducedTree->MuonsPt[i])>=0.1) continue;
	  if (ReducedTree->MuonsPt[i]<50) continue;
	  if (fabs(ReducedTree->MuonsEta[i])>=2.1) continue;
	    isCleanedJet = false;
	}
	*/
	if (isCleanedJet==false) continue;


		if(WWTree->event==evento) { std::cout<<"debug: "<<ReducedTree->Jets_bDiscriminatorICSV[i]<<" "<<ReducedTree->Jets_PtCorr[i]<<" "<<ReducedTree->JetsEta[i]<<std::endl; getchar(); }

	WWTree->njets++;
	//fill B-Tag info
        if (ReducedTree->Jets_bDiscriminatorICSV[i]>0.605)   WWTree->nBTagJet_loose++;
        if (ReducedTree->Jets_bDiscriminatorICSV[i]>0.890)   WWTree->nbtag++;
        if (ReducedTree->Jets_bDiscriminatorICSV[i]>0.970)   WWTree->nBTagJet_tight++;

	AK4.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[i],ReducedTree->JetsEta[i],ReducedTree->JetsPhi[i],ReducedTree->Jets_ECorr[i]);

	float deltaRlep = W.DeltaR(AK4);
	if (deltaRlep<oldDeltaRLep) indexCloserJetLep = i;

	float deltaR = HADW.DeltaR(AK4);
	if (deltaR<0.8) continue; //the vbf jets must be outside the had W cone
	
	if (WWTree->njets!=0) {
	  if (WWTree->jet2_pt!=0) {
	    WWTree->jet3_pt=ReducedTree->Jets_PtCorr[i];
	    WWTree->jet3_btag=ReducedTree->Jets_bDiscriminatorICSV[i];
	  }
	  else {
	    WWTree->jet2_pt=ReducedTree->Jets_PtCorr[i];
	    WWTree->jet2_btag=ReducedTree->Jets_bDiscriminatorICSV[i];
	  }
	}	
	if (deltaR<oldDeltaR)  indexCloserJet = i; //index of the closest jet to the AK8
	indexGoodJets.push_back(i); //save index of the "good" vbf jets candidate
      }
	
    if (indexGoodJets.size()<2)  fillVBF=false; //check if at least 2 jets are inside the collection

    if (indexCloserJet>=0) { //fill hadronic top mass
      AK4.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[indexCloserJet],ReducedTree->JetsEta[indexCloserJet],ReducedTree->JetsPhi[indexCloserJet],ReducedTree->Jets_ECorr[indexCloserJet]);
      WWTree->mass_ungroomedjet_closerjet  = (HADW + AK4).M();
    }
    if (indexCloserJetLep>=0) { //fill leptonic top mass
      AK4.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[indexCloserJetLep],ReducedTree->JetsEta[indexCloserJetLep],ReducedTree->JetsPhi[indexCloserJetLep],ReducedTree->Jets_ECorr[indexCloserJetLep]);
      WWTree->mass_leptonic_closerjet  = (W + AK4).M();
    }

    if (fillVBF) 
      {
	float tempPtMax=0.;
	int nVBF1=-1, nVBF2=-1; //position of the two vbf jets
    
	for (unsigned int i=0; i<indexGoodJets.size()-1; i++) {
	  for (unsigned int ii=i+1; ii<indexGoodJets.size(); ii++) {
	    VBF1.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[indexGoodJets.at(i)],ReducedTree->JetsEta[indexGoodJets.at(i)],ReducedTree->JetsPhi[indexGoodJets.at(i)],ReducedTree->Jets_ECorr[indexGoodJets.at(i)]);
	    VBF2.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[indexGoodJets.at(ii)],ReducedTree->JetsEta[indexGoodJets.at(ii)],ReducedTree->JetsPhi[indexGoodJets.at(ii)],ReducedTree->Jets_ECorr[indexGoodJets.at(ii)]);
	    TOT = VBF1 + VBF2;
	    if (TOT.Pt() < tempPtMax) continue;
	    tempPtMax = TOT.Pt(); //take the jet pair with largest Pt
	    nVBF1 = indexGoodJets.at(i); //save position of the 1st vbf jet
	    nVBF2 = indexGoodJets.at(ii); //save position of the 2nd vbf jet
	  }
	}
	
	if (nVBF1!=-1 && nVBF2!=-1) //save infos for vbf jet pair
	  {
	    VBF1.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[nVBF1],ReducedTree->JetsEta[nVBF1],ReducedTree->JetsPhi[nVBF1],ReducedTree->Jets_ECorr[nVBF1]);
	    VBF2.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[nVBF2],ReducedTree->JetsEta[nVBF2],ReducedTree->JetsPhi[nVBF2],ReducedTree->Jets_ECorr[nVBF2]);
	    TOT = VBF1 + VBF2;
	    
	    WWTree->vbf_maxpt_j1_pt = ReducedTree->Jets_PtCorr[nVBF1];
	    WWTree->vbf_maxpt_j1_eta = ReducedTree->JetsEta[nVBF1];
	    WWTree->vbf_maxpt_j1_phi = ReducedTree->JetsPhi[nVBF1];
	    WWTree->vbf_maxpt_j1_e = ReducedTree->Jets_ECorr[nVBF1];
	    WWTree->vbf_maxpt_j1_bDiscriminatorCSV = ReducedTree->Jets_bDiscriminatorICSV[nVBF1];
	    WWTree->vbf_maxpt_j2_pt = ReducedTree->Jets_PtCorr[nVBF2];
	    WWTree->vbf_maxpt_j2_eta = ReducedTree->JetsEta[nVBF2];
	    WWTree->vbf_maxpt_j2_phi = ReducedTree->JetsPhi[nVBF2];
	    WWTree->vbf_maxpt_j2_e = ReducedTree->Jets_ECorr[nVBF2];
	    WWTree->vbf_maxpt_j2_bDiscriminatorCSV = ReducedTree->Jets_bDiscriminatorICSV[nVBF2];
	    WWTree->vbf_maxpt_jj_pt = TOT.Pt();
	    WWTree->vbf_maxpt_jj_eta = TOT.Eta();
	    WWTree->vbf_maxpt_jj_phi = TOT.Phi();
	    WWTree->vbf_maxpt_jj_m = TOT.M();	
	  }
	
      }

    /////////////////MC Infos
    if (isMC)
      {
	TLorentzVector temp, temp2;
	//	std::cout<<"entry: "<<iEntry<<" "<<GenNuNum<<std::endl;
	double deltaPhiOld=100.;
	for (int i=0; i<ReducedTree->GenBosonNum; i++) {
	  double deltaPhi = getDeltaPhi(ReducedTree->GenBosonPhi[i],WWTree->W_phi);
	  if (abs(deltaPhi)>abs(deltaPhiOld))   continue;
	  //	  std::cout<<"bosone: "<<i<<" "<<ReducedTree->GenBosonPhi[i]<<" "<<W_phi<<std::endl;
	  temp.SetPtEtaPhiE(ReducedTree->GenBosonPt[i],ReducedTree->GenBosonEta[i],ReducedTree->GenBosonPhi[i],ReducedTree->GenBosonE[i]);
	  WWTree->W_pt_gen = ReducedTree->GenBosonPt[i];
	  WWTree->W_pz_gen = temp.Pz();
	  WWTree->W_rap_gen = temp.Rapidity();
	  deltaPhiOld = deltaPhi;
	}	
	if (ReducedTree->GenBosonNum==2) {
	  temp.SetPtEtaPhiE(ReducedTree->GenBosonPt[0],ReducedTree->GenBosonEta[0],ReducedTree->GenBosonPhi[0],ReducedTree->GenBosonE[0]);
	  temp2.SetPtEtaPhiE(ReducedTree->GenBosonPt[1],ReducedTree->GenBosonEta[1],ReducedTree->GenBosonPhi[1],ReducedTree->GenBosonE[1]);
	  WWTree->genGravMass=(temp+temp2).M();	
	}

	deltaPhiOld=100.;
       	for (int i=0; i<ReducedTree->GenNuNum; i++) {
	  double deltaPhi = getDeltaPhi(ReducedTree->GenNuPhi[i],WWTree->W_phi);
	  if (abs(deltaPhi)>abs(deltaPhiOld))   continue;	  
	  temp.SetPtEtaPhiE(ReducedTree->GenNuPt[i],ReducedTree->GenNuEta[i],ReducedTree->GenNuPhi[i],ReducedTree->GenNuE[i]);
	  WWTree->nu_pz_gen=temp.Pz();	  
	  WWTree->nu_pt_gen=temp.Pt();	  
	  WWTree->nu_phi_gen=temp.Phi();	  
	  WWTree->nu_eta_gen=temp.Eta();
	  deltaPhiOld = deltaPhi;
	}		
      }
    
    //fill the tree
    if(WWTree->event==evento) { std::cout<<"debug: "<<std::endl;
    std::cout<<WWTree->issignal <<" "<< 
      WWTree->W_pt <<" "<<
      WWTree->pfMET <<" "<<
      WWTree->l_pt <<" "<<
      WWTree->jet_pt<< " "<< 
      WWTree->nbtag <<" "<< 
      WWTree->jet_mass_pruned << " "<< 
      WWTree->jet_tau2tau1<<std::endl;
      }
    
    if(strcmp(leptonName.c_str(),"mu")==0 && WWTree->issignal==1 && WWTree->W_pt>200 && WWTree->pfMET>40 && WWTree->l_pt>50 && WWTree->jet_pt>200 && WWTree->nbtag <1 && WWTree->jet_mass_pruned > 40 && WWTree->jet_mass_pruned < 130 && WWTree->jet_tau2tau1<0.5) {
      outTree->Fill();
            if(WWTree->event==evento) { std::cout<<"debug: "<<std::endl; getchar(); }
    }
    if(strcmp(leptonName.c_str(),"el")==0 && WWTree->issignal==1 && WWTree->W_pt>200 && WWTree->pfMET>80 && WWTree->l_pt>90 && WWTree->jet_pt>200 && WWTree->nbtag <1 && WWTree->jet_mass_pruned > 40 && WWTree->jet_mass_pruned < 130 && WWTree->jet_tau2tau1<0.5) {
      outTree->Fill();
            if(WWTree->event==evento) { std::cout<<"debug: "<<std::endl; getchar(); }
    }
  }

  std::cout<<"lepton eff: "<<cutEff[0]<<std::endl
	   <<"met eff:    "<<cutEff[1]<<std::endl
	   <<"W eff:      "<<cutEff[2]<<std::endl
	   <<"AK8 found:  "<<cutEff[3]<<std::endl
	   <<"AK8 tight:  "<<cutEff[4]<<std::endl;

  //--------close everything-------------
  ReducedTree->fChain->Delete();
  outTree->Write();
  outROOT->Close();

  return(0);
}
