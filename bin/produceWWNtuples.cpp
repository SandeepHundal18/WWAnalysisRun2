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
#include "../interface/setOutputTree.h"
#include "../interface/METzCalculator.h"
#include "../interface/METzCalculator_Run2.h"
#include "../interface/analysisUtils.h"
#include "../interface/readJSONFile.h"

using namespace std;

//*****PU WEIGHT***************

vector<double> generate_weights(TH1* data_npu_estimated, int isForSynch){
  // see SimGeneral/MixingModule/python/mix_2015_25ns_Startup_PoissonOOTPU_cfi.pyy; copy and paste from there:
  const double npu_probs[52] = {
                        4.8551E-07,
                        1.74806E-06,
                        3.30868E-06,
                        1.62972E-05,
                        4.95667E-05,
                        0.000606966,
                        0.003307249,
                        0.010340741,
                        0.022852296,
                        0.041948781,
                        0.058609363,
                        0.067475755,
                        0.072817826,
                        0.075931405,
                        0.076782504,
                        0.076202319,
                        0.074502547,
                        0.072355135,
                        0.069642102,
                        0.064920999,
                        0.05725576,
                        0.047289348,
                        0.036528446,
                        0.026376131,
                        0.017806872,
                        0.011249422,
                        0.006643385,
                        0.003662904,
                        0.001899681,
                        0.00095614,
                        0.00050028,
                        0.000297353,
                        0.000208717,
                        0.000165856,
                        0.000139974,
                        0.000120481,
                        0.000103826,
                        8.88868E-05,
                        7.53323E-05,
                        6.30863E-05,
                        5.21356E-05,
                        4.24754E-05,
                        3.40876E-05,
                        2.69282E-05,
                        2.09267E-05,
                        1.5989E-05,
                        4.8551E-06,
                        2.42755E-06,
                        4.8551E-07,
                        2.42755E-07,
                        1.21378E-07,
                        4.8551E-08
};
  
  if (isForSynch==0) { //OFFICIAL RECIPE
    vector<double> result(52);
    double s = 0.0;
    for(int npu=0; npu<52; ++npu){
      double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));                              
      result[npu] = npu_estimated / npu_probs[npu];
      s += npu_estimated;
    }
    // normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
    for(int npu=0; npu<52; ++npu){
      result[npu] /= s;
    }
    return result;
  }

  else { //THIS IS FOR THE SYNCH ONLY. THIS IS NOT THE OFFICIAL RECIPE!
    vector<double> result(60);
    for(int npu=0; npu<60; ++npu){
      if (data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu))==NULL)
	result[npu] = 0.;
      else {
	double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));                              
	result[npu] = npu_estimated;
      }
    }
    return result;
  }

}


//*******MAIN*******************************************************************

int main (int argc, char** argv)
{ 
  std::string inputFolder = argv[1];
  std::string outputFile = argv[2];
  int isMC = atoi(argv[3]);
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
  float genMass = atof(argv[9]);
  int applyTrigger = atoi(argv[10]);
  std::string jsonFileName = argv[11];

  //  std::string jsonFileName="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt";
  std::map<int, std::vector<std::pair<int, int> > > jsonMap;
  jsonMap = readJSONFile(jsonFileName);
  std::cout<<"JSON file: "<<jsonFileName<<std::endl;

  // define map with events                                                                                                                                     
  std::map<std::pair<int,std::pair<int,int> >,int> eventsMap;

  //applyTrigger=false;
  std::cout<<"apply trigger: "<<applyTrigger<<std::endl;

  TLorentzVector W,MET,LEP;
  TLorentzVector NU0,NU1,NU2;
  TLorentzVector NU0_jes_up, NU0_jes_dn;
  TLorentzVector JET, HADW, AK4;
  TLorentzVector JET_jes_up, JET_jes_dn;
  TLorentzVector VBF1,VBF2,TOT;
  TLorentzVector ELE,MU;

  std::vector<TLorentzVector> tightMuon;
  std::vector<TLorentzVector> looseMuon;
  std::vector<TLorentzVector> tightEle;
  std::vector<TLorentzVector> looseEle;

  int ok=0, total=0;
  int runno=260627;
  int lumo=524;
  int evento=942309369;
  int count=0;
	
  setInputTree *ReducedTree = new setInputTree (inputTreeName.c_str());
  ReducedTree->Init();


  char command1[3000];
  sprintf(command1, "xrd eoscms dirlist %s/%s/  | awk '{print \"root://eoscms.cern.ch/\"$5}' > listTemp_%s.txt", (inputFolder).c_str(), (inputFile).c_str(), outputFile.c_str());
  std::cout<<command1<<std::endl;
  system(command1);
  char list1[2000];
  sprintf (list1, "listTemp_%s.txt", outputFile.c_str());
  ifstream rootList (list1);

  int fileCounter=0;
  Long64_t totalEntries=0;

  while (!rootList.eof())
    {
      char iRun_tW[700];
      rootList >> iRun_tW;
      ReducedTree->fChain->Add(iRun_tW);
      fileCounter++;
    }

  std::cout<<"number of files found: "<<fileCounter-2<<std::endl;
  std::cout<<"total entries: "<<ReducedTree->fChain->GetEntries()<<std::endl;
  totalEntries=ReducedTree->fChain->GetEntries();

  char command3[300];
  sprintf(command3, "rm listTemp_%s.txt", outputFile.c_str());
  system(command3);

  int cutEff[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  //--------pile up file -----------------
  //    TFile* pileupFile = TFile::Open("190456-208686-13Julv2_Prompt_Moriond2013.69400.observed.root");  
  //TH1F *pileupHisto = (TH1F*)pileupFile->Get("pileup");

  std::vector<double> weights_pu1; //these are made with our recipe
  std::vector<double> weights_pu2; //these are made with the official recipe

  TFile* pileupFile1 = TFile::Open("pileupDataRun2015D_72mb.root");  
  TH1F *pileupHisto1 = (TH1F*)pileupFile1->Get("pileup");  
  weights_pu1 = generate_weights(pileupHisto1,0);
  pileupFile1->Close();

  //  TFile* pileupFile2 = TFile::Open("puweights.root");  
  TFile* pileupFile2 = TFile::Open("PUxSynch.root");  
  TH1F *pileupHisto2 = (TH1F*)pileupFile2->Get("puweights");
  weights_pu2 = generate_weights(pileupHisto2,1);
  pileupFile2->Close();


  //---------output tree----------------
  TFile* outROOT = TFile::Open((std::string("output/output_")+leptonName+std::string("/")+outputFile+(".root")).c_str(),"recreate");
  outROOT->cd();
  TTree* outTree = new TTree("otree", "otree");
  outTree->SetDirectory(0);

  setOutputTree *WWTree = new setOutputTree(outTree);

  float top_NNLO_weight[2];

  std::ifstream badEventsFile;
  std::multimap<int,int> badEventsList;

  int run, lumi, evt;
  if (isMC==0) {
    if (strcmp(leptonName.c_str(),"el")==0)
      badEventsFile.open("SingleElectron_csc2015.txt");
    else
      badEventsFile.open("SingleMuon_csc2015.txt");      
    while(!badEventsFile.eof()) 
      {
	badEventsFile >> run >> lumi >> evt;
	badEventsList.insert(std::pair<int,int>(run,evt));
      }      
  }
  badEventsFile.close();

  //---------start loop on events------------
  Long64_t jentry2=0;
  for (Long64_t jentry=0; jentry<ReducedTree->fChain->GetEntries();jentry++,jentry2++) {
    //for (Long64_t jentry=531000; jentry<532000;jentry++,jentry2++) {

    Long64_t iEntry = ReducedTree->LoadTree(jentry);
    if (iEntry < 0) break;
    int nb = ReducedTree->fChain->GetEntry(jentry);   
    // if (Cut(ientry) < 0) continue;                                                                                                                           
    
    tightMuon.clear();
    tightEle.clear();
    looseMuon.clear();
    looseEle.clear();

    if(jentry2 % 1000 == 0)    
      std::cout << "read entry: " << jentry2 <<"/"<<totalEntries<<std:: endl;

    //*********************************                                                                                                                       
    // JSON FILE AND DUPLIACTES IN DATA                                                                                                                       
    WWTree->run   = ReducedTree->RunNum;
    WWTree->event = ReducedTree->EvtNum;
    WWTree->lumi = ReducedTree->LumiBlockNum;

    bool skipEvent = false;

    if( isMC==0 )         //apply json file
        {                                                                                                                                                        
	  if(AcceptEventByRunAndLumiSection(ReducedTree->RunNum,ReducedTree->LumiBlockNum,jsonMap) == false) skipEvent = true;                           
          std::pair<int,Long64_t> eventLSandID(ReducedTree->LumiBlockNum,ReducedTree->EvtNum);            
	  std::pair<int,std::pair<int,Long64_t> > eventRUNandLSandID(ReducedTree->RunNum,eventLSandID); 
          if( eventsMap[eventRUNandLSandID] == 1 ) skipEvent = true;                                             
          else eventsMap[eventRUNandLSandID] = 1;                                                          
          }
    if( skipEvent == true ) continue;    

    WWTree->initializeVariables(); //initialize all variables

    if (ReducedTree->passFilterHBHELooseRerun == 0) continue;
    if (ReducedTree->passFilterHBHEIsoRerun == 0) continue;
    if (ReducedTree->passFilterCSCHalo == 0) continue;
    if (ReducedTree->passFilterGoodVtx == 0) continue;
    if (ReducedTree->passFilterEEBadSC == 0) continue;
    
    WWTree->issignal = 0;
    WWTree->wSampleWeight = weight; //xsec/numberOfEntries
    WWTree->eff_and_pu_Weight = 1.; //temporary value
    WWTree->eff_and_pu_Weight_2 = 1.; //temporary value
    WWTree->eff_and_pu_Weight_3 = 1.; //temporary value
    WWTree->top1_NNLO_Weight = 1.;
    WWTree->top2_NNLO_Weight = 1.;
    WWTree->trig_eff_Weight = 1.;

    if (ReducedTree->genEventWeight>0)
      WWTree->genWeight=1.;
    else if (ReducedTree->genEventWeight<0)
      WWTree->genWeight=-1.;
    //    WWTree->genWeight = ReducedTree->genEventWeight;

    //PILE-UP WEIGHT
    if (isMC==1) {
      if(ReducedTree->NVtx<weights_pu1.size()){
	WWTree->eff_and_pu_Weight = weights_pu1[ReducedTree->npT]; //official pu recipe
      }
      else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
	std::cout<<"Warning! n_pu too big"<<std::endl;
	//	throw logic_error("n_pu too big");
        WWTree->eff_and_pu_Weight = 0.;
      }    

      if(ReducedTree->NVtx<weights_pu2.size()){
	WWTree->eff_and_pu_Weight_2 = weights_pu2[ReducedTree->NVtx]; //our pu recipe
      }
      else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
	std::cout<<"Warning! n_pu too big"<<std::endl;
	//	throw logic_error("n_pu too big");
        WWTree->eff_and_pu_Weight_2 = 0.;
      }    
      WWTree->eff_and_pu_Weight_3 = ReducedTree->PUWeight; //our pu recipe
    }    
    //require at least one lepton and one jet
    //    if ( strcmp(leptonName.c_str(),"el")==0 && ReducedTree->ElectronsNum==0) continue; 
    //    if ( strcmp(leptonName.c_str(),"mu")==0 && ReducedTree->MuonsNum==0) continue;      
        
    //save event variables
    WWTree->run   = ReducedTree->RunNum;
    WWTree->event = ReducedTree->EvtNum;
    WWTree->lumi = ReducedTree->LumiBlockNum;

    
   // WWTree->njets = ReducedTree->NJets;
    WWTree->nPV  = ReducedTree->NVtx;
    count=0;
    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug: "<<count<<std::endl; count++;
    
    /////////////////THE SELECTED LEPTON
    int nTightLepton=0;
    if (strcmp(leptonName.c_str(),"el")==0) {
      int passTrigger=0;
      float tempPt=0.;
      for (int i=0; i<ReducedTree->ElectronsNum; i++) {
	if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug ele: "<<i<<std::endl;
	if (applyTrigger==1)
	  for (int t=0; t<ReducedTree->TriggerProducerTriggerNames->size(); t++)
	    if(TString(ReducedTree->TriggerProducerTriggerNames->at(t)).Contains("HLT_Ele105_CaloIdVT_GsfTrkIdT") || 
	       TString(ReducedTree->TriggerProducerTriggerNames->at(t)).Contains("HLT_Ele115_CaloIdVT_GsfTrkIdT"))
	      if (ReducedTree->TriggerProducerTriggerPass->at(t)==1) passTrigger=1; //trigger
	if (passTrigger==0) continue;
	if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug ele: "<<i<<std::endl;
	//if (ReducedTree->TriggerProducerTriggerPass->at(0)==0) continue; //trigger
	if (ReducedTree->Electrons_isHEEP[i]==false) continue;       
	if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug ele: "<<i<<std::endl;
        if (ReducedTree->ElectronsPt[i]<=120) continue;
	if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug ele: "<<i<<std::endl;
	//        if (fabs(ReducedTree->ElectronsEta[i])>=2.5) continue; //this is already in the HEEP requirement
	//if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug ele: "<<i<<std::endl;
	if (ReducedTree->ElectronsPt[i]<tempPt) continue;
	if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug ele: "<<i<<std::endl;
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
      int passTrigger=0;
      float tempPt=0.;
      for (int i=0; i<ReducedTree->MuonsNum; i++) {
	if (applyTrigger==1)
	  for (int t=0; t<ReducedTree->TriggerProducerTriggerNames->size(); t++)
	    if(TString(ReducedTree->TriggerProducerTriggerNames->at(t)).Contains("HLT_Mu45_eta2p1") || 
	       TString(ReducedTree->TriggerProducerTriggerNames->at(t)).Contains("HLT_Mu50"))
	      if (ReducedTree->TriggerProducerTriggerPass->at(t)==1) passTrigger=1; //trigger
	if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug mu: "<<i<<std::endl;
	if (passTrigger==0) continue;
	if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug mu: "<<i<<std::endl;
	//if (ReducedTree->TriggerProducerTriggerPass->at(1)==0) continue; //trigger
	if (ReducedTree->Muons_isHighPt[i]==false) continue;
	if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug mu: "<<i<<std::endl;
	//	if (ReducedTree->Muons_isPFMuon[i]==false) continue; //not in the synch ntuple!!
        if ((ReducedTree->Muons_trackIso[i]/ReducedTree->MuonsPt[i])>=0.1) continue;
	if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug mu: "<<i<<std::endl;
        if (ReducedTree->MuonsPt[i]<53) continue;
	if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug mu: "<<i<<std::endl;
        if (fabs(ReducedTree->MuonsEta[i])>=2.1) continue;
	if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug mu: "<<i<<std::endl;
	MU.SetPtEtaPhiE(ReducedTree->MuonsPt[i],ReducedTree->MuonsEta[i],ReducedTree->MuonsPhi[i],ReducedTree->MuonsE[i]);
	tightMuon.push_back(MU);
	if (ReducedTree->MuonsPt[i]<tempPt) continue;
	if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug mu: "<<i<<std::endl;
	WWTree->l_pt  = ReducedTree->MuonsPt[i];
	WWTree->l_eta = ReducedTree->MuonsEta[i];
	WWTree->l_phi = ReducedTree->MuonsPhi[i];
	WWTree->l_e = ReducedTree->MuonsE[i];
	tempPt = WWTree->l_pt;
	nTightLepton++;
      }
    }
    if (nTightLepton==0) continue; //no leptons with required ID
    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug: "<<count<<std::endl; count++;

    if (strcmp(leptonName.c_str(),"mu")==0 && isMC==1) { //trigger SF for muon
      if ( fabs(WWTree->l_eta) < 0.9) {
	if      ( WWTree->l_pt> 50 && WWTree->l_pt<60) WWTree->trig_eff_Weight = 0.9764;
	else if ( WWTree->l_pt>=60) WWTree->trig_eff_Weight = 0.9693;
      }	
      else if ( fabs(WWTree->l_eta) >= 0.9 && fabs(WWTree->l_eta)<1.2 ) {
	if      ( WWTree->l_pt> 50 && WWTree->l_pt<60) WWTree->trig_eff_Weight = 0.9722;
	else if ( WWTree->l_pt>=60) WWTree->trig_eff_Weight = 0.9607;
      }	
      else {
	if      ( WWTree->l_pt> 50 && WWTree->l_pt<60) WWTree->trig_eff_Weight = 0.9631;
	else if ( WWTree->l_pt>=60) WWTree->trig_eff_Weight = 0.9677;
      }	      
    }

    //VETO ADDITIONAL LEPTONS
    int nLooseLepton=0;
    for (int i=0; i<ReducedTree->ElectronsNum; i++) {
    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug loose el: "<<i<<std::endl; count++;
      if (ReducedTree->Electrons_isHEEP[i]==false) continue;       
      if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug loose el: "<<i<<", pt: "<<ReducedTree->ElectronsPt[i]<<", eta: "<<ReducedTree->ElectronsEta[i]<<std::endl; count++;
      if (ReducedTree->ElectronsPt[i]<35) continue;       
      //    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug: "<<i<<std::endl; count++;
      //       if (fabs(ReducedTree->ElectronsEta[i])>=2.5) continue;
    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug loose el: "<<i<<std::endl; count++;
      ELE.SetPtEtaPhiE(ReducedTree->ElectronsPt[i],ReducedTree->ElectronsEta[i],ReducedTree->ElectronsPhi[i],ReducedTree->ElectronsE[i]);
      looseEle.push_back(ELE);      
      nLooseLepton++;
    }
    for (int i=0; i<ReducedTree->MuonsNum; i++) {
    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug loose mu: "<<i<<std::endl; count++;
      if (ReducedTree->Muons_isHighPt[i]==false) continue;
    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug loose mu: "<<i<<std::endl; count++;
      if ((ReducedTree->Muons_trackIso[i]/ReducedTree->MuonsPt[i])>=0.1) continue;
    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug loose mu: "<<i<<std::endl; count++;
      if (fabs(ReducedTree->MuonsEta[i])>=2.4) continue;
    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug loose mu: "<<i<<std::endl; count++;
      if (ReducedTree->MuonsPt[i]<20) continue;
    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug loose mu: "<<i<<std::endl; count++;
      MU.SetPtEtaPhiE(ReducedTree->MuonsPt[i],ReducedTree->MuonsEta[i],ReducedTree->MuonsPhi[i],ReducedTree->MuonsE[i]);
      looseMuon.push_back(MU);
      nLooseLepton++;
    }
    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo)     std::cout<<"nlooseleptons: "<<nLooseLepton<<std::endl;
    if (nLooseLepton!=1) continue; //no additional leptons
    cutEff[0]++;
    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug: "<<count<<std::endl; count++;

    //preselection on jet pt and met
    if (ReducedTree->METPt < 30) continue; 
    cutEff[1]++;
    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug: "<<count<<std::endl; count++;

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

    TLorentzVector W_mu, W_Met, W_Met_jes_up, W_Met_jes_dn;

    W_mu.SetPtEtaPhiE(WWTree->l_pt,WWTree->l_eta,WWTree->l_phi,WWTree->l_e);
    W_Met.SetPxPyPzE(ReducedTree->METPt * TMath::Cos(ReducedTree->METPhi), ReducedTree->METPt * TMath::Sin(ReducedTree->METPhi), 0., sqrt(ReducedTree->METPt*ReducedTree->METPt));
    W_Met_jes_up.SetPxPyPzE(ReducedTree->METPtUp * TMath::Cos(ReducedTree->METPhiUp), ReducedTree->METPtUp * TMath::Sin(ReducedTree->METPhiUp), 0., sqrt(ReducedTree->METPtUp*ReducedTree->METPtUp));
    W_Met_jes_dn.SetPxPyPzE(ReducedTree->METPtDown * TMath::Cos(ReducedTree->METPhiDown), ReducedTree->METPtDown * TMath::Sin(ReducedTree->METPhiDown), 0., sqrt(ReducedTree->METPtDown*ReducedTree->METPtDown));

    if(W_mu.Pt()<=0 || W_Met.Pt() <= 0 ){ std::cerr<<" Negative Lepton - Neutrino Pt "<<std::endl; continue ; }
    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug: "<<count<<std::endl; count++;

    // type0 calculation of neutrino pZ
    METzCalculator NeutrinoPz_type0;
    METzCalculator NeutrinoPz_type0_jes_up;
    METzCalculator NeutrinoPz_type0_jes_dn;
    METzCalculator_Run2 NeutrinoPz_run2;
    NeutrinoPz_type0.SetMET(W_Met);
    NeutrinoPz_type0.SetLepton(W_mu);
    NeutrinoPz_type0.SetLeptonType(leptonName.c_str());

    NeutrinoPz_type0_jes_up.SetMET(W_Met_jes_up);
    NeutrinoPz_type0_jes_up.SetLepton(W_mu);
    NeutrinoPz_type0_jes_up.SetLeptonType(leptonName.c_str());

    NeutrinoPz_type0_jes_dn.SetMET(W_Met_jes_dn);
    NeutrinoPz_type0_jes_dn.SetLepton(W_mu);
    NeutrinoPz_type0_jes_dn.SetLeptonType(leptonName.c_str());

    NeutrinoPz_run2.SetMET(W_Met);
    NeutrinoPz_run2.SetLepton(W_mu);
    NeutrinoPz_run2.SetLeptonType(leptonName.c_str());

    double pz1_type0 = NeutrinoPz_type0.Calculate(); // Default one -> according to type0
    double pz2_type0 = NeutrinoPz_type0.getOther(); // Default one

    double pz1_run2 = NeutrinoPz_run2.Calculate(); 

    double pz1_type0_jes_up = NeutrinoPz_type0_jes_up.Calculate(); // Default one -> according to type0
    double pz1_type0_jes_dn = NeutrinoPz_type0_jes_dn.Calculate(); // Default one -> according to type0

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
    WWTree->pfMET_jes_up   = sqrt(ReducedTree->METPtUp*ReducedTree->METPtUp);
    WWTree->pfMET_jes_dn   = sqrt(ReducedTree->METPtDown*ReducedTree->METPtDown);
    WWTree->pfMET_Phi = ReducedTree->METPhi;
    WWTree->nu_pz_type0 = pz1_type0;
    WWTree->nu_pz_type2 = pz1_type2;
    WWTree->nu_pz_run2 = pz1_run2;
    WWTree->nu_pz_isre = 1-NeutrinoPz_run2.IsComplex();
    WWTree->nu_pz_run2_oth = NeutrinoPz_run2.getOther();
    WWTree->nu_pz_run2_type = NeutrinoPz_run2.getType();


    /////////////////THE LEPTONIC W
    
    LEP.SetPtEtaPhiE(WWTree->l_pt,WWTree->l_eta,WWTree->l_phi,WWTree->l_e);
    NU0.SetPxPyPzE(ReducedTree->METPt*TMath::Cos(ReducedTree->METPhi),ReducedTree->METPt*TMath::Sin(ReducedTree->METPhi),WWTree->nu_pz_type0,TMath::Sqrt(WWTree->pfMET*WWTree->pfMET+WWTree->nu_pz_type0*WWTree->nu_pz_type0));

    NU0_jes_up.SetPxPyPzE(ReducedTree->METPtUp*TMath::Cos(ReducedTree->METPhiUp),ReducedTree->METPtUp*TMath::Sin(ReducedTree->METPhiUp),pz1_type0_jes_up,TMath::Sqrt(WWTree->pfMET_jes_up*WWTree->pfMET_jes_up+pz1_type0_jes_up*pz1_type0_jes_up));
    NU0_jes_dn.SetPxPyPzE(ReducedTree->METPtDown*TMath::Cos(ReducedTree->METPhiDown),ReducedTree->METPtDown*TMath::Sin(ReducedTree->METPhiDown),pz1_type0_jes_dn,TMath::Sqrt(WWTree->pfMET_jes_dn*WWTree->pfMET_jes_dn+pz1_type0_jes_dn*pz1_type0_jes_dn));

    NU2.SetPxPyPzE(ReducedTree->METPt*TMath::Cos(ReducedTree->METPhi),ReducedTree->METPt*TMath::Sin(ReducedTree->METPhi),WWTree->nu_pz_type2,TMath::Sqrt(WWTree->pfMET*WWTree->pfMET+WWTree->nu_pz_type2*WWTree->nu_pz_type2));
    NU1.SetPxPyPzE(ReducedTree->METPt*TMath::Cos(ReducedTree->METPhi),ReducedTree->METPt*TMath::Sin(ReducedTree->METPhi),WWTree->nu_pz_run2,TMath::Sqrt(WWTree->pfMET*WWTree->pfMET+WWTree->nu_pz_run2*WWTree->nu_pz_run2));
    W = LEP + NU2;
    
    WWTree->v_pt = W.Pt();
    WWTree->v_eta = W.Eta();
    WWTree->v_phi = W.Phi();
    WWTree->v_mt = TMath::Sqrt(2*LEP.Et()*NU2.Et()*(1-TMath::Cos(LEP.DeltaPhi(NU2))));
    //    W_mt = W.Mt();


    //FOR THE SYNCHORNIZATION!!! REMOVE IT FOR THE REAL ANALYSIS!!!!
    //    NU2.SetPtEtaPhiE(ReducedTree->METPt,0.,ReducedTree->METPhi,0.);
    //    W = NU2+LEP; 
    ////

    if (W.Pt()<100) continue;
    cutEff[2]++;
    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug: "<<count<<std::endl; count++;


    //    if (WWTree->v_pt < 150) continue;
//    if (WWTree->deltaR_lak8jet < (TMath::Pi()/2.0))   continue;

    ///////////THE FAT JET
    float tempPt=0., tempMass=0.;
    int nGoodAK8jets=0;
    int hadWpos = -1;
    int ttb_jet_position=-1; //position of AK8 jet in ttbar-topology
    if (ReducedTree->AK8JetsNum < 1 ) continue; 
    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug: "<<count<<std::endl; count++;
   
    for (unsigned int i=0; i<ReducedTree->AK8JetsNum; i++)
      {
	bool isCleanedJet = true;
	if (ReducedTree->AK8Jets_PtCorr[i]<100 || fabs(ReducedTree->AK8JetsEta[i])>2.4)  continue; //be careful: this is not inside the synchntuple code
	if (ReducedTree->AK8Jets_prunedMass[i]>tempMass) {
	  if ( (ReducedTree->AK8JetsEta[i]>0 && WWTree->l_eta<0) || 
	       (ReducedTree->AK8JetsEta[i]<0 && WWTree->l_eta>0)) { //jet and lepton in opposite hemisphere for ttb
	    ttb_jet_position=i; //save AK8 jet in ttbar topology
	    tempMass=ReducedTree->AK8Jets_prunedMass[i];
	  }
	}
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

	WWTree->ungroomed_jet_pt  = ReducedTree->AK8Jets_PtCorr[i];
	WWTree->ungroomed_jet_pt_jes_up = (ReducedTree->AK8Jets_PtCorr[i]/ReducedTree->AK8Jets_AK8correction[i])*ReducedTree->AK8Jets_AK8correctionUp[i];
	WWTree->ungroomed_jet_pt_jes_dn = (ReducedTree->AK8Jets_PtCorr[i]/ReducedTree->AK8Jets_AK8correction[i])*ReducedTree->AK8Jets_AK8correctionDown[i];
	WWTree->ungroomed_jet_eta = ReducedTree->AK8JetsEta[i];
	WWTree->ungroomed_jet_phi = ReducedTree->AK8JetsPhi[i];
	WWTree->ungroomed_jet_e   = ReducedTree->AK8Jets_ECorr[i];
	WWTree->jet_mass_pr   = ReducedTree->AK8Jets_prunedMass[i];
	WWTree->jet_mass_pr_jes_up = (ReducedTree->AK8Jets_prunedMass[i]/ReducedTree->AK8Jets_AK8massCorrection[i])*ReducedTree->AK8Jets_AK8massCorrectionUp[i];
	WWTree->jet_mass_pr_jes_dn = (ReducedTree->AK8Jets_prunedMass[i]/ReducedTree->AK8Jets_AK8massCorrection[i])*ReducedTree->AK8Jets_AK8massCorrectionDown[i];
        WWTree->jet_mass_so   = ReducedTree->AK8Jets_softDropMass[i];
        WWTree->jet_pt_so   = ReducedTree->AK8Jets_softDropPt[i];
	WWTree->jet_mass_tr   = ReducedTree->AK8Jets_trimmedMass[i];
	WWTree->jet_mass_fi   = ReducedTree->AK8Jets_filteredMass[i];
	WWTree->jet_tau2tau1   = ReducedTree->AK8Jets_tau2[i]/ReducedTree->AK8Jets_tau1[i];
	tempPt = WWTree->ungroomed_jet_pt;
	nGoodAK8jets++;
	hadWpos = i;
      }

    if (nGoodAK8jets==0) continue; //not found a good hadronic W candidate
    cutEff[3]++;
    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug: "<<count<<std::endl; count++;

    if (WWTree->ungroomed_jet_pt<100) continue;
    cutEff[4]++;
    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug: "<<count<<std::endl; count++;

    //////////////////ANGULAR VARIABLES
    JET.SetPtEtaPhiE(WWTree->ungroomed_jet_pt,WWTree->ungroomed_jet_eta,WWTree->ungroomed_jet_phi,WWTree->ungroomed_jet_e);
    WWTree->deltaR_lak8jet = JET.DeltaR(LEP);
    WWTree->deltaphi_METak8jet = JET.DeltaPhi(NU2);
    WWTree->deltaphi_Vak8jet = JET.DeltaPhi(W);
    if (WWTree->deltaR_lak8jet>(TMath::Pi()/2.0) && fabs(WWTree->deltaphi_METak8jet)>2.0 && fabs(WWTree->deltaphi_Vak8jet)>2.0)
      WWTree->issignal=1;
    JET_jes_up.SetPtEtaPhiE(WWTree->ungroomed_jet_pt*(ReducedTree->AK8Jets_AK8correctionUp[hadWpos]/ReducedTree->AK8Jets_AK8correction[hadWpos]),
		     WWTree->ungroomed_jet_eta,
		     WWTree->ungroomed_jet_phi,
		     WWTree->ungroomed_jet_e*(ReducedTree->AK8Jets_AK8correctionUp[hadWpos]/ReducedTree->AK8Jets_AK8correction[hadWpos]));
    JET_jes_dn.SetPtEtaPhiE(WWTree->ungroomed_jet_pt*(ReducedTree->AK8Jets_AK8correctionDown[hadWpos]/ReducedTree->AK8Jets_AK8correction[hadWpos]),
		     WWTree->ungroomed_jet_eta,
		     WWTree->ungroomed_jet_phi,
		     WWTree->ungroomed_jet_e*(ReducedTree->AK8Jets_AK8correctionDown[hadWpos]/ReducedTree->AK8Jets_AK8correction[hadWpos]));

    //FOUR-BODY INVARIANT MASS
    WWTree->mass_lvj_type0 = (LEP + NU0 + JET).M();
    WWTree->mass_lvj_type2 = (LEP + NU2 + JET).M();
    WWTree->mass_lvj_run2  = (LEP + NU1 + JET).M();
    WWTree->mass_lvj_type0_met_jes_up = (LEP + NU0_jes_up + JET_jes_up).M();
    WWTree->mass_lvj_type0_met_jes_dn = (LEP + NU0_jes_dn + JET_jes_dn).M();

    //--- ttbar topology ------
    if (ttb_jet_position>=0) {
      WWTree->ttb_ungroomed_jet_pt  = ReducedTree->AK8Jets_PtCorr[ttb_jet_position];
      WWTree->ttb_ungroomed_jet_eta = ReducedTree->AK8JetsEta[ttb_jet_position];
      WWTree->ttb_ungroomed_jet_phi = ReducedTree->AK8JetsPhi[ttb_jet_position];
      WWTree->ttb_ungroomed_jet_e   = ReducedTree->AK8Jets_ECorr[ttb_jet_position];
      WWTree->ttb_jet_mass_pr   = ReducedTree->AK8Jets_prunedMass[ttb_jet_position];
      WWTree->ttb_jet_mass_so   = ReducedTree->AK8Jets_softDropMass[ttb_jet_position];
      WWTree->ttb_jet_pt_so   = ReducedTree->AK8Jets_softDropPt[ttb_jet_position];
      WWTree->ttb_jet_mass_tr   = ReducedTree->AK8Jets_trimmedMass[ttb_jet_position];
      WWTree->ttb_jet_mass_fi   = ReducedTree->AK8Jets_filteredMass[ttb_jet_position];
      WWTree->ttb_jet_tau2tau1   = ReducedTree->AK8Jets_tau2[ttb_jet_position]/ReducedTree->AK8Jets_tau1[ttb_jet_position];

      WWTree->ttb_deltaeta_lak8jet = deltaEta(WWTree->ttb_ungroomed_jet_eta,WWTree->l_eta);
    }


    /////////VBF and b-tag section
    bool fillVBF = true;

    HADW.SetPtEtaPhiE(WWTree->ungroomed_jet_pt,WWTree->ungroomed_jet_eta,WWTree->ungroomed_jet_phi,WWTree->ungroomed_jet_e); //AK8 fat jet (hadronic W)
    std::vector<int> indexGoodJets;
    indexGoodJets.clear();
    if (indexGoodJets.size()!=0)  fillVBF=false;

    WWTree->njets=0;
    WWTree->nBTagJet_loose=0;
    WWTree->nBTagJet_medium=0;
    WWTree->nBTagJet_tight=0;

    float oldDeltaR = 1000.;
    float oldDeltaRLep = 1000.;
    int indexCloserJet = -1;
    int indexCloserJetLep = -1;


    //    LEP.SetPtEtaPhiE(WWTree->l_pt,WWTree->l_eta,WWTree->l_phi,WWTree->l_e);
    //    NU2.SetPxPyPzE(ReducedTree->METPt*TMath::Cos(ReducedTree->METPhi),ReducedTree->METPt*TMath::Sin(ReducedTree->METPhi),WWTree->nu_pz_type2,TMath::Sqrt(WWTree->pfMET*WWTree->pfMET+WWTree->nu_pz_type2*WWTree->nu_pz_type2));
    //    W = LEP + NU2;
    float deltaRbtag_prev=100.;
    float deltaRbtag_prev_loose=100.;

    for (unsigned int i=0; i<ReducedTree->JetsNum; i++) //loop on AK4 jet
      {
	bool isCleanedJet = true;
	if (ReducedTree->Jets_PtCorr[i]<=30 || ReducedTree->JetsPt[i]<=20 || fabs(ReducedTree->JetsEta[i])>=2.4)  continue;
	if (ReducedTree->Jets_isLooseJetId[i]==false) continue;

	//CLEANING
	if (deltaR(WWTree->ungroomed_jet_eta, WWTree->ungroomed_jet_phi,
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


	WWTree->njets++;

	AK4.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[i],ReducedTree->JetsEta[i],ReducedTree->JetsPhi[i],ReducedTree->Jets_ECorr[i]);

	//fill B-Tag info
	if (ReducedTree->Jets_bDiscriminatorICSV[i]>0.605) { 
	  WWTree->nBTagJet_loose++;
	  float deltaRbtag = HADW.DeltaR(AK4);
	  if (deltaRbtag>0.8 && deltaRbtag<deltaRbtag_prev_loose) {
	    WWTree->deltaR_AK8_closestBtagJet_loose = deltaRbtag;
	    deltaRbtag_prev_loose = deltaRbtag;
	  }	  
	}

	if (ReducedTree->Jets_bDiscriminatorICSV[i]>0.890) {  
	  WWTree->nBTagJet_medium++;
	  float deltaRbtag = HADW.DeltaR(AK4);
	  if (deltaRbtag>0.8 && deltaRbtag<deltaRbtag_prev) {
	    WWTree->deltaR_AK8_closestBtagJet = deltaRbtag;
	    deltaRbtag_prev = deltaRbtag;
	  }	  
	}
	if (ReducedTree->Jets_bDiscriminatorICSV[i]>0.970)   WWTree->nBTagJet_tight++;

	float deltaRlep = W.DeltaR(AK4);
	if (deltaRlep<oldDeltaRLep) indexCloserJetLep = i;

	float deltaR = HADW.DeltaR(AK4);
	if (deltaR<0.8) continue; //the vbf jets must be outside the had W cone
	
	if (WWTree->njets!=0) {
	  if (WWTree->jet2_pt!=0) {
	    WWTree->jet3_pt=ReducedTree->Jets_PtCorr[i];
	    WWTree->jet3_eta=ReducedTree->JetsEta[i];
	    WWTree->jet3_phi=ReducedTree->JetsPhi[i];
	    WWTree->jet3_e=ReducedTree->Jets_ECorr[i];
	    WWTree->jet3_btag=ReducedTree->Jets_bDiscriminatorICSV[i];
	  }
	  else {
	    WWTree->jet2_pt=ReducedTree->Jets_PtCorr[i];
	    WWTree->jet2_eta=ReducedTree->JetsEta[i];
	    WWTree->jet2_phi=ReducedTree->JetsPhi[i];
	    WWTree->jet2_e=ReducedTree->Jets_ECorr[i];
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
      WWTree->AK8_closerjet_pt = AK4.Pt();
      WWTree->AK8_closerjet_eta = AK4.Eta();
      WWTree->AK8_closerjet_phi = AK4.Phi();
      WWTree->AK8_closerjet_e = AK4.E();
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
    if (isMC==1)
      {
	TLorentzVector hadW, lepW, temp;
	int posWhad =-1, posWlep =-1, posTemp=-1, posGenJet=-1;
	//	std::cout<<"entry: "<<iEntry<<" "<<GenNuNum<<std::endl;
	double deltaPhiOld=100.;
	WWTree->genGravMass=100.;	

	for (int i=0; i<ReducedTree->GenBosonNum; i++) {
	  for (int j=i+1; j<ReducedTree->GenBosonNum; j++) {

	    hadW.SetPtEtaPhiE(ReducedTree->GenBosonPt[i],ReducedTree->GenBosonEta[i],ReducedTree->GenBosonPhi[i],ReducedTree->GenBosonE[i]);
	    lepW.SetPtEtaPhiE(ReducedTree->GenBosonPt[j],ReducedTree->GenBosonEta[j],ReducedTree->GenBosonPhi[j],ReducedTree->GenBosonE[j]);


	    if (fabs((hadW+lepW).M()-genMass)< fabs(WWTree->genGravMass-genMass)) { //found the gen graviton
	      WWTree->genGravMass=(hadW+lepW).M();	
	      posWhad=i; //save positions of the two W's in random order, will fix them later in the code
	      posWlep=j;

	    }

	  }	
	}

	if (posWhad!=-1 && posWlep!=-1) {

	  float oldDR=100.;
	  bool isWhadOk=true;

	  if(WWTree->event==127270357) std::cout<<"debug: "<<std::endl;

	  for (int i=0; i<ReducedTree->GenJetsAK8Num; i++) {
	      if (deltaR(ReducedTree->GenBosonEta[posWhad], ReducedTree->GenBosonPhi[posWhad],
			 ReducedTree->GenJetsAK8Eta[i], ReducedTree->GenJetsAK8Phi[i])< oldDR ) 
		{
		  posGenJet=i;		
		  oldDR = deltaR(ReducedTree->GenBosonEta[posWhad], ReducedTree->GenBosonPhi[posWhad],ReducedTree->GenJetsAK8Eta[i], ReducedTree->GenJetsAK8Phi[i]);
		  isWhadOk = true;
		  if(WWTree->event==127270357) std::cout<<"debug: had "<<ReducedTree->GenJetsAK8Phi[i]<<" "<<ReducedTree->GenBosonPhi[posWhad]<<" "<<oldDR<<std::endl;
		}
	      if (deltaR(ReducedTree->GenBosonEta[posWlep], ReducedTree->GenBosonPhi[posWlep],
			 ReducedTree->GenJetsAK8Eta[i], ReducedTree->GenJetsAK8Phi[i])< oldDR ) 
		{
		  posGenJet=i;
		  oldDR = deltaR(ReducedTree->GenBosonEta[posWlep], ReducedTree->GenBosonPhi[posWlep],ReducedTree->GenJetsAK8Eta[i], ReducedTree->GenJetsAK8Phi[i]);
		  isWhadOk = false;
		  if(WWTree->event==127270357) std::cout<<"debug: lep "<<ReducedTree->GenJetsAK8Phi[i]<<" "<<ReducedTree->GenBosonPhi[posWlep]<<" "<<oldDR<<std::endl;
         	}      
	  }
	  if(WWTree->event==127270357) std::cout<<ReducedTree->GenJetsAK8Phi[posGenJet]<<std::endl;
	  if (isWhadOk==false) //wrong W's positions saved, switch them
	    {
	      posTemp = posWhad;
	      posWhad = posWlep;
	      posWlep = posTemp;
	    }

	  hadW.SetPtEtaPhiE(ReducedTree->GenBosonPt[posWhad],ReducedTree->GenBosonEta[posWhad],ReducedTree->GenBosonPhi[posWhad],ReducedTree->GenBosonE[posWhad]);
	  lepW.SetPtEtaPhiE(ReducedTree->GenBosonPt[posWlep],ReducedTree->GenBosonEta[posWlep],ReducedTree->GenBosonPhi[posWlep],ReducedTree->GenBosonE[posWlep]);

	  WWTree->W_pt_gen = ReducedTree->GenBosonPt[posWlep];
	  WWTree->W_pz_gen = lepW.Pz();
	  WWTree->W_rap_gen = lepW.Rapidity();
	  
	  WWTree->hadW_pt_gen = ReducedTree->GenBosonPt[posWhad];
	  WWTree->hadW_eta_gen = ReducedTree->GenBosonEta[posWhad];
	  WWTree->hadW_phi_gen = ReducedTree->GenBosonPhi[posWhad];
	  WWTree->hadW_e_gen = ReducedTree->GenBosonE[posWhad];
	  WWTree->hadW_m_gen = hadW.M();

	  WWTree->lepW_pt_gen = ReducedTree->GenBosonPt[posWlep];
	  WWTree->lepW_eta_gen = ReducedTree->GenBosonEta[posWlep];
	  WWTree->lepW_phi_gen = ReducedTree->GenBosonPhi[posWlep];
	  WWTree->lepW_e_gen = ReducedTree->GenBosonE[posWlep];
	  WWTree->lepW_m_gen = lepW.M();

	  WWTree->AK8_pt_gen = ReducedTree->GenJetsAK8Pt[posGenJet];
	  WWTree->AK8_eta_gen = ReducedTree->GenJetsAK8Eta[posGenJet];
	  WWTree->AK8_phi_gen = ReducedTree->GenJetsAK8Phi[posGenJet];
	  WWTree->AK8_e_gen = ReducedTree->GenJetsAK8E[posGenJet];
	  WWTree->AK8_pruned_mass_gen = ReducedTree->GenJetsAK8_prunedMass[posGenJet];
	  WWTree->AK8_softdrop_mass_gen = ReducedTree->GenJetsAK8_softdropMass[posGenJet];
	  WWTree->AK8_softdrop_pt_gen = ReducedTree->GenJetsAK8_softdropPt[posGenJet];

          if (deltaR(ReducedTree->GenBosonEta[posWhad], ReducedTree->GenBosonPhi[posWhad],
                     WWTree->ungroomed_jet_eta, WWTree->ungroomed_jet_phi)<0.1)     ok++;
          total++;
	}

	deltaPhiOld=100.;
       	for (int i=0; i<ReducedTree->GenNuNum; i++) {
	  double deltaPhi = getDeltaPhi(ReducedTree->GenNuPhi[i],WWTree->v_phi);
	  if (abs(deltaPhi)>abs(deltaPhiOld))   continue;	  
	  temp.SetPtEtaPhiE(ReducedTree->GenNuPt[i],ReducedTree->GenNuEta[i],ReducedTree->GenNuPhi[i],ReducedTree->GenNuE[i]);
	  WWTree->nu_pz_gen=temp.Pz();	  
	  WWTree->nu_pt_gen=temp.Pt();	  
	  WWTree->nu_phi_gen=temp.Phi();	  
	  WWTree->nu_eta_gen=temp.Eta();
	  deltaPhiOld = deltaPhi;
	}		

	top_NNLO_weight[0] = 1.;
	top_NNLO_weight[1] = 1.;
        if ( TString(outputFile.c_str()).Contains("TTbar_powheg") && ReducedTree->GenTopNum > 1) { //look at here: http://arxiv.org/pdf/1511.00549.pdf
	  for (int i=0; i<2; i++) {
	    if (ReducedTree->GenTopPt[i]<60)                                           top_NNLO_weight[i] = 1./0.97;
	    else if (ReducedTree->GenTopPt[i]>=60 && ReducedTree->GenTopPt[i]<100)     top_NNLO_weight[i] = 1./0.995;
	    else if (ReducedTree->GenTopPt[i]>=100 && ReducedTree->GenTopPt[i]<150)    top_NNLO_weight[i] = 1.;
	    else if (ReducedTree->GenTopPt[i]>=150 && ReducedTree->GenTopPt[i]<200)    top_NNLO_weight[i] = 1./1.025;
	    else if (ReducedTree->GenTopPt[i]>=200 && ReducedTree->GenTopPt[i]<260)    top_NNLO_weight[i] = 1./1.045;
	    else if (ReducedTree->GenTopPt[i]>=260 && ReducedTree->GenTopPt[i]<320)    top_NNLO_weight[i] = 1./1.065;
	    else                                                                       top_NNLO_weight[i] = 1./1.095;
	  }
	}
	WWTree->top1_NNLO_Weight*=float(top_NNLO_weight[0]);
	WWTree->top2_NNLO_Weight*=float(top_NNLO_weight[1]);

	WWTree->gen_top1_pt = ReducedTree->GenTopPt[0];
	WWTree->gen_top2_pt = ReducedTree->GenTopPt[1];
      }
    
    WWTree->totalEventWeight = WWTree->genWeight*WWTree->eff_and_pu_Weight*WWTree->top1_NNLO_Weight*WWTree->top2_NNLO_Weight*WWTree->trig_eff_Weight;
    WWTree->totalEventWeight_2 = WWTree->genWeight*WWTree->eff_and_pu_Weight_2*WWTree->top1_NNLO_Weight*WWTree->top2_NNLO_Weight*WWTree->trig_eff_Weight;
    WWTree->totalEventWeight_3 = WWTree->genWeight*WWTree->eff_and_pu_Weight_3*WWTree->top1_NNLO_Weight*WWTree->top2_NNLO_Weight*WWTree->trig_eff_Weight;

    
    bool isBadEvent=false;
    if (isMC==0) {
      std::multimap<int,int>::iterator it = badEventsList.begin(); //filter bad events
      for (it=badEventsList.begin(); it!=badEventsList.end(); ++it) {
	if (it->first == WWTree->run && it->second == WWTree->event)
	  isBadEvent = true;
      }
    }
    if (isBadEvent)      continue;
    
    //fill the tree
    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"fill: "<<count<<std::endl; count++;
    outTree->Fill();
  }

  std::cout<<"matching: "<<(float)ok/(float)total<<std::endl;

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
