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
#include "TF1.h"
#include <TClonesArray.h>           

#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TAddJet.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"
#include "BaconAna/DataFormats/interface/TLHEWeight.hh"

// HEADER FILE FOR JES 
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// HEADER FOR B-TAG
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "../BtagUnc.hh"

#include "../interface/setOutputTree.h"
#include "../interface/METzCalculator.h"
#include "../interface/METzCalculator_Run2.h"
#include "../interface/analysisUtils.h"
#include "../interface/readJSONFile.h"
#include "../interface/Utils.hh"

using namespace std;

//*******MAIN*******************************************************************

int main (int argc, char** argv)
{ 
  
  int t0 = time(NULL);
  
  std::string inputFolder = argv[1];
  std::string outputFile = argv[2];
  int isMC = atoi(argv[3]);
  std::string cluster = argv[4];
  std::string inputTreeName = argv[5];
  std::string inputFile = argv[6];
  std::string xSecWeight = argv[7];
  std::string TotalNumberOfEntries = argv[8];
  std::string TotalNumberOfNegativeEntries = argv[9];
  int applyTrigger = atoi(argv[10]);
  std::string jsonFileName = argv[11];
  int isLocal = atoi(argv[12]);
  int VBFSel  = atoi(argv[13]);
  
  std::string leptonName;

  if ( VBFSel==1)	cout<<"==> VBF selection method : Select two highest pT jets"<<endl;
  else if ( VBFSel==2)	cout<<"==> VBF selection method : Select pair with highest mjj..."<<endl;
  else if ( VBFSel==3)	cout<<"==> VBF selection method : Select pair with highest DeltaEta..."<<endl;
  else {	cout<<"\n\nERROR:	Enter valid vbf selection criteria....\n\n"<<endl;
  		exit(0);  
	}
  
  std::string iHLTFile="${CMSSW_BASE}/src/BaconAna/DataFormats/data/HLTFile_25ns";
  const std::string cmssw_base = getenv("CMSSW_BASE");
  std::string cmssw_base_env = "${CMSSW_BASE}";
  size_t start_pos = iHLTFile.find(cmssw_base_env);
  if(start_pos != std::string::npos) {
  	iHLTFile.replace(start_pos, cmssw_base_env.length(), cmssw_base);
  }

  const baconhep::TTrigger triggerMenu(iHLTFile);  
  std::cout<<"Apply trigger: "<<applyTrigger<<std::endl;

//  TLorentzVector W_type0,W_type0_jes_up, W_type0_jes_dn, W_type0_jer_up, W_type0_jer_dn, W_type2, W_run2,W_puppi_type2, W_puppi_type0, W_puppi_run2;
  TLorentzVector LEP1, LEP2;// SJ1_PuppiAK8, SJ2_PuppiAK8, SJ1, SJ2;
  TLorentzVector NU0,NU1,NU2;//NU0_puppi,NU1_puppi,NU2_puppi;
  TLorentzVector NU0_jes_up, NU0_jes_dn;
 // TLorentzVector NU0_puppi_jes_up, NU0_puppi_jes_dn;
  TLorentzVector NU0_jer_up, NU0_jer_dn;
//  TLorentzVector JET, JET_PuppiAK8;
  TLorentzVector  AK4;
  TLorentzVector JET_jes_up, JET_jes_dn;// JET_PuppiAK8_jes_up, JET_PuppiAK8_jes_dn;
  TLorentzVector AK4_JET1,AK4_JET2;
  TLorentzVector AK4_JET1_jes_up, AK4_JET1_jes_dn;
  TLorentzVector AK4_JET2_jes_up, AK4_JET2_jes_dn;
 // TLorentzVector PuppiAK4_JET1,PuppiAK4_JET2;
//  TLorentzVector PuppiAK4_JET1_jes_up, PuppiAK4_JET1_jes_dn;
//  TLorentzVector PuppiAK4_JET2_jes_up, PuppiAK4_JET2_jes_dn;
  TLorentzVector VBF1,VBF2,TOT;
  TLorentzVector VBF1_jes_up, VBF1_jes_dn, VBF2_jes_up, VBF2_jes_dn;
  TLorentzVector ELE,MU;

  
  std::vector<TLorentzVector> tightMuon;
  std::vector<TLorentzVector> looseMuon;
  std::vector<TLorentzVector> tightEle;
  std::vector<TLorentzVector> looseEle;

  int ok=0, total=0;
  
  // Data structures to store info from TTrees
  baconhep::TEventInfo *info  	= new baconhep::TEventInfo();
  baconhep::TGenEventInfo *gen	= new baconhep::TGenEventInfo();
  TClonesArray *genPartArr 	= new TClonesArray("baconhep::TGenParticle");
  TClonesArray *muonArr    	= new TClonesArray("baconhep::TMuon");
  TClonesArray *electronArr	= new TClonesArray("baconhep::TElectron");
  TClonesArray *vertexArr	= new TClonesArray("baconhep::TVertex");
  TClonesArray *jetArr		= new TClonesArray("baconhep::TJet");
//  TClonesArray *vjetArrPuppi	= new TClonesArray("baconhep::TJet");
//  TClonesArray *vjetAddArrPuppi	= new TClonesArray("baconhep::TAddJet");
  TClonesArray *lheWgtArr	= new TClonesArray("baconhep::TLHEWeight");
  

  char command1[3000];
  char command2[3000];

  if ( cluster == "lxplus")
  	sprintf(command1, "eos find -f %s  | awk '!/log|fail/ {print $1}' | awk 'NF {print \"root://eoscms.cern.ch/\"$1}' > listTemp_%s.txt", (inputFolder).c_str(), outputFile.c_str());	// NF in awk command skips the blank line
  else 
	sprintf(command1,"xrdfs root://cmseos.fnal.gov ls %s | awk '{print \"root://cmseos.fnal.gov/\"$1}' > listTemp_%s.txt",(inputFolder).c_str(),  outputFile.c_str());
	//sprintf(command1,"eos root://cmseos.fnal.gov find -f %s | awk '!/log|fail/ {print $1}' | awk 'NF {print \"root://cmseos.fnal.gov/\"$1}' > listTemp_%s.txt",(inputFolder).c_str(),  outputFile.c_str());	// WORKS ONLY WITH INTERACTIVE NODE

  std::cout<<command1<<std::endl;
  sprintf(command2,"sed -i '/failed$/d' listTemp_OutPutRootFile.txt");
  system(command1);
  system(command2);
  char list1[2000];
  sprintf (list1, "listTemp_%s.txt", outputFile.c_str());
  ifstream rootList (list1);
  char command3[300];
  sprintf(command3, "rm listTemp_%s.txt", outputFile.c_str());
  system(command3);

  int fileCounter=0;

  vector<TString>  sampleName; 

  while (!rootList.eof())
  {
	char iRun_tW[700];
	rootList >> iRun_tW;
	if(!rootList.good())break;
	sampleName.push_back(iRun_tW);
	fileCounter++;
  }
TH1 *histo_aqgc[718];
  for(int j=0;j<718;j++)
    {
      stringstream ss;
      ss << j;
      string temp = ss.str();
      const char* name = temp.c_str();
      histo_aqgc[j] = new TH1D();
      histo_aqgc[j]->Sumw2();
}


  TFile *infile=0;
  TTree *eventTree=0;
  int cutEff[21]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  // Read and add pileup in histogram
  TFile* pileupFileMC = TFile::Open("puWeights_80x_37ifb.root");
  TH1D* puWeights = (TH1D*)pileupFileMC->Get("puWeights");
  TH1D* puWeightsUp = (TH1D*)pileupFileMC->Get("puWeightsUp");
  TH1D* puWeightsDown = (TH1D*)pileupFileMC->Get("puWeightsDown");
  puWeights->SetBins(75,0,75);
  puWeightsUp->SetBins(75,0,75);
  puWeightsDown->SetBins(75,0,75);

  //---------------- Root Files for ID, ISO, Trigger, GSF correctiosn for ELE and MU both: Starts -------------
  TFile* IDIsoEle = TFile::Open("egammaEffi_EGM2D_TightCutBasedIDSF.root","READ");
  TH1F *hIDIsoEle = (TH1F*)IDIsoEle->Get("EGamma_SF2D");

  TFile* GSFCorrEle = TFile::Open("egammaEffi_SF2D_GSF_tracking.root","READ");
  TH1F *hGSFCorrEle = (TH1F*)GSFCorrEle->Get("EGamma_SF2D");

  TFile* TriggerEle = TFile::Open("ElectronTrigger_SF.root","READ");
  TH1F* hTriggerEle = (TH1F*)TriggerEle->Get("HLT_Ele27");

  TFile* IDMuA = TFile::Open("MuonID_RunBCDEF_23SepReReco_19p72fb.root","READ");
  TH1F *hIDMuA = (TH1F*)IDMuA->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");

  TFile* IDMuB = TFile::Open("MuonID_RunGH_23SepReReco_16p146fb.root","READ");
  TH1F *hIDMuB = (TH1F*)IDMuB->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");

  TFile* IsoMuA = TFile::Open("MuonIso_RunBCDEF_23SepReReco_19p72fb.root","READ");
  TH1F *hIsoMuA = (TH1F*)IsoMuA->Get("TightISO_TightID_pt_eta/abseta_pt_ratio");

  TFile* IsoMuB = TFile::Open("MuonIso_RunGH_23SepReReco_16p146fb.root","READ");
  TH1F *hIsoMuB = (TH1F*)IsoMuB->Get("TightISO_TightID_pt_eta/abseta_pt_ratio");

  TFile* TriggerMuA = TFile::Open("MuonTrigger_RunBCDEF_23SepReReco_19p72fb.root","READ");
  TH1F* hTriggerMuA = (TH1F*)TriggerMuA->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");

  TFile* TriggerMuB = TFile::Open("MuonTrigger_RunGH_23SepReReco_16p146fb.root","READ");
  TH1F* hTriggerMuB = (TH1F*)TriggerMuB->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");

  //---------------- Root Files for ID, ISO, Trigger, GSF correctiosn for ELE and MU both: ENDS -------------
  
  TFile* file = TFile::Open( "puppiCorr.root","READ");
 // TF1* puppisd_corrGEN      = (TF1*)file->Get("puppiJECcorr_gen");
//  TF1* puppisd_corrRECO_cen = (TF1*)file->Get("puppiJECcorr_reco_0eta1v3");
//  TF1* puppisd_corrRECO_for = (TF1*)file->Get("puppiJECcorr_reco_1v3eta2v5");


  //---------output tree----------------
  TFile* outROOT = TFile::Open((outputFile+(".root")).c_str(),"recreate");
  outROOT->cd();
  TTree* outTree = new TTree("otree", "otree");
  setOutputTree* WWTree = new setOutputTree(outTree);

  int nEvents=0;	
  Long64_t jentry2=0;
  int count_genEvents=0;

  int nInputFiles = sampleName.size();

  if (isLocal==1) nInputFiles = 21;
  cout<<"==> Total number of input files : "<<nInputFiles<<endl;

  TH1D *MCpu = new TH1D("MCpu","",75,0,75);
  TH1D *MCpu_up = new TH1D("MCpu_up","",75,0,75);
  TH1D *MCpu_down = new TH1D("MCpu_down","",75,0,75);
  
  Long64_t TotalNumberOfEvents = 0, nNegEvents = 0; 


  // Set up b-tag scale factor readers
  BTagCalibration calib("csvv2", "CSVv2_Moriond17_B_H.csv");
  BTagCalibrationReader bTagReader(BTagEntry::OP_LOOSE,  // working point: can be OP_LOOSE, OP_MEDIUM, OP_TIGHT 
                                   "central",             // label for the central value (see the scale factor file)
                                   {"up","down"});        // vector of labels for systematics
  bTagReader.load(calib, BTagEntry::FLAV_B, "comb");      // use the "comb" measurements for b-jets
  bTagReader.load(calib, BTagEntry::FLAV_C, "comb");      // use the "comb" measurements for c-jets
  bTagReader.load(calib, BTagEntry::FLAV_UDSG, "incl");   // use the "incl" measurements for light jets
  

  vector<const baconhep::TJet*> goodJetsv;
  
  // Loop on input files
  for(int i=0;i<nInputFiles;i++)
  {//loop10 begins
     infile = TFile::Open(sampleName[i]);
     eventTree = (TTree*)infile->Get("Events");
     
     TotalNumberOfEvents+=eventTree->GetEntries();
     if(isMC)
     { 
        eventTree->SetBranchAddress("Info", &info);    TBranch *infoBr = eventTree->GetBranch("Info");
  	TBranch *genBr=0;
     	eventTree->SetBranchAddress("GenEvtInfo", &gen); genBr = eventTree->GetBranch("GenEvtInfo");
	for (Long64_t jentry=0; jentry<eventTree->GetEntries();jentry++,jentry2++)
	{
	  //eventTree->GetEntry(jentry);
	    genBr->GetEntry(jentry);
    	    infoBr->GetEntry(jentry);	    
	    MCpu->Fill(info->nPUmean);
	    MCpu_up->Fill(info->nPUmeanp);
	    MCpu_down->Fill(info->nPUmeanm);
	    if (jentry2%50000 == 0) std::cout << "\t File no. " << i << "; Neg Event Count; read entry: " << jentry2 <<"/"<<TotalNumberOfEvents<<std:: endl;
	    if (gen->weight<0)	nNegEvents++;
	}
     }
     delete infile;
     infile=0, eventTree=0;
  }  ///loop10 ends
  
  
  cout<<"==> Total number of events : "<<TotalNumberOfEvents<<endl;
  cout<<"==> Total number of negative events : "<<nNegEvents<<endl;

  float weight = std::atof(xSecWeight.c_str())/(std::atof(TotalNumberOfEntries.c_str()) - 2*std::atof(TotalNumberOfNegativeEntries.c_str()));
  cout<<"Weight of cross-sec/events = "<<weight<<endl;
  int totalEntries=0;


  JetCorrectorParameters paramAK4puppi("Summer16_23Sep2016V4_MC_JEC/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFPuppi.txt");
  JetCorrectorParameters paramAK4chs("Summer16_23Sep2016V4_MC_JEC/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt");
  JetCorrectorParameters paramAK8chs("Summer16_23Sep2016V4_MC_JEC/Summer16_23Sep2016V4_MC_Uncertainty_AK8PFchs.txt");
  JetCorrectorParameters paramAK8puppi("Summer16_23Sep2016V4_MC_JEC/Summer16_23Sep2016V4_MC_Uncertainty_AK8PFPuppi.txt");
  JetCorrectionUncertainty *fJetUnc_AK4chs = new JetCorrectionUncertainty(paramAK4chs);
//  JetCorrectionUncertainty *fJetUnc_AK4puppi = new JetCorrectionUncertainty(paramAK4puppi);
 // JetCorrectionUncertainty *fJetUnc_AK8chs = new JetCorrectionUncertainty(paramAK8chs);
//  JetCorrectionUncertainty *fJetUnc_AK8puppi = new JetCorrectionUncertainty(paramAK8puppi);

  //---------start loop on events------------
  std::cout << "---------start loop on events------------" << std::endl;
  jentry2=0;
  for(int i=0;i<nInputFiles;i++)
  {
  cout<<"\n\n=====	Processing File Number : "<<i<<"/"<<nInputFiles<<"\n\t"<<sampleName[i]<<"\n-------"<<endl;

  infile = TFile::Open(sampleName[i]);
  eventTree = (TTree*)infile->Get("Events");
  
  totalEntries+=eventTree->GetEntries();

  nEvents=eventTree->GetEntries();

  cout<<"\t==> Entries = "<<nEvents<<endl;



  eventTree->SetBranchAddress("Info", &info);    TBranch *infoBr = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Muon", &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("PV",   &vertexArr); TBranch *vertexBr = eventTree->GetBranch("PV");
  eventTree->SetBranchAddress("AK4CHS",   &jetArr); TBranch *jetBr = eventTree->GetBranch("AK4CHS");    
 // eventTree->SetBranchAddress("AK8Puppi",   &vjetArrPuppi); TBranch *vjetBrPuppi = eventTree->GetBranch("AK8Puppi");  
 // eventTree->SetBranchAddress("AddAK8Puppi",   &vjetAddArrPuppi); TBranch *vjetAddBrPuppi = eventTree->GetBranch("AddAK8Puppi");  
  TBranch *genBr=0, *genPartBr=0, *lhePartBr=0;
  if(isMC)
     {  //loop11 begins 
       eventTree->SetBranchAddress("GenEvtInfo", &gen); genBr = eventTree->GetBranch("GenEvtInfo");
       eventTree->SetBranchAddress("GenParticle",&genPartArr); genPartBr = eventTree->GetBranch("GenParticle");
       if(eventTree->GetListOfBranches()->FindObject("LHEWeight"))
       {
       eventTree->SetBranchAddress("LHEWeight",&lheWgtArr); lhePartBr = eventTree->GetBranch("LHEWeight");	       }
     } //loop11 ends

  //for (Long64_t jentry=0; jentry<172; jentry++,jentry2++)
  for (Long64_t jentry=0; jentry<eventTree->GetEntries();jentry++,jentry2++)
  {
   // if (jentry2 != 87 && jentry2 != 113) continue;	// For debug
//if (jentry2>100) continue;
    infoBr->GetEntry(jentry);	    


    int GenPassCut = 0;

    tightMuon.clear();
    tightEle.clear();
    looseMuon.clear();
    looseEle.clear();
    
//cout<<"check1--all good"<<endl;
    if (jentry2%10000 == 0) std::cout << "\tread entry: " << jentry2 <<"/"<<TotalNumberOfEvents<<std:: endl;
    
    //*********************************
    WWTree->initializeVariables(); //initialize all variables

    WWTree->run   = info->runNum;
    WWTree->event = info->evtNum;
    WWTree->lumi  = info->lumiSec;


    /////////////////MC Info
    if (isMC==1)
    {//loop14 begins
      lheWgtArr->Clear();
      if(lhePartBr)
	{
	  lhePartBr->GetEntry(jentry);
	}
      genPartArr->Clear();
      genBr->GetEntry(jentry);
      genPartBr->GetEntry(jentry);
      TLorentzVector hadW, lepW, VBFJ1, VBFJ2, VBFJ, temp;
      TLorentzVector genLep, genNeutrino, genWquarks, genVBFquarks;
      std::vector<TLorentzVector> v_GEN_hadW, v_GEN_lepW, v_GEN_VBFJ1, v_GEN_VBFJ2, v_GEN_VBFJ, v_GEN_temp;
      std::vector<TLorentzVector> v_genLep, v_genNeutrino, v_genWquarks, v_genVBFquarks;

      v_GEN_hadW.clear();	v_GEN_lepW.clear();	v_GEN_VBFJ1.clear();	v_GEN_VBFJ2.clear();
      v_GEN_VBFJ.clear();	v_GEN_temp.clear();	v_genLep.clear();	v_genNeutrino.clear();
      v_genWquarks.clear();//	v_genVBFquarks.clear();
 //  cout<<"check2--all good"<<endl;
   
     for (int i = 0; i<genPartArr->GetEntries();i++)
	{  //loop12 begins
	  const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*) ((*genPartArr)[i]);
	  Int_t parentPdg=dynamic_cast<baconhep::TGenParticle *>(genPartArr->At(genloop->parent>-1 ? genloop->parent : 0))->pdgId;
	  if( (abs(genloop->pdgId) == 11 || abs(genloop->pdgId) == 13 ) && abs(parentPdg) == 24)
	    {
	      genLep.SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	      v_genLep.push_back(genLep);
	    }
	  if( (abs(genloop->pdgId) == 12 || abs(genloop->pdgId) == 14 ) && abs(parentPdg) == 24)
	    {
	      genNeutrino.SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	      v_genNeutrino.push_back(genNeutrino);
	    }
	
	  if( (abs(genloop->pdgId) == 1 || abs(genloop->pdgId) == 3 || abs(genloop->pdgId) == 5 || abs(genloop->pdgId) == 2 || abs(genloop->pdgId) == 4 || abs(genloop->pdgId) == 6) && genloop->status == 23 && abs(parentPdg) != 24)
	    {
	      genVBFquarks.SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
	      v_genVBFquarks.push_back(genVBFquarks);
	    }
	} //loop12 ends 
      if (v_genLep.size()==2 && v_genNeutrino.size()==2 && v_genVBFquarks.size()==2)
	{ //loop13 begins
	  WWTree->isGen           = 1;
	  WWTree->lep1_pt_gen      = v_genLep[0].Pt();
	  WWTree->lep1_eta_gen     = v_genLep[0].Eta();
	  WWTree->lep1_m_gen       = (v_genLep[0]).M();
          WWTree->nu1_pz_gen	  = v_genNeutrino[0].Pz();
	  WWTree->nu1_pt_gen	  = v_genNeutrino[0].Pt();
	  WWTree->nu1_eta_gen	  = v_genNeutrino[0].Eta();
          WWTree->lep2_pt_gen      = (v_genLep[1]).Pt();
	  WWTree->lep2_eta_gen     = (v_genLep[1]).Eta();
          WWTree->lep2_m_gen     = (v_genLep[1]).M();
	  WWTree->nu2_pz_gen	  = (v_genNeutrino[1]).Pz();
	  WWTree->nu2_pt_gen	  = (v_genNeutrino[1]).Pt();
	  WWTree->nu2_eta_gen	  = (v_genNeutrino[1]).Eta();


          WWTree->dilep_m_gen    = (v_genLep[0]+v_genLep[1]).M();    

	  
	  WWTree->lepW1_pt_gen     = (v_genLep[0]+v_genNeutrino[0]).Pt();
	  WWTree->lepW1_eta_gen    = (v_genLep[0]+v_genNeutrino[0]).Eta();
	  WWTree->lepW1_phi_gen    = (v_genLep[0]+v_genNeutrino[0]).Phi();
	  WWTree->lepW1_e_gen      = (v_genLep[0]+v_genNeutrino[0]).E();
	  WWTree->lepW1_m_gen      = (v_genLep[0]+v_genNeutrino[0]).M();

          WWTree->lepW2_pt_gen     = (v_genLep[1]+v_genNeutrino[1]).Pt();
	  WWTree->lepW2_eta_gen    = (v_genLep[1]+v_genNeutrino[1]).Eta();
	  WWTree->lepW2_phi_gen    = (v_genLep[1]+v_genNeutrino[1]).Phi();
	  WWTree->lepW2_e_gen      = (v_genLep[1]+v_genNeutrino[1]).E();
	  WWTree->lepW2_m_gen      = (v_genLep[1]+v_genNeutrino[1]).M();
        
          WWTree->WW_mass_gen	= (v_genLep[0] + v_genNeutrino[0] + v_genLep[1] + v_genNeutrino[1]).M();
	  WWTree->WW_mT_gen	= (v_genLep[0] + v_genNeutrino[0] + v_genLep[1] + v_genNeutrino[1]).Mt();
	  WWTree->WW_pT_gen	= (v_genLep[0] + v_genNeutrino[0] + v_genLep[1] + v_genNeutrino[1]).Pt();


	  WWTree->AK4_1_pt_gen	= v_genVBFquarks[0].Pt();
	  WWTree->AK4_1_eta_gen	= v_genVBFquarks[0].Eta();
	  WWTree->AK4_1_phi_gen	= v_genVBFquarks[0].Phi();
	  WWTree->AK4_1_e_gen	= v_genVBFquarks[0].E();
	  WWTree->AK4_1_mass_gen	= v_genVBFquarks[0].M();
	  
	  WWTree->AK4_2_pt_gen	= v_genVBFquarks[1].Pt();
	  WWTree->AK4_2_eta_gen	= v_genVBFquarks[1].Eta();
	  WWTree->AK4_2_phi_gen	= v_genVBFquarks[1].Phi();
	  WWTree->AK4_2_e_gen	= v_genVBFquarks[1].E();
	  WWTree->AK4_2_mass_gen	= v_genVBFquarks[1].M();
	  
	  WWTree->AK4_jj_DeltaEta_gen = abs(v_genVBFquarks[0].Eta() - v_genVBFquarks[1].Eta());
	  WWTree->AK4_jj_mass_gen = (v_genVBFquarks[0] + v_genVBFquarks[1]).M();
	 
         WWTree->Zeppen1=abs(WWTree->lep1_eta_gen -(WWTree->AK4_1_eta_gen + WWTree->AK4_2_eta_gen)/2)/WWTree->AK4_jj_DeltaEta_gen ;
         WWTree->Zeppen2=abs(WWTree->lep2_eta_gen -(WWTree->AK4_1_eta_gen + WWTree->AK4_2_eta_gen)/2)/WWTree->AK4_jj_DeltaEta_gen ;




 
	  count_genEvents++;
	}//loop13 ends
//cout<<"check3--all good"<<endl;

	if ( WWTree->lep1_pt_gen > 25 && WWTree->lep2_pt_gen >20 && WWTree->dilep_m_gen >20   && abs(WWTree->AK4_jj_DeltaEta_gen)>2.5 && WWTree->AK4_1_pt_gen >30 && WWTree->AK4_2_pt_gen>30 && WWTree->AK4_jj_mass_gen > 500   && abs((WWTree->dilep_m_gen) -91.2)>15 && max(WWTree->Zeppen1,WWTree->Zeppen2)<0.75 ) 
         {
	        GenPassCut = 1;
	 }
    	for (int i = 0; i<lheWgtArr->GetEntries();i++)     // Note that i is starting from 446.
	{
		const baconhep::TLHEWeight *lhe = (baconhep::TLHEWeight*)((*lheWgtArr)[i]);
		//WWTree->LHEid[i] = lhe->id;
		WWTree->LHEWeight[i] = lhe->weight;
	}
    }//loop14 ends

    
    WWTree->issignal = 0;
    WWTree->wSampleWeight = weight; //xsec/TotalNumberOfEntries
    WWTree->pu_Weight = 1.; //temporary value
    WWTree->pu_Weight_up = 1.; //temporary value
    WWTree->pu_Weight_down = 1.; //temporary value
    WWTree->top1_NNLO_Weight = 1.;
    WWTree->top2_NNLO_Weight = 1.;
    WWTree->id_eff_Weight = 1.;
    WWTree->id_eff_Weight2 = 1.;
    WWTree->trig_eff_Weight = 1.;
    WWTree->trig_eff_Weight2 = 1.;

    if (gen->weight>0)
      WWTree->genWeight=1.;
    else if (gen->weight<0) {  //loop15 begins
      WWTree->genWeight=-1.;
      //nNegEvents++;
    }  //loop15 ends
    cutEff[0]++;

    if (isMC==1)
    {   //loop16 begins
    if (GenPassCut == 1)   cutEff[1]++;
    } //loop16ends

    
    vertexArr->Clear();
    vertexBr->GetEntry(jentry);
    WWTree->nPV = vertexArr->GetEntries();
//  cout<<"check4--all good"<<endl;

    //PILE-UP WEIGHT
    if (isMC==1) { //loop17 begins
       if(int(info->nPUmean)<75){
           WWTree->pu_Weight = puWeights->GetBinContent(info->nPUmean); //our pu recipe
           WWTree->pu_Weight_up = puWeightsUp->GetBinContent(info->nPUmean); //our pu recipe
           WWTree->pu_Weight_down = puWeightsDown->GetBinContent(info->nPUmean); //our pu recipe
       }
       else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
           std::cout<<"Warning! n_pu too big"<<std::endl;
	   // throw logic_error("n_pu too big");
	   WWTree->pu_Weight = 0.;
	   WWTree->pu_Weight_up = 0.;
	   WWTree->pu_Weight_down = 0.;
       } 
    }//loop17 ends

//PILE-UP WEIGHT END

   if(applyTrigger==1)
     if(!(triggerMenu.pass("HLT_IsoMu24_v*",info->triggerBits) || triggerMenu.pass("HLT_IsoTkMu24_v*",info->triggerBits) ||  triggerMenu.pass("HLT_Ele27_WPTight_Gsf_v*",info->triggerBits))) continue;
  
    /////////////////THE SELECTED LEPTON
    int nTightEle=0, nLooseEle=0;
    int nTightMu=0, nLooseMu=0;
    double pt_cut = 20;//pvrevious=25
    double leadelept_cut = 25;
    double leadmupt_cut = 25;
//  double sublead
//electron part begins
    electronArr->Clear();
    electronBr->GetEntry(jentry);
    const baconhep::TElectron *leadele = NULL;
    const baconhep::TElectron *subele = NULL;
    double leadeleE=-999, subeleE=-999;
    double iso = 1.5;
    for (int i=0; i<electronArr->GetEntries(); i++) {    //loop1 on electron begins
      const baconhep::TElectron *ele = (baconhep::TElectron*)((*electronArr)[i]);
      if (ele->pt<=pt_cut) continue;
      if (fabs(ele->eta)>=2.5) continue;
      if(!passEleLooseSel(ele,info->rhoIso)) continue;
      nLooseEle++;
      if(!passEleTightSel(ele,info->rhoIso)) continue;
      ELE.SetPtEtaPhiM(ele->pt,ele->eta,ele->phi,0.0005109989461);
      tightEle.push_back(ELE);
      nTightEle++;
      iso = ele->chHadIso + TMath::Max( 0.0,(ele->gammaIso + ele->neuHadIso - info->rhoIso*eleEffArea(ele->eta)) );
      if(!leadele || ele->pt>leadele->pt)
	{
	  if(!(ele->pt>leadelept_cut)) continue;
	  subele = leadele;
	  leadele = ele;
	  leadeleE = ELE.E();
	  WWTree->l_iso1 = iso/ele->pt;
	}
      else if (!subele || ele->pt > subele->pt)
	{
	  subele = ele;
	  subeleE = ELE.E();
	  WWTree->l_iso2 = iso/ele->pt;
	}
    } //loop1 on electron ends
  

//muon part begins
    muonArr->Clear();
    muonBr->GetEntry(jentry);
    const baconhep::TMuon *leadmu = NULL;
    const baconhep::TMuon *submu = NULL;
    double leadmue=-999, submue = -999;
    iso = 1.5;
    for(Int_t i=0; i<muonArr->GetEntries(); i++) { //loop on muon begins
      const baconhep::TMuon *mu = (baconhep::TMuon*)((*muonArr)[i]);
      if (mu->pt<pt_cut) continue;
      if (fabs(mu->eta)>=2.4) continue;
      if(!passMuonLooseSel(mu)) continue;
      nLooseMu++;
      if(!passMuonTightSel(mu)) continue;
      nTightMu++;
      MU.SetPtEtaPhiM(mu->pt,mu->eta,mu->phi,0.1056583745);
      tightMuon.push_back(MU);
      iso = mu->chHadIso + TMath::Max(mu->neuHadIso + mu->gammaIso - 0.5*(mu->puIso), double(0));
      if(!leadmu || mu->pt>leadmu->pt)
	{
	  if(!(mu->pt>leadmupt_cut)) continue;
	  submu = leadmu;
	  leadmu = mu;
	  leadmue = MU.E();
	  WWTree->l_iso1 = iso/mu->pt;
	}
      else if (!submu || mu->pt > submu->pt)
	{
	  submu = mu;
	  submue = MU.E();
	  WWTree->l_iso2 = iso/mu->pt;
	}
    }   //loop on muon end

//if (!(nTightMu+nTightEle)==2) continue;
//cout<<"Tight muons="<<nTightMu<<"Tight Electrons="<<nTightEle<<endl;
if((nTightMu+nTightEle)==2)
{
if (nTightMu==2)
{



    if(leadmu)
      {
	WWTree->l_pt1  = leadmu->pt;
	WWTree->l_eta1 = leadmu->eta;
	WWTree->l_phi1 = leadmu->phi;	
	WWTree->l_e1 = leadmue;	
	WWTree->l_charge1 = leadmu->q;
      }
    if(submu)
      {
	WWTree->l_pt2  = submu->pt;
	WWTree->l_eta2 = submu->eta;
	WWTree->l_phi2 = submu->phi;	
	WWTree->l_e2 = submue;	
	WWTree->l_charge2 = submu->q;
      }
}
 else if (nTightEle==2)
{   if(leadele)
      {
        WWTree->l_pt1  = leadele->pt;
        WWTree->l_eta1 = leadele->eta;
        WWTree->l_phi1 = leadele->phi;
        WWTree->l_e1 = leadeleE;
        WWTree->l_charge1 = leadele->q;
      }
    if(subele)
      {
        WWTree->l_pt2  = subele->pt;
        WWTree->l_eta2 = subele->eta;
        WWTree->l_phi2 = subele->phi;
        WWTree->l_e2 = subeleE;
        WWTree->l_charge2 = subele->q;
      }   //electron part ends

}

 else if (nTightEle==1 && nTightMu==1)
{
 if(leadmu)
{
        WWTree->l_pt1  = leadmu->pt;
        WWTree->l_eta1 = leadmu->eta;
        WWTree->l_phi1 = leadmu->phi;
        WWTree->l_e1 = leadmue;
        WWTree->l_charge1 = leadmu->q;

}
if (leadele)

{
        WWTree->l_pt2  = leadele->pt;
        WWTree->l_eta2 = leadele->eta;
        WWTree->l_phi2 = leadele->phi;
        WWTree->l_e2 = leadeleE;
        WWTree->l_charge2 = leadele->q;
}


}
}


if (WWTree->l_charge1!=WWTree->l_charge2) continue;
 if(!(WWTree->l_pt1>0)&&(WWTree->l_pt2>0)) continue;
    if ((nTightMu+nTightEle)<2) continue;
    if((nLooseEle+nLooseMu)>2) continue;
 //   if(nTightMu>1 && nLooseEle>1) continue;                                                                                                 
 //   if(nTightEle>1 && nLooseMu>1) continue;                                                                                                   
 //   if(nTightMu==1 && nLooseMu>1) continue;
 //   if(nTightEle==1 && nLooseEle>1) continue;
  //

//    if(!(WWTree->l_pt1>0)&&(WWTree->l_pt2>0)) continue;
  
//cout<<"check51--all good"<<endl;
cutEff[2]++;

 /*   if(nTightMu>0){
      WWTree->type=0;
      leptonName = "mu";	
    }else{
      WWTree->type=1;
      leptonName = "el";
    }*/
 //cout<<"check52--all good"<<endl;  
    LEP1.SetPtEtaPhiE(WWTree->l_pt1,WWTree->l_eta1,WWTree->l_phi1,WWTree->l_e1);
    LEP2.SetPtEtaPhiE(WWTree->l_pt2,WWTree->l_eta2,WWTree->l_phi2,WWTree->l_e2);
   if(WWTree->l_pt1>0&&WWTree->l_pt2>0)
     {
	WWTree->dilep_pt  = (LEP1+LEP2).Pt();
	WWTree->dilep_eta = (LEP1+LEP2).Eta();
	WWTree->dilep_phi = (LEP1+LEP2).Phi();	
	WWTree->dilep_m = (LEP1+LEP2).M();
}
//cout<<"check53--all good"<<endl;
for (int j=0;j<718;j++)
                    {
                      histo_aqgc[j]->Fill(dilep_m,(LHEWeight[j+446]));
}

if(WWTree->dilep_m <20) continue;
cutEff[3]++;

float Zmass = 91.2;
//cout<<"check54--all good"<<endl;
if((nTightEle==2) && (abs((WWTree->dilep_m) - Zmass) < 15))	continue;
//cout<<"check55--all good"<<endl;
cutEff[4]++;
//cout<<"check6--all good"<<endl;

     // }

    //ID & GSF efficiency SF for electrons (https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Electron_efficiencies_and_scale)
    if (strcmp(leptonName.c_str(),"el")==0 && isMC==1) {//loop 3 begins
	//  apply ID, ISO SF's
	WWTree->id_eff_Weight = GetSFs_Lepton(WWTree->l_pt1, WWTree->l_eta1, hIDIsoEle);	// Get Scale factor corresponding to the pt and eta.
	
	// apply GSF/RECO SF's for electrons
	WWTree->id_eff_Weight = WWTree->id_eff_Weight*GetSFs_Lepton(WWTree->l_pt1, WWTree->l_eta1, hGSFCorrEle);
	WWTree->trig_eff_Weight = GetSFs_Lepton(WWTree->l_pt1, WWTree->l_eta1, hTriggerEle);
    	
//	if(WWTree->l_pt2>0){
	   //  apply ID, ISO SF's
	   WWTree->id_eff_Weight2 = GetSFs_Lepton(WWTree->l_pt2, WWTree->l_eta2, hIDIsoEle);	// Get Scale factor corresponding to the pt and eta.

	   // apply GSF/RECO SF's for electrons
	   WWTree->id_eff_Weight2 = WWTree->id_eff_Weight2*GetSFs_Lepton(WWTree->l_pt2, WWTree->l_eta2, hGSFCorrEle);
	   WWTree->trig_eff_Weight2 = GetSFs_Lepton(WWTree->l_pt2, WWTree->l_eta2, hTriggerEle);
//	}
    }//loop 3 ends

    //ID&ISO efficiency SF for muons (https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults)
    if (strcmp(leptonName.c_str(),"mu")==0 && isMC==1) { //loop4 begins
	//  apply ID SF's
	if (WWTree->run<278820){
		WWTree->id_eff_Weight = GetSFs_Lepton(WWTree->l_pt1, abs(WWTree->l_eta1), hIDMuA);
		if (WWTree->l_pt2>0) WWTree->id_eff_Weight2 = GetSFs_Lepton(WWTree->l_pt2, abs(WWTree->l_eta2), hIDMuA);}
	else{
		WWTree->id_eff_Weight = GetSFs_Lepton(WWTree->l_pt1, abs(WWTree->l_eta1), hIDMuB);
		if (WWTree->l_pt2>0) WWTree->id_eff_Weight2 = GetSFs_Lepton(WWTree->l_pt2, abs(WWTree->l_eta2), hIDMuB);}

	//  apply ISO SF's
	if (WWTree->run<278820){
		WWTree->id_eff_Weight = WWTree->id_eff_Weight*GetSFs_Lepton(WWTree->l_pt1, abs(WWTree->l_eta1), hIsoMuA);
		WWTree->id_eff_Weight2 = WWTree->id_eff_Weight2*GetSFs_Lepton(WWTree->l_pt2, abs(WWTree->l_eta2), hIsoMuA);}
	else{
		WWTree->id_eff_Weight = WWTree->id_eff_Weight*GetSFs_Lepton(WWTree->l_pt1, abs(WWTree->l_eta1), hIsoMuB);
		WWTree->id_eff_Weight2 = WWTree->id_eff_Weight2*GetSFs_Lepton(WWTree->l_pt2, abs(WWTree->l_eta2), hIsoMuB);}
	
	// apply Trigger SF's
	if (WWTree->run<278820){
		WWTree->trig_eff_Weight  = GetSFs_Lepton(WWTree->l_pt1, abs(WWTree->l_eta1), hTriggerMuA);
		WWTree->trig_eff_Weight2 = GetSFs_Lepton(WWTree->l_pt2, abs(WWTree->l_eta2), hTriggerMuA);}
	else{
		WWTree->trig_eff_Weight  = GetSFs_Lepton(WWTree->l_pt1, abs(WWTree->l_eta1), hTriggerMuB);
 	WWTree->trig_eff_Weight2 = GetSFs_Lepton(WWTree->l_pt2, abs(WWTree->l_eta2), hTriggerMuB);}
    } //loop4 ends
	
    //////////////THE MET
    
    // //preselection on met
     if (info->pfMETC < 40) continue;   //Et(miss)>40GeV
    // if ((WWTree->l_pt1>0)&&(WWTree->l_pt2>0))
    cutEff[5]++;
  // cout<<"check7--all good"<<endl;

    
    TLorentzVector W_Met;
    WWTree->pfMET_Corr = info->pfMETC;
    WWTree->pfMET_Corr_phi = info->pfMETCphi;

  
    TLorentzVector W_Met_jes_up, W_Met_jes_dn, W_Met_jer_up, W_Met_jer_dn, AK4Up, AK4Down;//, AK4Up_Puppi, AK4Down_Puppi;
  
    W_Met.SetPxPyPzE(info->pfMETC * TMath::Cos(info->pfMETCphi), info->pfMETC * TMath::Sin(info->pfMETCphi), 0., sqrt(info->pfMETC*info->pfMETC));
    W_Met_jer_up.SetPxPyPzE(info->pfMETCjerup * TMath::Cos(info->pfMETCphijerup), info->pfMETCjerup * TMath::Sin(info->pfMETCphijerup), 0., sqrt(info->pfMETCjerup*info->pfMETCjerup) );
    W_Met_jer_dn.SetPxPyPzE(info->pfMETCjerdn * TMath::Cos(info->pfMETCphijerdn), info->pfMETCjerdn * TMath::Sin(info->pfMETCphijerdn), 0., sqrt(info->pfMETCjerdn*info->pfMETCjerdn) );

    ////////////////////////////////////////////////////////////////
    //		
    //		MET JES Calculate
    //
    ////////////////////////////////////////////////////////////////
    jetArr->Clear();
    jetBr->GetEntry(jentry);
    for ( int i=0; i<jetArr->GetEntries(); i++) //loop on AK4 jet
    { //loop5 on AK4 jets begins
      const baconhep::TJet *jet = (baconhep::TJet*)((*jetArr)[i]);
      TLorentzVector AK4_LV_temp, AK4_LV_temp2;

      // Get uncertanity
      double unc = func(jet->pt, jet->eta, fJetUnc_AK4chs); 

      // Get AK4 LorentzVector 
      AK4_LV_temp.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);

      // calculate Up variation
      AK4_LV_temp2.SetPtEtaPhiM((1.+unc)*jet->pt, jet->eta, jet->phi, (1.+unc)*jet->mass);
      AK4Up += AK4_LV_temp2 - AK4_LV_temp;

      // calculate Down variation
      AK4_LV_temp2.SetPtEtaPhiM((1.-unc)*jet->pt, jet->eta, jet->phi, (1.-unc)*jet->mass);
      AK4Down += AK4_LV_temp2 - AK4_LV_temp;
    }  //loop5 on AK4 jets end
    W_Met_jes_up = W_Met + AK4Up;
    W_Met_jes_dn = W_Met + AK4Down;
    //////////////////////////////////////// END: MET JES Calculate
   
    /////////VBF and b-tag section
    WWTree->njets=0;
    WWTree->nBTagJet_loose=0;
    WWTree->nBTagJet_medium=0;
    WWTree->nBTagJet_tight=0;
    
    WWTree->njets_unmerged=0;
    WWTree->nBTagJet_loose_unmerged=0;
    WWTree->nBTagJet_medium_unmerged=0;
    WWTree->nBTagJet_tight_unmerged=0;

    int OnlyTwoVBFTypeJets = 0;
    
    std::vector<int> indexGoodVBFJets;


    jetArr->Clear();
    jetBr->GetEntry(jentry);
    for ( int i=0; i<jetArr->GetEntries(); i++) //loop on AK4 jet
    {//loop r begins
      const baconhep::TJet *jet = (baconhep::TJet*)((*jetArr)[i]);
      bool isCleaned = true;
   //   bool isCleanedFromFatJet = true;
//28feb        double unc = func(jet->pt, jet->eta, fJetUnc_AK4chs); 
	// calculate Up variation
//28feb	VBF1_jes_up.SetPtEtaPhiM((1.+unc)*jet->pt, jet->eta, jet->phi, (1.+unc)*jet->mass);	
	// calculate Down variation
//28feb	VBF1_jes_dn.SetPtEtaPhiM((1.-unc)*jet->pt, jet->eta, jet->phi, (1.-unc)*jet->mass);	
      
      if (jet->pt<=30) continue;
      if (!passJetLooseSel(jet)) continue;

      //fill B-Tag info
      
      if (fabs(jet->eta) < 2.4 && jet->pt>30){
      		if (jet->csv>0.5426)  WWTree->nBTagJet_loose++;
      		if (jet->csv>0.8484)  WWTree->nBTagJet_medium++;
      		if (jet->csv>0.9535)  WWTree->nBTagJet_tight++;

		goodJetsv.push_back(jet);
      }
    //if (WWTree->nBTagJet_medium>0) continue;
    
//  if (jet->csv>0.8484)   continue;
      //CLEANING FROM LEPTONS
      for ( std::size_t j=0; j<tightEle.size(); j++) {
        if (deltaR(tightEle.at(j).Eta(), tightEle.at(j).Phi(),
                   jet->eta,   jet->phi) < 0.3) {
          isCleaned = false;
        }
      }
      for ( std::size_t j=0; j<tightMuon.size(); j++) {
        if (deltaR(tightMuon.at(j).Eta(), tightMuon.at(j).Phi(),
                   jet->eta,   jet->phi) < 0.3) {
          isCleaned = false;
        }
      }
      
      if (isCleaned==false) continue;
  //    if (isCleanedFromFatJet==false) continue;
       if (fabs(jet->eta)>=5.0) continue;//2.4 check change
      if (jet->pt<=30) continue;
//if (jet->csv>0.8484) continue;
         indexGoodVBFJets.push_back(i); //save index of the "good" vbf jets candidates
      
//      if (fabs(jet->eta)>=5.0) continue;//2.4 check change
      
      WWTree->njets++;
      AK4.SetPtEtaPhiM(jet->pt,jet->eta,jet->phi,jet->mass);
      
      
      
      //------------------------------
      // !!! VBF non-Puppi missing !!!
      //------------------------------
    }//loop r ends
//	cout<<"check8--all good"<<indexGoodVBFJets.size()<<endl;
//=========================================================================================================================================================================================================================================================================================================================================================================
//	if (WWTree->nBTagJet_medium>0) continue;
if (indexGoodVBFJets.size()<2) continue;
if (WWTree->nBTagJet_medium>=1) continue;
cutEff[6]++;
    if (indexGoodVBFJets.size()>=2) 
    {//loop y begins
       
      
      float tempPtMax=0.;
      float DRvbf;
      int nVBF1=-1, nVBF2=-1; //position of the two vbf jets
      
      for (std::size_t i=0; i<indexGoodVBFJets.size()-1; i++) {//loop p
        for ( std::size_t ii=i+1; ii<indexGoodVBFJets.size(); ii++) {//loop q 
	  const baconhep::TJet *jet1 = (baconhep::TJet*)((*jetArr)[indexGoodVBFJets.at(i)]);
	  const baconhep::TJet *jet2 = (baconhep::TJet*)((*jetArr)[indexGoodVBFJets.at(ii)]);
	  if (jet1->pt < 30) continue;
	  if (jet2->pt < 30) continue;
          VBF1.SetPtEtaPhiM(jet1->pt,jet1->eta,jet1->phi,jet1->mass);
          VBF2.SetPtEtaPhiM(jet2->pt,jet2->eta,jet2->phi,jet2->mass);
          TOT = VBF1 + VBF2;

	  // calculate Up/Down variation for leading jet
     //28feb     double unc1 = func(jet1->pt, jet1->eta, fJetUnc_AK4chs); 
//28feb	  VBF1_jes_up.SetPtEtaPhiM((1.+unc1)*jet1->pt, jet1->eta, jet1->phi, (1.+unc1)*jet1->mass);	
//28feb	  VBF1_jes_dn.SetPtEtaPhiM((1.-unc1)*jet1->pt, jet1->eta, jet1->phi, (1.-unc1)*jet1->mass);	

	  // calculate Up/Down variation for sub-leading jet
      //28feb    double unc2 = func(jet2->pt, jet2->eta, fJetUnc_AK4chs); 
//28feb	  VBF2_jes_up.SetPtEtaPhiM((1.+unc2)*jet2->pt, jet2->eta, jet2->phi, (1.+unc2)*jet2->mass);	
//28feb	  VBF2_jes_dn.SetPtEtaPhiM((1.-unc2)*jet2->pt, jet2->eta, jet2->phi, (1.-unc2)*jet2->mass);	

	 //if (TOT.M()<500) continue;
	  if ( VBFSel==1)
	  {
		if (TOT.Pt() < tempPtMax) continue;
		tempPtMax = TOT.Pt(); //take the jet pair with largest Pt
	  }
	  else if ( VBFSel==2)
	  {
		if (TOT.M() < tempPtMax) continue;
		tempPtMax = TOT.M(); //take the jet pair with largest mass
	  }
	  else if ( VBFSel==3)
	  {
		DRvbf = abs(VBF1.Eta()-VBF2.Eta());
		if (DRvbf < tempPtMax) continue;
		tempPtMax = DRvbf; //take the jet pair with largest dEta_jj
	  }
	  else
	  {
	  	cout<<"\n\nERROR:	Enter valid vbf selection criteria....\n\n"<<endl;
		exit(0);
	  }
          nVBF1 = indexGoodVBFJets.at(i); //save position of the 1st vbf jet
          nVBF2 = indexGoodVBFJets.at(ii); //save position of the 2nd vbf jet
//cout<<"     "<<nVBF1<<"     "<<nVBF2<<endl;
        }//loop q ends
      }//loop p ends
      if (nVBF1!=-1 && nVBF2 !=-1) OnlyTwoVBFTypeJets=1;
//      cutEff[7]++;
//if (OnlyTwoVBFTypeJets!=1) continue;
//cutEff[7]++;
//cout<<nVBF1<<"               "<<nVBF2<<endl;

//
      if (nVBF1!=-1 && nVBF2!=-1) //save infos for vbf jet pair
      {//loop x begins
        const baconhep::TJet *jet1 = (baconhep::TJet*)((*jetArr)[nVBF1]);
	const baconhep::TJet *jet2 = (baconhep::TJet*)((*jetArr)[nVBF2]);
        VBF1.SetPtEtaPhiM(jet1->pt,jet1->eta,jet1->phi,jet1->mass);
        VBF2.SetPtEtaPhiM(jet2->pt,jet2->eta,jet2->phi,jet2->mass);
        TOT = VBF1 + VBF2;
//cout<<VBF1<<"               "<<VBF2<<endl;
	WWTree->vbf_maxpt_j1_pt = jet1->pt;
	WWTree->vbf_maxpt_j1_eta = jet1->eta;
	WWTree->vbf_maxpt_j1_phi = jet1->phi;
	WWTree->vbf_maxpt_j1_e = VBF1.E();
	WWTree->vbf_maxpt_j1_mass = VBF1.M();
	WWTree->vbf_maxpt_j1_bDiscriminatorCSV = jet1->csv;
	WWTree->vbf_maxpt_j1_charge = jet1->q;

	WWTree->vbf_maxpt_j2_pt = jet2->pt;
	WWTree->vbf_maxpt_j2_eta = jet2->eta;
	WWTree->vbf_maxpt_j2_phi = jet2->phi;
	WWTree->vbf_maxpt_j2_e = VBF2.E();
	WWTree->vbf_maxpt_j2_mass = VBF2.M();
	WWTree->vbf_maxpt_j2_bDiscriminatorCSV = jet2->csv;
	WWTree->vbf_maxpt_j2_charge = jet2->q;
//cout<<"check9--all good"<<endl;

	// calculate Up/Down variation for leading jet
//28feb        double unc1 = func(jet1->pt, jet1->eta, fJetUnc_AK4chs); 
//28feb	VBF1_jes_up.SetPtEtaPhiM((1.+unc1)*jet1->pt, jet1->eta, jet1->phi, (1.+unc1)*jet1->mass);	
//28feb	VBF1_jes_dn.SetPtEtaPhiM((1.-unc1)*jet1->pt, jet1->eta, jet1->phi, (1.-unc1)*jet1->mass);	
	
	// calculate Up/Down variation for sub-leading jet
//28feb        double unc2 = func(jet2->pt, jet2->eta, fJetUnc_AK4chs); 
//28feb	VBF2_jes_up.SetPtEtaPhiM((1.+unc2)*jet2->pt, jet2->eta, jet2->phi, (1.+unc2)*jet2->mass);
//28feb	VBF2_jes_dn.SetPtEtaPhiM((1.-unc2)*jet2->pt, jet2->eta, jet2->phi, (1.-unc2)*jet2->mass);	
//28feb
/*	WWTree->vbf_maxpt_j1_pt_jes_up 	= VBF1_jes_up.Pt();
	WWTree->vbf_maxpt_j1_eta_jes_up = VBF1_jes_up.Eta();
	WWTree->vbf_maxpt_j1_phi_jes_up = VBF1_jes_up.Phi();
	WWTree->vbf_maxpt_j1_mass_jes_up= VBF1_jes_up.M();
	WWTree->vbf_maxpt_j1_pt_jes_dn 	= VBF1_jes_dn.Pt();
	WWTree->vbf_maxpt_j1_eta_jes_dn = VBF1_jes_dn.Eta();
	WWTree->vbf_maxpt_j1_phi_jes_dn = VBF1_jes_dn.Phi();
	WWTree->vbf_maxpt_j1_mass_jes_dn= VBF1_jes_dn.M();

	WWTree->vbf_maxpt_j2_pt_jes_up 	= VBF2_jes_up.Pt();
	WWTree->vbf_maxpt_j2_eta_jes_up = VBF2_jes_up.Eta();
	WWTree->vbf_maxpt_j2_phi_jes_up = VBF2_jes_up.Phi();
	WWTree->vbf_maxpt_j2_mass_jes_up= VBF2_jes_up.M();
	WWTree->vbf_maxpt_j2_pt_jes_dn 	= VBF2_jes_dn.Pt();
	WWTree->vbf_maxpt_j2_eta_jes_dn = VBF2_jes_dn.Eta();
	WWTree->vbf_maxpt_j2_phi_jes_dn = VBF2_jes_dn.Phi();
	WWTree->vbf_maxpt_j2_mass_jes_dn= VBF2_jes_dn.M();*/

        WWTree->vbf_maxpt_jj_pt = TOT.Pt();
//28feb        WWTree->vbf_maxpt_jj_pt_jes_up = (VBF1_jes_up + VBF2_jes_up).Pt();
//28feb        WWTree->vbf_maxpt_jj_pt_jes_dn = (VBF1_jes_dn + VBF2_jes_dn).Pt();

        WWTree->vbf_maxpt_jj_eta = TOT.Eta();
        WWTree->vbf_maxpt_jj_phi = TOT.Phi();

        WWTree->vbf_maxpt_jj_m = TOT.M();	
//28feb        WWTree->vbf_maxpt_jj_m_jes_up = (VBF1_jes_up + VBF2_jes_up).M();
//28feb        WWTree->vbf_maxpt_jj_m_jes_dn = (VBF1_jes_dn + VBF2_jes_dn).M();

	WWTree->vbf_maxpt_jj_Deta = abs(VBF1.Eta() - VBF2.Eta());
//28feb	WWTree->vbf_maxpt_jj_Deta_jes_up = abs(VBF1_jes_up.Eta() - VBF2_jes_up.Eta());
//28feb	WWTree->vbf_maxpt_jj_Deta_jes_dn = abs(VBF1_jes_dn.Eta() - VBF2_jes_dn.Eta());

//28feb	WWTree->AK4_DR_GENRECO_11 = abs(deltaR(WWTree->AK4_1_eta_gen, WWTree->AK4_1_phi_gen, WWTree->vbf_maxpt_j1_eta, WWTree->vbf_maxpt_j1_phi));
//28feb	WWTree->AK4_DR_GENRECO_12 = abs(deltaR(WWTree->AK4_1_eta_gen, WWTree->AK4_1_phi_gen, WWTree->vbf_maxpt_j2_eta, WWTree->vbf_maxpt_j2_phi));
//28feb	WWTree->AK4_DR_GENRECO_21 = abs(deltaR(WWTree->AK4_2_eta_gen, WWTree->AK4_2_phi_gen, WWTree->vbf_maxpt_j1_eta, WWTree->vbf_maxpt_j1_phi));
//28feb	WWTree->AK4_DR_GENRECO_22 = abs(deltaR(WWTree->AK4_2_eta_gen, WWTree->AK4_2_phi_gen, WWTree->vbf_maxpt_j2_eta, WWTree->vbf_maxpt_j2_phi));
      }//loop x ends
    }//loop y ends
    indexGoodVBFJets.clear();
//cout<<jentry<<"    mass of 2 jets=="<<WWTree->vbf_maxpt_jj_m<<endl;
    ///////////////////////////////////////////////////////////////////////////////////////
    if (OnlyTwoVBFTypeJets == 1) WWTree->isVBF=1;
    if (OnlyTwoVBFTypeJets == 0 ) continue;
cutEff[7]++;
 if (WWTree->vbf_maxpt_jj_m<500)  continue;
//if (WWTree->vbf_maxpt_jj_Deta  <2.5) continue;
 // if (WWTree->l_pt2<0) 
  cutEff[8]++;
//cout<<"check10-all good"<<endl;

   // else cutEff[18]++;
    /////////////////////////////////////////////////////////////////////////////////////////
    //
    //		Calculate b-tag weight
    //		
    /////////////////////////////////////////////////////////////////////////////////////////
//28feb
 /*   vector<float> btagSFv       = getJetSFs(goodJetsv, bTagReader, "central", "central");
    vector<float> btagSFUpHFv   = getJetSFs(goodJetsv, bTagReader, "up",      "central");
    vector<float> btagSFDownHFv = getJetSFs(goodJetsv, bTagReader, "down",    "central");
    vector<float> btagSFUpLFv   = getJetSFs(goodJetsv, bTagReader, "central", "up");
    vector<float> btagSFDownLFv = getJetSFs(goodJetsv, bTagReader, "central", "down"); 
    
    WWTree->btag0Wgt       = getBtagEventReweightEtaBin(0, goodJetsv, btagSFv);
    WWTree->btag1Wgt       = getBtagEventReweightEtaBin(1, goodJetsv, btagSFv);
    WWTree->btag2Wgt       = getBtagEventReweightEtaBin(2, goodJetsv, btagSFv);
    
    if (isnan(WWTree->btag0Wgt) == 1 || WWTree->btag0Wgt == 0) { 
    cout<<"********************* btag weight is nan"<<endl;
   
    }

    WWTree->btag0WgtUpHF   = getBtagEventReweightEtaBin(0, goodJetsv, btagSFUpHFv);
    WWTree->btag0WgtDownHF = getBtagEventReweightEtaBin(0, goodJetsv, btagSFDownHFv);
    WWTree->btag0WgtUpLF   = getBtagEventReweightEtaBin(0, goodJetsv, btagSFUpLFv);
    WWTree->btag0WgtDownLF = getBtagEventReweightEtaBin(0, goodJetsv, btagSFDownLFv);
    WWTree->btag1WgtUpHF   = getBtagEventReweightEtaBin(1, goodJetsv, btagSFUpHFv);
    WWTree->btag1WgtDownHF = getBtagEventReweightEtaBin(1, goodJetsv, btagSFDownHFv);
    WWTree->btag1WgtUpLF   = getBtagEventReweightEtaBin(1, goodJetsv, btagSFUpLFv);
    WWTree->btag1WgtDownLF = getBtagEventReweightEtaBin(1, goodJetsv, btagSFDownLFv);
   
    btagSFv.clear();	btagSFUpHFv.clear();	btagSFDownHFv.clear();
    btagSFUpLFv.clear();	btagSFDownLFv.clear();
    /////////////////////////////////////////////////////////////////////////////////////////	END btag weight calculation
  
        WWTree->totalEventWeight_2Lep = WWTree->genWeight*WWTree->pu_Weight*WWTree->top1_NNLO_Weight*WWTree->top2_NNLO_Weight*WWTree->trig_eff_Weight*WWTree->id_eff_Weight*WWTree->trig_eff_Weight2*WWTree->id_eff_Weight2;
*/
    
    WWTree->nEvents = TotalNumberOfEvents;
    WWTree->nNegEvents = nNegEvents;
    WWTree->nTotEvents = std::atof(TotalNumberOfEntries.c_str());
    WWTree->nTotNegEvents = std::atof(TotalNumberOfNegativeEntries.c_str());

//if(VBF1.
if (fabs(VBF1.Eta() - VBF2.Eta())  < 2.5) continue;
cutEff[9]++;

WWTree->ZeppenfeldW1 =abs(LEP1.Eta() - (VBF1.Eta() + VBF2.Eta())/2.0)/abs(VBF1.Eta() - VBF2.Eta());
WWTree->ZeppenfeldW2 =abs(LEP2.Eta() - (VBF1.Eta() + VBF2.Eta())/2.0)/abs(VBF1.Eta() - VBF2.Eta());



if (max(WWTree->ZeppenfeldW1,WWTree->ZeppenfeldW2) >0.75) continue;

cutEff[10]++;

WWTree->dilep_m_allcuts = (LEP1+LEP2).M();
WWTree->vbf_maxpt_jj_m_allcuts = TOT.M();
    outTree->Fill();

    goodJetsv.clear();
    }//loop on entries end
    delete infile;
    infile=0, eventTree=0;
    /////////////////FILL THE TREE
  }//loop on file ends
  //delete puWeight;	delete puWeight_up;	delete puWeight_down;
  delete MCpu;	delete MCpu_up;	delete MCpu_down;
  delete puWeightsDown;	delete puWeightsUp;	delete puWeights;
  
  pileupFileMC->Close();
  file->Close();
  std::cout << "---------end loop on events------------" << std::endl;
  std::cout << std::endl;
  std::cout << "GEN events = " << count_genEvents << std::endl;


  
  std::cout << "----------------------" << std::endl;
  std::cout << " SUMMARY" << std::endl;
  std::cout << "----------------------" << std::endl;
  std::cout << std::endl;
  std::cout<<"MC matching: "<<(float)ok/(float)total<<std::endl;
  std::cout<<"negative events: "<<nNegEvents<<std::endl;
  std::cout << std::endl;
  std::cout<<"(0) all events:        "<<cutEff[0]<<"\t:\t"<<((float)cutEff[0]*100.0)/(float)cutEff[0]<<std::endl
  	   <<"(1) Gen Events:        "<<cutEff[1]<<"\t:\t"<<((float)cutEff[1]*100.0)/(float)cutEff[0]<<std::endl
	   <<"(2) tight lepton:      "<<cutEff[2]<<"\t:\t"<<((float)cutEff[2]*100.0)/(float)cutEff[0]<<std::endl
	   <<"(3) DileptonMass:               "<<cutEff[3]<<"\t:\t"<<((float)cutEff[3]*100.0)/(float)cutEff[2]<<std::endl
	   <<"(4) (Dilepton-Z)Mass:  "<<cutEff[4]<<"\t:\t"<<((float)cutEff[4]*100.0)/(float)cutEff[3]<<std::endl
	   <<"(5) MET:               "<<cutEff[5]<<"\t:\t"<<((float)cutEff[5]*100.0)/(float)cutEff[4]<<std::endl
	   <<"(6) >=2 good VBF jets: "<<cutEff[6]<<"\t:\t"<<((float)cutEff[6]*100.0)/(float)cutEff[5]<<std::endl
	   <<"(7) =2 VBF jets      : "<<cutEff[7]<<"\t:\t"<<((float)cutEff[7]*100.0)/(float)cutEff[6]<<std::endl
	   <<"(8) Highest pt VBF jets with mjj>500:  "<<cutEff[8]<<"\t:\t"<<((float)cutEff[8]*100.)/(float)cutEff[7]<<std::endl
           <<"(9) Events with VBFjets delta eta>2.5"<<cutEff[9]<<"\t:\t"<<((float)cutEff[9]*100.)/(float)cutEff[8]<<std::endl
           <<"(9) ZeppenCut:                       "<<cutEff[10]<<"\t:\t"<<((float)cutEff[10]*100.)/(float)cutEff[9]<<std::endl;
std::cout << "Yield =  " << cutEff[10]*weight*35900<<endl;
  //--------close everything-------------
  delete info; delete gen;
  delete genPartArr; delete muonArr; delete electronArr; delete vertexArr;
  delete jetArr;



// delete vjetArrPuppi; delete vjetAddArrPuppi; 
delete lheWgtArr;
  outROOT->Write();
  outROOT->Close();
  int t1 = time(NULL);
  printf ("\n==> time to run this code = %0.3f min\n", (float)(t1 - t0)/60.0);
  return(0);
}
