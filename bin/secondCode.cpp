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
/*  TLorentzVector LEP1, LEP2;// SJ1_PuppiAK8, SJ2_PuppiAK8, SJ1, SJ2;
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
*/
  int ok=0, total=0;
  
  // Data structures to store info from TTrees
//  baconhep::TEventInfo *info  	= new baconhep::TEventInfo();
  baconhep::TGenEventInfo *gen	= new baconhep::TGenEventInfo();
  TClonesArray *genPartArr 	= new TClonesArray("baconhep::TGenParticle");
//  TClonesArray *muonArr    	= new TClonesArray("baconhep::TMuon");
//  TClonesArray *electronArr	= new TClonesArray("baconhep::TElectron");
//  TClonesArray *vertexArr	= new TClonesArray("baconhep::TVertex");
//  TClonesArray *jetArr		= new TClonesArray("baconhep::TJet");
//  TClonesArray *vjetArrPuppi	= new TClonesArray("baconhep::TJet");
//  TClonesArray *vjetAddArrPuppi	= new TClonesArray("baconhep::TAddJet");
//  TClonesArray *lheWgtArr	= new TClonesArray("baconhep::TLHEWeight");
  

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
/*TH1 *histo_aqgc[718];
  for(int j=0;j<718;j++)
    {
      stringstream ss;
      ss << j;
      string temp = ss.str();
      const char* name = temp.c_str();
      histo_aqgc[j] = new TH1D();
      histo_aqgc[j]->Sumw2();
}
*/

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
/*  TFile* IDIsoEle = TFile::Open("egammaEffi_EGM2D_TightCutBasedIDSF.root","READ");
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

*/
  //---------output tree----------------
  TFile* outROOT = TFile::Open((outputFile+(".root")).c_str(),"recreate");
  outROOT->cd();
  TTree* outTree = new TTree("otree", "otree");
  setOutputTree* WWTree = new setOutputTree(outTree);

  int nEvents=0;	
  Long64_t jentry2=0;
  int count_genEvents=0;

  int nInputFiles = sampleName.size();
cout <<"number of files =" <<nInputFiles<<endl;

  if (isLocal==1) nInputFiles = 21;
  cout<<"==> Total number of input files : "<<nInputFiles<<endl;

  TH1D *MCpu = new TH1D("MCpu","",75,0,75);
  TH1D *MCpu_up = new TH1D("MCpu_up","",75,0,75);
  TH1D *MCpu_down = new TH1D("MCpu_down","",75,0,75);
  
  Long64_t TotalNumberOfEvents = 0, nNegEvents = 0; 


  // Set up b-tag scale factor readers
 /* BTagCalibration calib("csvv2", "CSVv2_Moriond17_B_H.csv");
  BTagCalibrationReader bTagReader(BTagEntry::OP_LOOSE,  // working point: can be OP_LOOSE, OP_MEDIUM, OP_TIGHT 
                                   "central",             // label for the central value (see the scale factor file)
                                   {"up","down"});        // vector of labels for systematics
  bTagReader.load(calib, BTagEntry::FLAV_B, "comb");      // use the "comb" measurements for b-jets
  bTagReader.load(calib, BTagEntry::FLAV_C, "comb");      // use the "comb" measurements for c-jets
  bTagReader.load(calib, BTagEntry::FLAV_UDSG, "incl");   // use the "incl" measurements for light jets
  
*/
  vector<const baconhep::TJet*> goodJetsv;
  
  // Loop on input files
  for(int i=0;i<nInputFiles;i++)
  {//loop10 begins
     infile = TFile::Open(sampleName[i]);
     eventTree = (TTree*)infile->Get("Events");
     
     TotalNumberOfEvents+=eventTree->GetEntries();
     if(isMC)
     { 
//        eventTree->SetBranchAddress("Info", &info);    TBranch *infoBr = eventTree->GetBranch("Info");
  	TBranch *genBr=0;
     	eventTree->SetBranchAddress("GenEvtInfo", &gen); genBr = eventTree->GetBranch("GenEvtInfo");
	for (Long64_t jentry=0; jentry<eventTree->GetEntries();jentry++,jentry2++)
	{
	  //eventTree->GetEntry(jentry);
	    genBr->GetEntry(jentry);
//    	    infoBr->GetEntry(jentry);	    
//	    MCpu->Fill(info->nPUmean);
//	    MCpu_up->Fill(info->nPUmeanp);
//	    MCpu_down->Fill(info->nPUmeanm);
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


/*  JetCorrectorParameters paramAK4puppi("Summer16_23Sep2016V4_MC_JEC/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFPuppi.txt");
  JetCorrectorParameters paramAK4chs("Summer16_23Sep2016V4_MC_JEC/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt");
  JetCorrectorParameters paramAK8chs("Summer16_23Sep2016V4_MC_JEC/Summer16_23Sep2016V4_MC_Uncertainty_AK8PFchs.txt");
  JetCorrectorParameters paramAK8puppi("Summer16_23Sep2016V4_MC_JEC/Summer16_23Sep2016V4_MC_Uncertainty_AK8PFPuppi.txt");
*/ // JetCorrectionUncertainty *fJetUnc_AK4chs = new JetCorrectionUncertainty(paramAK4chs);
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



//  eventTree->SetBranchAddress("Info", &info);    TBranch *infoBr = eventTree->GetBranch("Info");
//  eventTree->SetBranchAddress("Muon", &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");
//  eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
//  eventTree->SetBranchAddress("PV",   &vertexArr); TBranch *vertexBr = eventTree->GetBranch("PV");
//  eventTree->SetBranchAddress("AK4CHS",   &jetArr); TBranch *jetBr = eventTree->GetBranch("AK4CHS");    
 // eventTree->SetBranchAddress("AK8Puppi",   &vjetArrPuppi); TBranch *vjetBrPuppi = eventTree->GetBranch("AK8Puppi");  
 // eventTree->SetBranchAddress("AddAK8Puppi",   &vjetAddArrPuppi); TBranch *vjetAddBrPuppi = eventTree->GetBranch("AddAK8Puppi");  
  TBranch *genBr=0, *genPartBr=0, *lhePartBr=0;
  if(isMC)
     {  //loop11 begins 
       eventTree->SetBranchAddress("GenEvtInfo", &gen); genBr = eventTree->GetBranch("GenEvtInfo");
       eventTree->SetBranchAddress("GenParticle",&genPartArr); genPartBr = eventTree->GetBranch("GenParticle");
//       if(eventTree->GetListOfBranches()->FindObject("LHEWeight"))
  //     {
   //    eventTree->SetBranchAddress("LHEWeight",&lheWgtArr); lhePartBr = eventTree->GetBranch("LHEWeight");	       }
     } //loop11 ends

  //for (Long64_t jentry=0; jentry<172; jentry++,jentry2++)
  for (Long64_t jentry=0; jentry<eventTree->GetEntries();jentry++,jentry2++)
  {
   // if (jentry2 != 87 && jentry2 != 113) continue;	// For debug
//if (jentry2>100) continue;
  //  infoBr->GetEntry(jentry);	    


    int GenPassCut = 0;

/*    tightMuon.clear();
    tightEle.clear();
    looseMuon.clear();
    looseEle.clear();
  */  
//cout<<"check1--all good"<<endl;
    if (jentry2%10000 == 0) std::cout << "\tread entry: " << jentry2 <<"/"<<TotalNumberOfEvents<<std:: endl;
    
    //*********************************
    WWTree->initializeVariables(); //initialize all variables

  //  WWTree->run   = info->runNum;
  //  WWTree->event = info->evtNum;
 //   WWTree->lumi  = info->lumiSec;


    /////////////////MC Info
    if (isMC==1)
    {//loop14 begins
   //   lheWgtArr->Clear();
     // if(lhePartBr)
//	{
//	  lhePartBr->GetEntry(jentry);
//	}
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
   /* 	for (int i = 0; i<lheWgtArr->GetEntries();i++)     // Note that i is starting from 446.
	{
		const baconhep::TLHEWeight *lhe = (baconhep::TLHEWeight*)((*lheWgtArr)[i]);
		//WWTree->LHEid[i] = lhe->id;
		WWTree->LHEWeight[i] = lhe->weight;
	}*/
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

    
  //  vertexArr->Clear();
 //   vertexBr->GetEntry(jentry);
 //   WWTree->nPV = vertexArr->GetEntries();
//  cout<<"check4--all good"<<endl;

    //PILE-UP WEIGHT
    if (isMC==1) { //loop17 begins
   /*    if(int(info->nPUmean)<75){
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
       } */
    }//loop17 ends

     // }

	
  
    //		MET JES Calculate
    //
    ////////////////////////////////////////////////////////////////
    ///// END: MET JES Calculate
   //<endl;

//
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
 // file->Close();
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
  	   <<"(1) Gen Events:        "<<cutEff[1]<<"\t:\t"<<((float)cutEff[1]*100.0)/(float)cutEff[0]<<std::endl;
	   
std::cout << "Yield =  " << cutEff[1]*weight*35900<<endl;
  //--------close everything-------------
//  delete info;
 delete gen;
  delete genPartArr;
// delete muonArr; delete electronArr; delete vertexArr;
//  delete jetArr;



// delete vjetArrPuppi; delete vjetAddArrPuppi; 
//delete lheWgtArr;
  outROOT->Write();
  outROOT->Close();
  int t1 = time(NULL);
  printf ("\n==> time to run this code = %0.3f min\n", (float)(t1 - t0)/60.0);
  return(0);
}
