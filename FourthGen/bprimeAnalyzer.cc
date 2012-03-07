#define BPRIMEANALYZER_00_00_00_CXX
#include "bprimeAnalyzer.h"
#include "mathutils.h"

#include "HEADERDIR/format.h" 
#include "HEADERDIR/TriggerBooking.h" 
#include "BTagWeight.h" 

#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <assert.h>

#include <TLorentzVector.h> 

using namespace std; 

const int beginRunNo[11] = {
160431,   
161217,   
163270,   
165088,   
165970,   
165970,   
170826,   
173236,   
175832,   
178420,   
179959  
} ;

const int endRunno[11] = {
161176,  
163261,  
163869,  
165633,  
166967,  
166967,  
173198,  
173692,  
178380,  
179889,  
180252 
} ; 

const int hltBitIndex[11] = {
HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1,  
HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2,  
HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3,  
HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4,  
HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5,  
HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6,  
HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7,  
HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8,  
HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8,  
HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9,  
HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10 
} ; 

const char* hltPaths[11] = {
"HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1",  
"HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2",  
"HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3",  
"HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4",  
"HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5",  
"HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6",  
"HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7",  
"HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8",  
"HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8",  
"HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9",  
"HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10" 
} ; 

void generateWeightsPU(int PUmode, double PUWeight3D[50][50][50]) {

  TFile *f1 = new TFile("/afs/cern.ch/user/d/devdatta/Analysis/FourthGen/PUWeight3D.root"); 
  TH3D *WHist=(TH3D *)f1->Get("WHist");

  if(PUmode==1) {
    for(int ii=0;ii<50;ii++)
      for(int jj=0;jj<50;jj++)
        for(int kk=0;kk<50;kk++)
          PUWeight3D[ii][jj][kk]=WHist->GetBinContent(ii+1,jj+1,kk+1);
   }

  return ; 

}

double getWeightPU (int pv1, int pv2, int pv3, double PUWeight3D[50][50][50]) {

  int npm1 = min(pv1,49);
  int np0  = min(pv2,49);
  int npp1 = min(pv3,49);
  return (npm1>-1 && np0>-1 && npp1>-1) ? PUWeight3D[npm1][np0][npp1] : 0; 

}

double getWeightBtag (int& bq,int& cq,int& lq,std::vector<double>& v_bpt,std::vector<double>& v_cpt,std::vector<double>& v_lpt,std::vector<double>& v_leta) { 
  double weightBtag(0.) ;
  weightBtag =  (bq+cq+lq)*v_bpt.size()+v_cpt[0]+v_lpt[1]+v_leta.size() ; 
  return weightBtag ; 
}

Histograms::Histograms () {
  hPUWeight_      = new TH1F("hPUWeight","evt weight for PU",100,0.,2.)	; 
  hnPU_           = new TH1F("hnPU","no. of PU",50,-0.5,49.5) ; 
  hnGoodVtxs_     = new TH1F("hnGoodVtxs","no. of good 1ry vtxs",51,-0.5,50.5) ; 
  hlepId_         = new TH1F("hlepId","Lepton PDG ID",30,-15.,15.) ; 
  hnRecoEl_       = new TH1F("hnRecoEl","no. of reco el.",11,-.5,10.5) ; 
  hnPFEl_         = new TH1F("hnPFEl","no. of PF el.",11,-.5,10.5) ; 
  h1stElPt_       = new TH1F("h1stElPt","p_{T}^{1st e} (GeV/c)",1000,0.,1000.) ; 
  h2ndElPt_       = new TH1F("h2ndElPt","p_{T}^{2nd e} (GeV/c)",1000,0.,1000.) ; 
  h1stElEta_      = new TH1F("h1stElEta","#eta^{1st e} ",200,-4.,4.) ; 
  h2ndElEta_      = new TH1F("h2ndElEta","#eta^{2nd e} ",200,-4.,4.) ; 
  hnZcands_       = new TH1F("hnZcands","no. of Z candidates per event",20,-0.5,19.5) ; 
  hZPt_           = new TH1F("hZPt","p_{T}^{Z} (GeV/c)",1000,0.,1000.) ; 
  hZEta_          = new TH1F("hZEta","#eta^{Z}",200,-4.,4.) ; 
  hZPhi_          = new TH1F("hZPhi","#phi^{Z}",200,-4.,4.) ; 
  hmee_           = new TH1F("hmee","M_{e^{+}e^{-}} (GeV)",200,0.,200.) ; 
  hnJets_         = new TH1F("hnJets","no. of jets",21,-.5,20.5) ; 
  hallJetsPt_     = new TH1F("hallJetsPt","p_{T}^{all jets} (GeV/c)",1000,0.,1000.) ; 
  hallJetsEta_    = new TH1F("hallJetsEta","#eta^{all jets} ",200,-4.,4.) ; 
  hnBjets_        = new TH1F("hnBjets","no. of b-jets",11,-.5,10.5) ; 
  hallBJetsPt_    = new TH1F("hallBJetsPt","p_{T}^{Bjet} (GeV/c)",1000,0.,1000.) ; 
  hallBJetsEta_   = new TH1F("hallBJetsEta","#eta^{Bjets} ",200,-4.,4.) ; 
  hallBJetsPhi_   = new TH1F("hallBJetsPhi","#phi^{Bjets} ",200,-4.,4.) ; 
  hleadingBJetPt_ = new TH1F("hleadingBJetPt","p_{T}^{b-jet} (GeV/c)",1000,0.,1000.) ; 
  hleadingBJetEta_= new TH1F("hleadingBJetEta","#eta^{b-jet} ",200,-4.,4.) ; 
  hleadingBJetPhi_= new TH1F("hleadingBJetPhi","#phi^{b-jet} ",200,-4.,4.) ; 
  hnBprimes_      = new TH1F("hnBprimes","no. of B' candidates per event",20,-0.5,19.5) ; 
  hmeeBjet_       = new TH1F("hmeeBjet","M_{e^{+}e^{-}j} (GeV)",1000,0.,1000.) ; 
  heeBjetPt_      = new TH1F("heeBjetPt","p_{T}^{eeBjet} (GeV/c)",1000,0.,1000.) ; 
  heeBjetEta_     = new TH1F("heeBjetEta","#eta^{eeBjJet} ",200,-4.,4.) ; 
  heeBjetPhi_     = new TH1F("heeBjetPhi","#phi^{eeBjet} ",200,-4.,4.) ; 
}

void Histograms::fill () {

} ;

/**\ Default constructor */ 
bprimeAnalyzer :: bprimeAnalyzer (const std::string& infiles,const bool& isData,const double& crosssection, const double& weight,const double& lumiInt,const double& events,const std::string&outfile) :
  inputfiles_ (infiles),
  outfile_ ((char*)outfile.c_str(),"RECREATE"),
  isData_ (isData),
  thisChain_ (0),   	
  bprimeEvtsTree_(0), 
  evtProc_ (0), 
  evtAccept_ (0), 
  nzcands_ (0),
  nzPtGt95_(0),
  nBprimeCands_ (0), 
  el1Index_ (0),
  el2Index_ (0), 
  jetIndex_ (0),  
  bjetIndex_ (0)  
{

  dataset_.cs = crosssection;
  dataset_.wt = weight ;
  dataset_.lumi = lumiInt ;
  dataset_.evts = events ;

  setHLTBits_(); 

  initHists_(histBeforeSelection_,"Before selection");
  initHists_(histElSelection_,"El. selection"); 
  initHists_(histZCand_,"ZCand");  
  initHists_(histZPtGt95_,"ZPtGt95");  
  initHists_(histBprimeCand_,"BprimeCand");  

}

 /**\ Default destructor */ 
bprimeAnalyzer :: ~bprimeAnalyzer () {
	outfile_.Write() ; 
	outfile_.Close() ; 
        std::cout << "\n total events "  << dataset_.evts << std::endl ; 
        std::cout << "\n evt proc "    << evtProc_ << std::endl ;
        std::cout << "\n evt accepted "    << evtAccept_ << std::endl ;
        std::cout << "\n nZcands "    << nzcands_ << std::endl ;
        std::cout << "\n nZPtGt95 "    << nzPtGt95_ << std::endl ;
        std::cout << "\n nBprimeCands "    << nBprimeCands_ << std::endl ;

	std::cout << " Events after all cuts for " << dataset_.lumi << "/pb int. lumi = " << (dataset_.cs*dataset_.wt*dataset_.lumi*evtAccept_)/dataset_.evts << std::endl ; 
	std::cout << " Events with Z-candidates for " << dataset_.lumi << "/pb int. lumi = " << (dataset_.cs*dataset_.wt*dataset_.lumi*nzcands_)/dataset_.evts << std::endl ; 
	std::cout << " Events with Z with pTZ>95GeV for " << dataset_.lumi << "/pb int. lumi = " << (dataset_.cs*dataset_.wt*dataset_.lumi*nzPtGt95_)/dataset_.evts << std::endl ; 
	std::cout << " No. of B' candidates for " << dataset_.lumi << "/pb int. lumi = " << (dataset_.cs*dataset_.wt*dataset_.lumi*nBprimeCands_)/dataset_.evts << std::endl ; 
	std::cout << " Ending this instance of bprimeAnalyzer\n" ;
}

/**\ Process files */
void bprimeAnalyzer :: process_ () {

  char inputFile[200] ; 
  ifstream infile ; 
  infile.open ((char*)inputfiles_.c_str(),ifstream::in) ; 
  assert(!infile.fail()) ;
  if (strstr((char*)inputfiles_.c_str(),"Reduced")!=0) {
    thisChain_ = new TChain("root") ;   
  } else {  
    thisChain_ = new TChain("bprimeKit/root") ; 
  }
  //infile >> inputFile ; 
  //if (strncmp(inputFile,"#",1)!=0 || strncmp(inputFile,"%",1)!=0)  { 
  //  std::cout << " Rootfile read " << inputFile << std::endl; 
  //  thisChain_->Add(inputFile) ; 
  //  evtProc_ += thisChain_->GetEntries();
  //  std::cout << "\ninfile " << inputFile << std::endl ; 
  //  evtLoop_() ;
  //  thisChain_->Reset() ; 
  //}
  while (infile.good()) { 
    infile >> inputFile ; 
    if(strncmp(inputFile,"#",1)==0) continue; 
    if(strncmp(inputFile,"%",1)==0) break ; 
    thisChain_->Add(inputFile); 
    std::cout << "\ninfile " << inputFile << std::endl ; 
  }
  infile.close() ;   
  evtProc_ = thisChain_->GetEntries(); 
  evtLoop_();
  thisChain_->Reset() ;

  if ( !isData_ ) { 
    normalizeHists_(histBeforeSelection_) ; 
    normalizeHists_(histElSelection_) ; 
    normalizeHists_(histZCand_) ; 
    normalizeHists_(histZPtGt95_) ; 
    normalizeHists_(histBprimeCand_) ; 
  }	  

  return ; 

}

void bprimeAnalyzer::initHists_(Histograms*& h,std::string str) {

  outfile_.cd();
  outfile_.mkdir((char*)str.c_str());
  outfile_.cd((char*)str.c_str()) ; 
  h = new Histograms();
  outfile_.cd();
  
  return; 

}

void bprimeAnalyzer::fillHists_ (Histograms* H) { return; } 

void bprimeAnalyzer::normalizeHists_ (Histograms* H) {

  const double weight = (dataset_.cs*dataset_.wt*dataset_.lumi)/dataset_.evts ; 	
        
  H->hPUWeight_      ->Scale(weight) ; 
  H->hnPU_           ->Scale(weight) ; 
  H->hnGoodVtxs_     ->Scale(weight) ; 
  H->hlepId_         ->Scale(weight) ; 
  H->h1stElPt_       ->Scale(weight) ; 
  H->h2ndElPt_       ->Scale(weight) ; 
  H->h1stElEta_      ->Scale(weight) ; 
  H->h2ndElEta_      ->Scale(weight) ; 
  H->hnZcands_       ->Scale(weight) ; 
  H->hZPt_           ->Scale(weight) ; 
  H->hZEta_          ->Scale(weight) ; 
  H->hZPhi_          ->Scale(weight) ; 
  H->hmee_           ->Scale(weight) ; 
  H->hnJets_         ->Scale(weight) ; 
  H->hallJetsPt_     ->Scale(weight) ; 
  H->hallJetsEta_    ->Scale(weight) ; 
  H->hnBjets_        ->Scale(weight) ; 
  H->hallBJetsPt_    ->Scale(weight) ; 
  H->hallBJetsEta_   ->Scale(weight) ; 
  H->hallBJetsPhi_   ->Scale(weight) ; 
  H->hleadingBJetPt_ ->Scale(weight) ; 
  H->hleadingBJetEta_->Scale(weight) ; 
  H->hleadingBJetPhi_->Scale(weight) ; 
  H->hnBprimes_      ->Scale(weight) ; 
  H->hmeeBjet_       ->Scale(weight) ; 
  H->heeBjetPt_      ->Scale(weight) ; 
  H->heeBjetEta_     ->Scale(weight) ; 
  H->heeBjetPhi_     ->Scale(weight) ; 

  return ; 

}

/**\ Set HLT bits */
void bprimeAnalyzer::setHLTBits_ () {

  for (int ii=0;ii<11;++ii) {
    std::cout << " trig bit " << hltBitIndex[ii] << std::endl ; 
  }

  return ; 	
}

/**\ Event loop Apply Cuts fill Histograms */
void bprimeAnalyzer::evtLoop_ () { 

  GenInfoBranches    GenInfo ; 
  EvtInfoBranches    EvtInfo ; 
  VertexInfoBranches VtxInfo ; 
  LepInfoBranches    LepInfo ; 
  JetInfoBranches    JetInfo ; 

  GenInfo.Register(thisChain_) ; 
  EvtInfo.Register(thisChain_) ; 
  VtxInfo.Register(thisChain_) ; 
  LepInfo.Register(thisChain_,"PFLepInfo") ; 
  JetInfo.Register(thisChain_,"PFJetInfo") ; 

  double PUWeight3D[50][50][50];
  for(int ii=0;ii<50;ii++)
    for(int jj=0;jj<50;jj++)
      for(int kk=0;kk<50;kk++) 
        PUWeight3D[ii][jj][kk]=0.;      
  generateWeightsPU(1, PUWeight3D);

  /**\ Event selection flags */
  bool hasGoodEPpair(false) ; 
  bool hasZCand(false) ; 
  bool hasZCandPtGt95(false) ; 
  bool hasBjetJetPtGt30(false) ; 
  bool hasBJetPtGt65(false) ; 
  bool hasBprimeCand(false) ; 

  /**\ Event variables */ 
  double weightPU(0.) ; 
  double weightBtagging(0.) ; 
  double weightBtaggingErrorP(0.) ; 
  double weightBtaggingErrorM(0.) ; 
  int nZCands(0) ; 
  double mee(0.) ; 
  double ptZ(0.) ; 
  double etaZ(0.) ; 
  double phiZ(0.) ;  
  double njets(0) ; 
  int    nbjets(0) ; 
  int nBprimeCands(0.) ; 
  double meeBjet(0.) ;  
  double pt_eeBjet(0.) ; 
  double eta_eeBjet(0.) ; 
  double phi_eeBjet(0.) ; 
  int nGoodVtxs(0) ; 
  int nZPtGt95(0) ;  

  //TFile* bprimef = new TFile("/nasdata2/devdatta/bprimeAnalyzerSkims/Data/bprimeCandSkims.root","RECREATE"); 
  TFile* bprimef = new TFile("/nasdata2/devdatta/bprimeAnalyzerSkims/MC/Signal/bprimeCandSkims.root","RECREATE"); 
  //TFile* bprimef = new TFile("/nasdata2/devdatta/bprimeAnalyzerSkims/MC/Signal/bprimeCandSkims.root","RECREATE"); 

  TTree* zCandTree = thisChain_->CloneTree(0);
  TBranch *brzweightPU             = zCandTree->Branch("weightPU",             &weightPU,             "weightPU/D"); 
  zCandTree->SetName("zCandTree") ; 

  TTree* bprimeCandTree = thisChain_->CloneTree(0); 
  TBranch *brbprimeweightPU             = bprimeCandTree->Branch("weightPU",             &weightPU,             "weightPU/D"); 
  TBranch *brbprimeweightBtagging       = bprimeCandTree->Branch("weightBtagging",       &weightBtagging      , "weightBtagging/D"); 
  TBranch *brbprimeweightBtaggingErrorP = bprimeCandTree->Branch("weightBtaggingErrorP", &weightBtaggingErrorP, "weightBtaggingErrorP/D"); 
  TBranch *brbprimeweightBtaggingErrorM = bprimeCandTree->Branch("weightBtaggingErrorM", &weightBtaggingErrorM, "weightBtaggingErrorM/D"); 
  bprimeCandTree->SetName("bprimeCandTree") ; 

  /**\ Looping over events */ 
  for (int iEntry=0;iEntry<thisChain_->GetEntries();++iEntry) { 

    hasGoodEPpair = false ; 
    hasZCand = false ; 
    hasZCandPtGt95 = false ; 
    hasBjetJetPtGt30 = false ; 
    hasBJetPtGt65 = false ; 
    hasBprimeCand = false ; 

    weightPU = 0. ; 
    weightBtagging = 1. ; 
    weightBtaggingErrorP = 0. ; 
    weightBtaggingErrorM = 0. ; 
    nZCands = 0 ; 
    mee = 0. ;
    ptZ = 0. ;
    etaZ = 0. ;
    phiZ = 0. ; 
    njets = 0 ; 
    nbjets = 0 ; 
    nBprimeCands = 0. ;
    meeBjet = 0.; 
    pt_eeBjet = 0. ; 
    eta_eeBjet = 0. ; 
    phi_eeBjet = 0. ; 
    nGoodVtxs = 0 ; 
    nZPtGt95 = 0 ;  

    el1Index_.clear() ; 
    el2Index_.clear() ; 
    jetIndex_.clear() ; 
    bjetIndex_.clear() ; 
    std::vector<JetVariables> JetCollection; 

    /**\ Btag weight variables */
    int bquarks(0) ;
    int cquarks(0) ;
    int lquarks(0) ;
    std::vector<double>v_bpt; 
    std::vector<double>v_cpt; 
    std::vector<double>v_lpt; 
    std::vector<double>v_leta; 
    bool getBtagWt(true) ; 

    thisChain_->GetEntry(iEntry) ; 

    /**\ Get PU weight of event for MC */ 
    if (!isData_) {
      int nm1 = -1; int n0 = -1; int np1 = -1;
      for(int i=0;i<EvtInfo.nBX;i++) {
        if(EvtInfo.BXPU[i]==-1)  { nm1=EvtInfo.nPU[i]; } 
        if(EvtInfo.BXPU[i]==0)   { n0=EvtInfo.nPU[i] ; } 
        if(EvtInfo.BXPU[i]==1)   { np1=EvtInfo.nPU[i]; } 
      }
      weightPU = getWeightPU(nm1,n0,np1,PUWeight3D); 
    } else { weightPU = 1; } 

    /**\ Appy HLT requirement */
    bool isPassHLT(false) ; 
    int runno = EvtInfo.RunNo ;  
    if (isData_) {
      int hltBlock(-1) ; 
      for (int iHLT=0;iHLT<11;++iHLT) {
        if (runno>=beginRunNo[iHLT] && runno<=endRunno[iHLT]) { hltBlock=iHLT; break; } 
      }
      if (hltBlock>=0 && EvtInfo.TrgBook[hltBitIndex[hltBlock]]==1) { 
        isPassHLT = true ; 
      } else {isPassHLT = false ; }  
    } else { // isMC 
      for (int iHLT=0;iHLT<11;++iHLT) {
        if (EvtInfo.TrgBook[hltBitIndex[iHLT]]==1 ) { 
	  isPassHLT = true ; 
	  break ; 
	}
	else {isPassHLT = false ;  }  
      } 
    }

    if (isPassHLT==false) continue ;  

    /**\ Select good vertices */
    for (int iVtx=0; iVtx < VtxInfo.Size; ++iVtx) { 
      if (VtxInfo.Type[iVtx]==1 && VtxInfo.isFake[iVtx]==false && VtxInfo.Ndof[iVtx]>4 && VtxInfo.Rho[iVtx]<2. && VtxInfo.z[iVtx]<24.) { ++nGoodVtxs ; } 
    }

    if (isData_ && nGoodVtxs < 1)  { continue ; } 
            
    /**\ Fill histBeforeSelection_ */ 
    histBeforeSelection_->hPUWeight_->Fill(weightPU) ; 
    for (int iBX=0; iBX < EvtInfo.nBX; ++iBX) {
      histBeforeSelection_->hnPU_->Fill(EvtInfo.nPU[iBX],weightPU) ; 
    }
    histBeforeSelection_->hnGoodVtxs_->Fill(nGoodVtxs,weightPU) ; 

    /**\ Start reconstruction of B' candidates */ 
    for (int iEle=0; iEle < LepInfo.Size; ++iEle) {
      if (//isGoodElectron_(iEle)
              TMath::Abs(LepInfo.LeptonType[iEle])==11 
          &&  (TMath::Abs(LepInfo.simpleEleId80cIso[iEle]-7) < 0.1 || TMath::Abs(LepInfo.simpleEleId80cIso[iEle]-5) < 0.1)
          &&  TMath::Abs(LepInfo.ElTrackDxy_BS[iEle]) < 0.04
          &&  LepInfo.ChargeGsf[iEle] == LepInfo.ChargeCtf[iEle]
          &&  LepInfo.ChargeGsf[iEle] == LepInfo.ChargeScPix[iEle]
          &&  LepInfo.TrackIso[iEle]/LepInfo.Et[iEle] < 0.2
          &&  LepInfo.EcalIso[iEle]/LepInfo.Et[iEle] < 0.2
          &&  LepInfo.HcalIso[iEle]/LepInfo.Et[iEle] < 0.2 
          &&  (LepInfo.ChargedHadronIso[iEle]+LepInfo.NeutralHadronIso[iEle]+LepInfo.PhotonIso[iEle])/TMath::Abs(LepInfo.Pt[iEle]) < 0.15 
          &&  LepInfo.Pt[iEle] > 25. 
          &&  ( TMath::Abs(LepInfo.Eta[iEle])<1.4442 || (TMath::Abs(LepInfo.Eta[iEle])>1.566 && TMath::Abs(LepInfo.Eta[iEle])<2.4) ) 
	 ) {
        el1Index_.push_back(iEle) ; 
        /**\ Select second electron */ 
        for (int jEle=iEle+1; jEle < LepInfo.Size; ++jEle) { 
          if(    LepInfo.ChargeGsf[iEle]*LepInfo.ChargeGsf[jEle] < 0
              //&& isGoodElectron_(jEle) 
	      && TMath::Abs(LepInfo.LeptonType[jEle])==11 
              && (TMath::Abs(LepInfo.simpleEleId80cIso[jEle]-7) < 0.1 || TMath::Abs(LepInfo.simpleEleId80cIso[jEle]-5) < 0.1)  
              && TMath::Abs(LepInfo.ElTrackDxy_BS[jEle]) < 0.04 
              && LepInfo.ChargeGsf[jEle] == LepInfo.ChargeCtf[jEle] 
              && LepInfo.ChargeGsf[jEle] == LepInfo.ChargeScPix[jEle] 
              && LepInfo.TrackIso[jEle]/LepInfo.Et[jEle] < 0.2 
	      && LepInfo.EcalIso[jEle]/LepInfo.Et[jEle] < 0.2 
              && LepInfo.HcalIso[jEle]/LepInfo.Et[jEle] < 0.2 
              && (LepInfo.ChargedHadronIso[jEle]+LepInfo.NeutralHadronIso[jEle]+LepInfo.PhotonIso[jEle])/TMath::Abs(LepInfo.Pt[jEle]) < 0.15 
              && LepInfo.Pt[jEle] > 25.
              && ( TMath::Abs(LepInfo.Eta[jEle])<1.4442 || (TMath::Abs(LepInfo.Eta[jEle])>1.566 && TMath::Abs(LepInfo.Eta[jEle])<2.4) ) 
	    ) {  
            el2Index_.push_back(jEle) ; 
            hasGoodEPpair = true ; 		  
	    TLorentzVector ele14v (LepInfo.Px[iEle],LepInfo.Py[iEle],LepInfo.Pz[iEle],LepInfo.Energy[iEle]) ; 
	    TLorentzVector ele24v (LepInfo.Px[jEle],LepInfo.Py[jEle],LepInfo.Pz[jEle],LepInfo.Energy[jEle]) ; 
            TLorentzVector zcand4v (ele14v+ele24v) ;
	    mee = zcand4v.Mag() ; 
            ptZ = zcand4v.Pt() ; 
	    etaZ =  zcand4v.Eta() ; 
	    phiZ = zcand4v.Phi() ; 

	    /**\ Fill histElSelection_ */ 
            for (int iBX=0; iBX < EvtInfo.nBX; ++iBX) {
              histElSelection_->hnPU_->Fill(EvtInfo.nPU[iBX],weightPU) ; 
            }
            histElSelection_->hnGoodVtxs_->Fill(nGoodVtxs,weightPU) ; 
            histElSelection_->h1stElPt_->Fill(LepInfo.Pt[iEle]) ; 
            histElSelection_->h2ndElPt_->Fill(LepInfo.Pt[jEle]) ; 
            histElSelection_->h1stElEta_->Fill(LepInfo.Eta[iEle]) ; 
            histElSelection_->h2ndElEta_->Fill(LepInfo.Eta[jEle]) ; 
            histElSelection_->hmee_->Fill(mee,weightPU); 
            histElSelection_->hZPt_->Fill(ptZ,weightPU); 
            histElSelection_->hZEta_->Fill(etaZ,weightPU); 
            histElSelection_->hZPhi_->Fill(phiZ,weightPU); 

            /**\ Sel 1: Select Z-candidates */
            if (mee>60. && mee<120.) { 
	      hasZCand = true ; 
	      ++nZCands ; 

	      /**\ Fill histZCand_ */ 
              histZCand_->h1stElPt_->Fill(LepInfo.Pt[iEle]) ; 
              histZCand_->h2ndElPt_->Fill(LepInfo.Pt[jEle]) ; 
              histZCand_->h1stElEta_->Fill(LepInfo.Eta[iEle]) ; 
              histZCand_->h2ndElEta_->Fill(LepInfo.Eta[jEle]) ; 
              histZCand_->hmee_->Fill(mee,weightPU); 
              histZCand_->hZPt_->Fill(ptZ,weightPU); 
              histZCand_->hZEta_->Fill(etaZ,weightPU); 
              histZCand_->hZPhi_->Fill(phiZ,weightPU); 

	      /**\ Sel 2: Z-cands above 95 GeV/c */ 
	      if (ptZ > 95.) { 
                hasZCandPtGt95 = true ; 		        
		++nZPtGt95 ; 

	        /**\ Loop over jets */
                for (int iJet=0;iJet < JetInfo.Size;++iJet) {

                  /**\ Remove jets overlapping with leptons (e/mu) */ 
	          bool jetLepOverlap(false) ; 
                  for (int iLep=0; iLep < LepInfo.Size; ++iLep) {

                  if(   LepInfo.LeptonType[iLep]==13
                     && LepInfo.Pt[iLep]>20.
                     && TMath::Abs(LepInfo.Eta[iLep])<2.1
                     && LepInfo.MuIDGlobalMuonPromptTight[iLep]==1
                     && LepInfo.MuInnerTrackNHits[iLep]>10
                     && TMath::Abs(LepInfo.MuInnerTrackDxy_BS[iLep])<0.2
                     && LepInfo.MuNPixelLayers[iLep]>=1
                     && LepInfo.MuNChambersMatchesSegment[iLep]>=2
                     && LepInfo.MuNMuonhits[iLep]>=1
                     && LepInfo.MuType[iLep]%8>=6
                     && (( LepInfo.TrackIso[iLep]+LepInfo.EcalIso[iLep]+LepInfo.HcalIso[iLep])/LepInfo.Pt[iLep]<0.15)
                    ) {
                      TLorentzVector lep4v (LepInfo.Px[iLep],LepInfo.Py[iLep],LepInfo.Pz[iLep],LepInfo.Energy[iLep]) ;
                      TLorentzVector jet4v (JetInfo.Px[iJet],JetInfo.Py[iJet],JetInfo.Pz[iJet],JetInfo.Energy[iJet]) ; 
                      if (lep4v.DeltaR(jet4v)<.5) jetLepOverlap = true; 
		  }

                  if (    TMath::Abs(LepInfo.LeptonType[iLep])==11 
                      &&  (TMath::Abs(LepInfo.simpleEleId80cIso[iLep]-7) < 0.1 || TMath::Abs(LepInfo.simpleEleId80cIso[iLep]-5) < 0.1)
                      &&  TMath::Abs(LepInfo.ElTrackDxy_BS[iLep]) < 0.04
                      &&  LepInfo.ChargeGsf[iLep] == LepInfo.ChargeCtf[iLep]
                      &&  LepInfo.ChargeGsf[iLep] == LepInfo.ChargeScPix[iLep]
                      &&  LepInfo.TrackIso[iLep]/LepInfo.Et[iLep] < 0.2
                      &&  LepInfo.EcalIso[iLep]/LepInfo.Et[iLep] < 0.2
                      &&  LepInfo.HcalIso[iLep]/LepInfo.Et[iLep] < 0.2 
                      &&  (LepInfo.ChargedHadronIso[iLep]+LepInfo.NeutralHadronIso[iLep]+LepInfo.PhotonIso[iLep])/TMath::Abs(LepInfo.Pt[iLep]) < 0.15 
                      &&  LepInfo.Pt[iLep] > 20. 
                      &&  ( TMath::Abs(LepInfo.Eta[iLep])<1.4442 || (TMath::Abs(LepInfo.Eta[iLep])>1.566 && TMath::Abs(LepInfo.Eta[iLep])<2.5) ) 
                     ) { 
                      TLorentzVector lep4v (LepInfo.Px[iLep],LepInfo.Py[iLep],LepInfo.Pz[iLep],LepInfo.Energy[iLep]) ;
                      TLorentzVector jet4v (JetInfo.Px[iJet],JetInfo.Py[iJet],JetInfo.Pz[iJet],JetInfo.Energy[iJet]) ; 
                      if (lep4v.DeltaR(jet4v)<.5) jetLepOverlap = true; 
		    }

	          }
	          if (jetLepOverlap) continue ;  // Reject jets overlapping with leptons 
            
                  if (JetInfo.Pt[iJet]<30. ||
                      TMath::Abs(JetInfo.Eta[iJet])>2.4 || 
                      JetInfo.CHF[iJet]<=0.00 ||
                      JetInfo.NHF[iJet]>=0.99 ||
                      JetInfo.CEF[iJet]>=0.99 ||
                      JetInfo.NEF[iJet]>=0.99 ||
                      JetInfo.NCH[iJet]<=0.00  
                  ) {
	            continue; 
	          }  // apply jet ID 
                  
                  ++njets ; 			
	          jetIndex_.push_back(iJet) ; 

	          /**\ Fill histZPtGt95_ */ 
	          histZPtGt95_->hallJetsPt_->Fill(JetInfo.Pt[iJet]) ; 
	          histZPtGt95_->hallJetsEta_->Fill(JetInfo.Eta[iJet]) ; 

	          /**\ Pick b-jets here */
	          if (JetInfo.TrackCountHiPurBJetTags[iJet]<=1.93) { continue ; } // apply b-tagging criteria 

                  ++nbjets ; 
	          bjetIndex_.push_back(iJet) ; 		    

	          /**\ Fill histZPtGt95_ */ 
                  histZPtGt95_->hallBJetsPt_->Fill(JetInfo.Pt[iJet],weightPU); 
	          histZPtGt95_->hallBJetsEta_->Fill(JetInfo.Eta[iJet],weightPU);
	          histZPtGt95_->hallBJetsPhi_->Fill(JetInfo.Phi[iJet],weightPU);
            
                } // Loop over all jets 

	        /**\ Fill histZPtGt95_ */ 
                for (int iBX=0; iBX < EvtInfo.nBX; ++iBX) {
                  histZPtGt95_->hnPU_->Fill(EvtInfo.nPU[iBX],weightPU) ; 
                }
                histZPtGt95_->hnGoodVtxs_->Fill(nGoodVtxs,weightPU) ; 
                histZPtGt95_->h1stElPt_->Fill(LepInfo.Pt[iEle]) ; 
                histZPtGt95_->h2ndElPt_->Fill(LepInfo.Pt[jEle]) ; 
                histZPtGt95_->hmee_->Fill(mee,weightPU); 
                histZPtGt95_->hZPt_->Fill(ptZ,weightPU); 
                histZPtGt95_->hnBjets_->Fill(nbjets,weightPU) ; 
                if (nbjets>0) histZPtGt95_->hleadingBJetPt_ ->Fill(JetInfo.Pt [bjetIndex_[0]],weightPU) ; 
                if (nbjets>0) histZPtGt95_->hleadingBJetEta_->Fill(JetInfo.Eta[bjetIndex_[0]],weightPU) ; 
                if (nbjets>0) histZPtGt95_->hleadingBJetPhi_->Fill(JetInfo.Phi[bjetIndex_[0]],weightPU) ; 

	        /**\ Sel 3: Requiring at least one b-jet and one other jet */ 
	        if (njets >=2 && nbjets >=1) { 

                  hasBjetJetPtGt30 = true ; 			

		  /**\ Getting Btag weight */
		  if (getBtagWt==true) {
                    getBtagWt = false ;			  
	            for (unsigned iJet=0;iJet< jetIndex_.size();++iJet) { 
                      bool btag(false) ; 
                      int flavour = TMath::Abs(JetInfo.GenPdgID[jetIndex_[iJet]]) ; 
                      double ptj = JetInfo.Pt[jetIndex_[iJet]] ;
		      double etaj = JetInfo.Eta[jetIndex_[iJet]] ;
	              for (unsigned iBjet=0;iBjet< bjetIndex_.size();++iBjet) { 
                        if ( jetIndex_[iJet]==bjetIndex_[iBjet] ) { 
                          btag = true ; 
		          break ; 
		        }
		      }
		      JetVariables Jet(btag,flavour,ptj,etaj) ; 
		      JetCollection.push_back(Jet) ; 
		    }
                    BTagWeight bTagWeight(JetCollection,isData_,1.93,0,1,nbjets) ; 
                    int ntags(1) ;  
                    weightBtagging = bTagWeight.weight(ntags) ; 
		  }

	          for (unsigned int iBjet=0;iBjet< bjetIndex_.size();++iBjet) { 

		    /**\ Compute b-tagging SF weight */ 
                    if (TMath::Abs(JetInfo.GenPdgID[bjetIndex_[iBjet]])==5) {
                      ++bquarks ; 
		      v_bpt.push_back(JetInfo.Pt[bjetIndex_[iBjet]]) ; 
		    }
                    if (TMath::Abs(JetInfo.GenPdgID[bjetIndex_[iBjet]])==4) {
		      ++cquarks ; 
		      v_cpt.push_back(JetInfo.Pt[bjetIndex_[iBjet]]) ; 
                    }
                    if (TMath::Abs(JetInfo.GenPdgID[bjetIndex_[iBjet]]) <4) { 
		      ++lquarks ; 
		      v_lpt.push_back(JetInfo.Pt[bjetIndex_[iBjet]]) ; 
		      v_leta.push_back(JetInfo.Eta[bjetIndex_[iBjet]]) ; 
		    }

                    /**\ Sel 4: Select Bprime candidates */ 
                    if (JetInfo.Pt[bjetIndex_[iBjet]] > 65.) { 
                      hasBJetPtGt65 = true ; 			    

	              /**\ Reconstruct Bprime candidate */
                      double meeBjetEnergy = LepInfo.Energy[iEle]+LepInfo.Energy[jEle]+JetInfo.Energy[bjetIndex_[iBjet]] ; 
                      double meeBjetPx = LepInfo.Px[iEle]+LepInfo.Px[jEle]+JetInfo.Px[bjetIndex_[iBjet]] ; 
                      double meeBjetPy = LepInfo.Py[iEle]+LepInfo.Py[jEle]+JetInfo.Py[bjetIndex_[iBjet]] ; 
                      double meeBjetPz = LepInfo.Pz[iEle]+LepInfo.Pz[jEle]+JetInfo.Pz[bjetIndex_[iBjet]] ; 

	              TLorentzVector bprime4v (meeBjetPx,meeBjetPy,meeBjetPz,meeBjetEnergy) ; 

                      meeBjet = bprime4v.Mag() ; 
                      pt_eeBjet = bprime4v.Pt() ; 
                      eta_eeBjet = bprime4v.Eta() ; 
                      phi_eeBjet = bprime4v.Phi() ; 

                      hasBprimeCand = true ; 			
                      ++nBprimeCands ; 	
                      ++nBprimeCands_ ; 

	              /**\ Fill histBprimeCand_ */ 
                      histZPtGt95_->hallBJetsPt_->Fill(JetInfo.Pt[bjetIndex_[iBjet]],weightPU); 
	              histZPtGt95_->hallBJetsEta_->Fill(JetInfo.Eta[bjetIndex_[iBjet]],weightPU); 
	              histZPtGt95_->hallBJetsPhi_->Fill(JetInfo.Phi[bjetIndex_[iBjet]],weightPU); 
                      histBprimeCand_->hmeeBjet_->Fill(meeBjet,weightPU) ;  
                      histBprimeCand_->heeBjetPt_->Fill(pt_eeBjet,weightPU) ;  
                      histBprimeCand_->heeBjetEta_->Fill(eta_eeBjet,weightPU) ;  
                      histBprimeCand_->heeBjetPhi_->Fill(phi_eeBjet,weightPU) ;  

	            } // Bprime candidates accepted 
	          } // loop over all b-jets 
	          if (hasBprimeCand) {  
	              /**\ Fill histBprimeCand_ */ 
                      histBprimeCand_->h1stElPt_->Fill(LepInfo.Pt[iEle]) ; 
                      histBprimeCand_->h2ndElPt_->Fill(LepInfo.Pt[jEle]) ; 
                      histBprimeCand_->hmee_->Fill(mee,weightPU); 
                      histBprimeCand_->hZPt_->Fill(ptZ,weightPU); 
                      histBprimeCand_->hleadingBJetPt_ ->Fill(JetInfo.Pt [bjetIndex_[0]],weightPU) ; 
                      histBprimeCand_->hleadingBJetEta_->Fill(JetInfo.Eta[bjetIndex_[0]],weightPU) ; 
                      histBprimeCand_->hleadingBJetPhi_->Fill(JetInfo.Phi[bjetIndex_[0]],weightPU) ; 
	          }
		} // Njets >= 2 && Nbjets >= 1 
	      } // pTZ > 95 GeV/c 
            } // Z-candidate accepted: MZ =[60,120] GeV/c2 
          } // 2nd ele accepted 
        } // loop over all ele: 2 
      } // 1st ele accepted 
    } // loop over all ele: 1 

    if (hasZCand) { 
      brzweightPU->Fill(); 	    
      zCandTree->Fill(); 
      ++nzcands_ ; 
      for (int iBX=0; iBX < EvtInfo.nBX; ++iBX) {
        histZCand_->hnPU_->Fill(EvtInfo.nPU[iBX],weightPU) ; 
      }
      histZCand_->hPUWeight_->Fill(weightPU) ; 
      histZCand_->hnGoodVtxs_->Fill(nGoodVtxs,weightPU) ;  
      histZCand_->hnZcands_->Fill(nZCands,weightPU) ;  
      histZCand_->hnJets_->Fill(njets,weightPU) ;  
      if (hasZCandPtGt95) ++nzPtGt95_ ; 
    } 
    if (hasBprimeCand) {
      brbprimeweightPU             ->Fill() ;  
      brbprimeweightBtagging       ->Fill() ; 
      brbprimeweightBtaggingErrorP ->Fill() ; 
      brbprimeweightBtaggingErrorM ->Fill() ; 
      bprimeCandTree->Fill(); 
      ++evtAccept_; 
      for (int iBX=0; iBX < EvtInfo.nBX; ++iBX) {
        histBprimeCand_->hnPU_->Fill(EvtInfo.nPU[iBX],weightPU) ; 
      }
      histBprimeCand_->hPUWeight_->Fill(weightPU) ; 
      histBprimeCand_->hnGoodVtxs_->Fill(nGoodVtxs,weightPU) ;  
      histBprimeCand_->hnZcands_->Fill(nZCands,weightPU) ;  
      histBprimeCand_->hnJets_->Fill(njets,weightPU) ;  
      histBprimeCand_->hnBjets_->Fill(nbjets,weightPU) ;  
      histBprimeCand_->hnBprimes_->Fill(nBprimeCands,weightPU) ; 
    }

  } // event loop 

  zCandTree->Write("",TObject::kOverwrite);
  delete zCandTree; 
  bprimeCandTree->Write("",TObject::kOverwrite); 
  delete bprimeCandTree; 
  delete bprimef; 

}

bool bprimeAnalyzer:: isGoodElectron_ (int& iEle) { 

  bool isGoodElectron(false) ; 

  LepInfoBranches    LepInfo ;
  LepInfo.Register(thisChain_,"PFLepInfo") ; 

  if (    TMath::Abs(LepInfo.LeptonType[iEle])==11 
      &&  (TMath::Abs(LepInfo.simpleEleId80cIso[iEle]-7) < 0.1 || TMath::Abs(LepInfo.simpleEleId80cIso[iEle]-5) < 0.1)
      &&  TMath::Abs(LepInfo.ElTrackDxy_BS[iEle]) < 0.04
      &&  LepInfo.ChargeGsf[iEle] == LepInfo.ChargeCtf[iEle]
      &&  LepInfo.ChargeGsf[iEle] == LepInfo.ChargeScPix[iEle]
      &&  LepInfo.TrackIso[iEle]/LepInfo.Et[iEle] < 0.2
      &&  LepInfo.EcalIso[iEle]/LepInfo.Et[iEle] < 0.2
      &&  LepInfo.HcalIso[iEle]/LepInfo.Et[iEle] < 0.2 
      &&  (LepInfo.ChargedHadronIso[iEle]+LepInfo.NeutralHadronIso[iEle]+LepInfo.PhotonIso[iEle])/TMath::Abs(LepInfo.Pt[iEle]) < 0.15 
      &&  LepInfo.Pt[iEle] > 20. 
      &&  ( TMath::Abs(LepInfo.Eta[iEle])<1.4442 || (TMath::Abs(LepInfo.Eta[iEle])>1.566 && TMath::Abs(LepInfo.Eta[iEle])<2.5) ) 
     ) { 
   //if(   LepInfo.LeptonType[iEle]==11
   //   && LepInfo.Pt[iEle]>20.
   //   && ( TMath::Abs(LepInfo.Eta[iEle])<1.4442 || (TMath::Abs(LepInfo.Eta[iEle])>1.566 && TMath::Abs(LepInfo.Eta[iEle])<2.4) ) 
   //   && LepInfo.simpleEleId80relIso[iEle] == 7
   //   && TMath::Abs(LepInfo.ElTrackDxy_BS[iEle])<0.04
   //  ) {
     isGoodElectron = true ;	     
   } else isGoodElectron = false ; 

  return isGoodElectron ; 	

}

bool bprimeAnalyzer:: isGoodMuon_ (int& iMu) { 

  bool isGoodMuon(false) ; 

  LepInfoBranches    LepInfo ;
  LepInfo.Register(thisChain_,"PFLepInfo") ; 

   if(   LepInfo.LeptonType[iMu]==13
      && LepInfo.Pt[iMu]>20.
      && TMath::Abs(LepInfo.Eta[iMu])<2.1
      && LepInfo.MuIDGlobalMuonPromptTight[iMu]==1
      && LepInfo.MuInnerTrackNHits[iMu]>10
      && TMath::Abs(LepInfo.MuInnerTrackDxy_BS[iMu])<0.2
      && LepInfo.MuNPixelLayers[iMu]>=1
      && LepInfo.MuNChambersMatchesSegment[iMu]>=2
      && LepInfo.MuNMuonhits[iMu]>=1
      && LepInfo.MuType[iMu]%8>=6
      && (( LepInfo.TrackIso[iMu]+LepInfo.EcalIso[iMu]+LepInfo.HcalIso[iMu])/LepInfo.Pt[iMu]<0.15)
     ) {
    isGoodMuon = true ; 	   
   } else isGoodMuon = false ; 

  return isGoodMuon ; 	

}

int main (int argc, char **argv) { 

  time_t start, stop ; 
  double time_elapsed ; 
  time(&start) ; 

  std::string ip ; 
  std::string op ; 
  double cs(0.) ;
  double wt(0.) ; 
  bool isData(0) ;
  double lumi(0) ; 
  double events(0);

  std::cout << " enter input\n" ;
  std::cin >> ip ;
  std::cout << " enter output\n" ; 
  std::cin >> op ; 
  std::cout << " enter isData (=1 for data =0 for MC)\n" ;
  std::cin >> isData ;
  if (!isData) {
    std::cout << " enter sample cross-section\n" ;
    std::cin >> cs ; 
    std::cout << " enter sample weight\n" ;
    std::cin >> wt ; 
    std::cout << " enter integrated luminosity \n" ;
    std::cin >> lumi ; 
    std::cout << " enter total number of events\n" ;
    std::cin >> events;
  }

  bprimeAnalyzer* bprimeanalyzer = new bprimeAnalyzer(ip,isData,cs,wt,lumi,events,op) ; 
  bprimeanalyzer->process_() ; 
  delete bprimeanalyzer; 
  
  time(&stop) ; 
  time_elapsed = difftime(stop, start) ; 
  std::cout << "\n Time taken for program " << time_elapsed << " seconds \n" ;

  return 0 ; 

}

