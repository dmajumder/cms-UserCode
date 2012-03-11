#define REWEIGHTMC_CXX

#include "TRandom.h"

#include "reweightMC.h"
#include "HEADERDIR/format.h" 
#include "BTagWeight.h" 

#include <iostream>
#include <cmath>
#include <fstream>
#include <assert.h>

#include <TLorentzVector.h> 

using namespace std;

reweightMC :: reweightMC (const std::string& infiles,const bool& isData,const double& crosssection, const double& weight,const double& lumiInt,const double& events,const std::string&outfile) 
	:  bprimeAnalyzer(infiles,isData,crosssection,weight,lumiInt,events,outfile), 
  outfile_ ((char*)outfile.c_str(),"RECREATE")  
{
}

reweightMC :: ~reweightMC () { 
  outfile_.Write() ; 
  outfile_.Close() ; 
}

void reweightMC::initHists_ () {

  outfile_.cd();
  hmeeBjet_["hmeeBjet_JESHigh"]	= new TH1F ("hmeeBjet_JESHigh","",1000,0.,1000.) ; 
  hmeeBjet_["hmeeBjet_JESMean"]	= new TH1F ("hmeeBjet_JESMean","",1000,0.,1000.) ; 
  hmeeBjet_["hmeeBjet_JESLow"]	= new TH1F ("hmeeBjet_JESLow","",1000,0.,1000.) ; 
  hmeeBjet_["hmeeBjet_JERHigh"]	= new TH1F ("hmeeBjet_JERHigh","",1000,0.,1000.) ; 
  hmeeBjet_["hmeeBjet_JERMean"]	= new TH1F ("hmeeBjet_JERMean","",1000,0.,1000.) ; 
  hmeeBjet_["hmeeBjet_JERLow"]	= new TH1F ("hmeeBjet_JERLow","",1000,0.,1000.) ; 
	
}

void reweightMC::process_() {
  	
  initHists_() ; 	

  char inputFile[200] ; 
  ifstream infile ; 
  infile.open ((char*)inputfiles_.c_str(),ifstream::in) ; 
  assert(!infile.fail()) ;
  if (strstr((char*)inputfiles_.c_str(),"Reduced")!=0) {
    thisChain_ = new TChain("root") ;   
  } else {  
    thisChain_ = new TChain("zCandTree") ; 
  }
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

  return ; 	
}

void reweightMC :: evtLoop_ () {

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

  /**\ Looping over events */ 
  for (int iEntry=0;iEntry<thisChain_->GetEntries();++iEntry) { 

    LepInfoBranches* lepInfo = &LepInfo ; 
    JetInfoBranches* jetInfo = &JetInfo ; 

    double mee(0.) ; 
    double ptZ(0.) ; 

    thisChain_->GetEntry(iEntry) ; 

    /**\ Start reconstruction of B' candidates */ 
    for (int iEle=0; iEle < LepInfo.Size; ++iEle) { 
      //if (isGoodElectron_(lepInfo,iEntry,iEle) )
      { 
        /**\ Select second electron */ 
        for (int jEle=iEle+1; jEle < LepInfo.Size; ++jEle) { 
          if(   LepInfo.ChargeGsf[iEle]*LepInfo.ChargeGsf[jEle] < 0
            // && isGoodElectron_(lepInfo,iEntry,jEle) 
	    ) {  

	    TLorentzVector ele14v (LepInfo.Px[iEle],LepInfo.Py[iEle],LepInfo.Pz[iEle],LepInfo.Energy[iEle]) ; 
	    TLorentzVector ele24v (LepInfo.Px[jEle],LepInfo.Py[jEle],LepInfo.Pz[jEle],LepInfo.Energy[jEle]) ; 
            TLorentzVector zcand4v (ele14v+ele24v) ;
	    mee = zcand4v.Mag() ; 
            ptZ = zcand4v.Pt() ; 

            /**\ Sel 1: Select Z-candidates */
            if (   mee>60. && mee<120. 
	      /**\ Sel 2: Z-cands above 95 GeV/c */ 
	      && ptZ > 95.) { 

              TString jetModify = "JESHigh" ; 
              jetAndBprimeSel_(lepInfo,jetInfo,iEntry,iEle,jEle,jetModify) ; 

              jetModify = "JESMean" ; 
              jetAndBprimeSel_(lepInfo,jetInfo,iEntry,iEle,jEle,jetModify) ; 

              jetModify = "JESLow" ; 
              jetAndBprimeSel_(lepInfo,jetInfo,iEntry,iEle,jEle,jetModify) ; 

              jetModify = "JERHigh" ; 
              jetAndBprimeSel_(lepInfo,jetInfo,iEntry,iEle,jEle,jetModify) ; 

              jetModify = "JERMean" ; 
              jetAndBprimeSel_(lepInfo,jetInfo,iEntry,iEle,jEle,jetModify) ; 

              jetModify = "JERLow" ; 
              jetAndBprimeSel_(lepInfo,jetInfo,iEntry,iEle,jEle,jetModify) ; 

            } // Z-candidate accepted: MZ =[60,120] GeV/c2 and pTZ > 95 GeV/c 
          } // 2nd ele accepted 
        } // loop over all ele: 2 
      } // 1st ele accepted 
    } // loop over all ele: 1 

  } // event loop 

  return ; 

}

inline bool reweightMC::jetAndBprimeSel_(LepInfoBranches* LepInfo,JetInfoBranches* JetInfo,int& thisEntry,int& ele1,int& ele2, TString& jetModify) { 

  //thisChain_->GetEntry(thisEntry) ; 	

  //LepInfoBranches    LepInfo ; 
  //JetInfoBranches    JetInfo ; 

  //LepInfo->Register(thisChain_,"PFLepInfo") ; 
  //JetInfo->Register(thisChain_,"PFJetInfo") ; 

  int njets(0) ; 
  int nbjets(0) ; 
  
  /**\ Loop over jets */
  for (int iJet=0;iJet < JetInfo->Size;++iJet) {

    /**\ Remove jets overlapping with leptons (e/mu) */ 
    bool jetLepOverlap(false) ; 
    for (int iLep=0; iLep < LepInfo->Size; ++iLep) {
      if(bprimeAnalyzer::isGoodMuon_(LepInfo,thisEntry,iLep) || bprimeAnalyzer:: isGoodElectron_(LepInfo,thisEntry,iLep)) { 
        TLorentzVector lep4v (LepInfo->Px[iLep],LepInfo->Py[iLep],LepInfo->Pz[iLep],LepInfo->Energy[iLep]) ;
        TLorentzVector jet4v (JetInfo->Px[iJet],JetInfo->Py[iJet],JetInfo->Pz[iJet],JetInfo->Energy[iJet]) ; 
        if (lep4v.DeltaR(jet4v)<.5) { jetLepOverlap = true; break ; }
      }
    }
    if (jetLepOverlap) continue ;  // Reject jets overlapping with leptons 

    double jetPtNew(0.) ; 

    if (jetModify=="JESHigh") {
      jetPtNew = JetInfo->Pt[iJet]*1.05 ; 
    }
  
    if (jetModify=="JESMean") {
      jetPtNew = JetInfo->Pt[iJet]*1. ; 
    }
  
    if (jetModify=="JESLow") {
      jetPtNew = JetInfo->Pt[iJet]*.95 ; 
    }
  
    if (jetModify=="JERHigh") { 
      TRandom* gausJER = new TRandom() ; 	    
      jetPtNew = gausJER->Gaus(JetInfo->Pt[iJet],.2*JetInfo->Pt[iJet]) ; 
    }
  
    if (jetModify=="JERMean") {
      TRandom* gausJER = new TRandom() ; 	    
      jetPtNew = gausJER->Gaus(JetInfo->Pt[iJet],.1*JetInfo->Pt[iJet]) ; 
    }
  
    if (jetModify=="JERLow") {
      jetPtNew = JetInfo->Pt[iJet] ; 
    }
  
    if (jetPtNew < 30. || 
        TMath::Abs(JetInfo->Eta[iJet])>2.4 || 
        JetInfo->CHF[iJet]<=0.00 ||
        JetInfo->NHF[iJet]>=0.99 ||
        JetInfo->CEF[iJet]>=0.99 ||
        JetInfo->NEF[iJet]>=0.99 ||
        JetInfo->NCH[iJet]<=0.00  
    ) { continue; }  // apply jet ID 
    
    ++njets; 

    /**\ Pick b-jets here */
    if (JetInfo->TrackCountHiPurBJetTags[iJet]<=1.93) { continue ; } // apply b-tagging criteria 

    ++nbjets ; 

    /**\ Sel 3: Requiring at least one b-jet and one other jet */ 
    if (   njets >=2 && nbjets >=1 
        /**\ Sel 4: Select Bprime candidates */ 
	&& jetPtNew  > 65.) {  
    
        /**\ Reconstruct Bprime candidate */
        double meeBjetEnergy = LepInfo->Energy[ele1]+LepInfo->Energy[ele2]+JetInfo->Energy[iJet] ; 
        double meeBjetPx = LepInfo->Px[ele1]+LepInfo->Px[ele2]+JetInfo->Px[iJet] ; 
        double meeBjetPy = LepInfo->Py[ele1]+LepInfo->Py[ele2]+JetInfo->Py[iJet] ; 
        double meeBjetPz = LepInfo->Pz[ele1]+LepInfo->Pz[ele2]+JetInfo->Pz[iJet] ; 
    
        TLorentzVector bprime4v (meeBjetPx,meeBjetPy,meeBjetPz,meeBjetEnergy) ; 
	double bprimeMass = bprime4v.Mag() ; 
	TString histname = "hmeeBjet_"+jetModify ; 
	hmeeBjet_[histname]->Fill(bprimeMass) ; 
    
    } // Njets >= 2 && Nbjets >= 1 Bprime candidates accepted  
  } // Loop over all jets 

  return true ; 

}

