#ifndef selectLeptons_cxx 
#define selectLeptons_cxx
#include "bprimeAnalyzer.h"
#include "mathutils.h"

#include "HEADERDIR/format.h"
#include "HEADERDIR/HitFitInfoBranches.h"

#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <assert.h>

#include <TLorentzVector.h> 

class selectLeptons : public bprimeAnalyzer {

	public:

	selectLeptons () ; 
	~selectLeptons () ; 

	private:

	TH1F* hnpfEle_ = new TH1F("npfEle","No. of PF electrons",21,0-0.5,20.5); 
	TH1F* hngsfEle_ = new TH1F("ngsfEle","No. of GSF electrons",21,0-0.5,20.5); 
	void evtLoop_ () ; 

};

selectleptons::selectleptons () {

  bprimeAnalyzer* bprimeanalyzer = new bprimeAnalyzer(ip,isData,cs,wt,lumi,events,op) ; 
  bprimeanalyzer->process_() ; 

}

selectleptons::~selectleptons () { 
}

/**\ Event loop Apply Cuts fill Histograms */
void selectLeptons:: evtLoop_ () { 

  GenInfoBranches GenInfo; 
  EvtInfoBranches    EvtInfo ; 
  VertexInfoBranches VtxInfo ; 
  LepInfoBranches    LepInfo ;
  LepInfoBranches    PFLepInfo ;

  GenInfo.Register(thisChain_) ; 
  EvtInfo.Register(thisChain_) ; 
  VtxInfo.Register(thisChain_) ; 
  LepInfo.Register(thisChain_) ; 
  LepInfo.Register(thisChain_,"PFLepInfo") ; 

  Double_t Weight3D[50][50][50];
  for(int ii=0;ii<50;ii++)
    for(int jj=0;jj<50;jj++)
      for(int kk=0;kk<50;kk++) 
        Weight3D[ii][jj][kk]=0.;      
  generate_weights(1, Weight3D);

  /**\ Looping over events */ 
  for (int iEntry=0;iEntry<thisChain_->GetEntries();++iEntry) { 

    /**\========print number of events done == == == == == == */  
    //if (iEntry%1000 == 0) std::cout << " processing " << iEntry << "th entry\n" ; 

    thisChain_->GetEntry(iEntry) ; 

    /**\ Select good vertices */
    for (int iVtx=0; iVtx < VtxInfo.Size; ++iVtx) { 
      if (VtxInfo.Type[iVtx]==1 && VtxInfo.isFake[iVtx]==false && VtxInfo.isValid[iVtx] && VtxInfo.Ndof[iVtx]>4 && VtxInfo.Rho[iVtx]<2. && VtxInfo.z[iVtx]<24.) { 
         ++nGoodVtxs ; 
      } 
    }

    histBeforeSelection_->hnGoodVtxs_->Fill(nGoodVtxs,weightPU) ; 

    if (isData_) {
      if (nGoodVtxs < 1)  continue ; 
    }

    /**\ Apply trigger */ 
    int TrgBit0 = HLT_DoubleEle15_SW_L1R ; 
    int TrgBit1 = HLT_DoubleEle17_SW_L1R ; 

    int npfEle(0) ;
    int ngsfEle(0) ;
            
    /**\ Select GSF electron */ 
    for (int iEle=0; iEle < LepInfo.Size; ++iEle) {
      if (    TMath::Abs(LepInfo.LeptonType[iEle])==11 
          &&  (TMath::Abs(LepInfo.simpleEleId80cIso[iEle]-7) < 0.1 || TMath::Abs(LepInfo.simpleEleId80cIso[iEle]-5) < 0.1)
          &&  TMath::Abs(LepInfo.ElTrackDxy_BS[iEle]) < 0.04
          &&  LepInfo.ChargeGsf[iEle] == LepInfo.ChargeCtf[iEle]
          &&  LepInfo.ChargeGsf[iEle] == LepInfo.ChargeScPix[iEle]
          &&  LepInfo.Pt[iEle] > 25. 
          &&  ( TMath::Abs(LepInfo.Eta[iEle])<1.4442 || (TMath::Abs(LepInfo.Eta[iEle])>1.566 && TMath::Abs(LepInfo.Eta[iEle])<2.4) ) 
	 ) {
         ++ngsfEle ;	      
      } // 1st ele accepted 
    } // loop over all ele: 1

    hngsfEle_->Fill(ngsfEle) ; 

    /**\ Select PF electron */ 
    for (int iEle=0; iEle < PFLepInfo.Size; ++iEle) {
      if (    TMath::Abs(PFLepInfo.LeptonType[iEle])==11 
          &&  (TMath::Abs(PFLepInfo.simpleEleId80cIso[iEle]-7) < 0.1 || TMath::Abs(PFLepInfo.simpleEleId80cIso[iEle]-5) < 0.1)
          &&  TMath::Abs(PFLepInfo.ElTrackDxy_BS[iEle]) < 0.04
          &&  PFLepInfo.ChargeGsf[iEle] == PFLepInfo.ChargeCtf[iEle]
          &&  PFLepInfo.ChargeGsf[iEle] == PFLepInfo.ChargeScPix[iEle]
          &&  PFLepInfo.TrackIso[iEle]/PFLepInfo.Et[iEle] < 0.2
          &&  PFLepInfo.EcalIso[iEle]/PFLepInfo.Et[iEle] < 0.2
          &&  PFLepInfo.HcalIso[iEle]/PFLepInfo.Et[iEle] < 0.2 
          &&  (PFLepInfo.ChargedHadronIso[iEle]+PFLepInfo.NeutralHadronIso[iEle]+PFLepInfo.PhotonIso[iEle])/TMath::Abs(PFLepInfo.Pt[iEle]) < 0.15 
          &&  PFLepInfo.Pt[iEle] > 25. 
          &&  ( TMath::Abs(PFLepInfo.Eta[iEle])<1.4442 || (TMath::Abs(PFLepInfo.Eta[iEle])>1.566 && TMath::Abs(PFLepInfo.Eta[iEle])<2.4) ) 
	 ) {
         ++npfEle ;	      
      } // 1st ele accepted 
    } // loop over all ele: 1

    hnpfEle_->Fill(npfEle) ; 

  } // event loop 

}

#endif 
