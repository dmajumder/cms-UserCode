#define READBPRIMEKITNTUPLE_CXX 
#include "/afs/cern.ch/work/d/devdatta/Analysis/FourthGen/CMSSW_5_3_6_patch1/src/BprimeKit/Headers/interface/format.h"
#include "/afs/cern.ch/work/d/devdatta/Analysis/FourthGen/CMSSW_5_3_6_patch1/src/BprimeKit/Headers/interface/TriggerBooking.h"
#include "/afs/cern.ch/work/d/devdatta/Analysis/FourthGen/CMSSW_5_3_6_patch1/src/BprimeKit/Headers/interface/SelectedJet.h" 
#include "/afs/cern.ch/work/d/devdatta/Analysis/FourthGen/CMSSW_5_3_6_patch1/src/BprimeKit/Headers/interface/JetCollection.h" 



#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <assert.h>
#include <TLorentzVector.h> 
#include <TChain.h> 
#include <TFile.h>
#include <TH1D.h>

using namespace std; 

const int doubleElectronHLTBitIndex[10] = {
	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5,
	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6,
	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7,
	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8,
	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9,
	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10, 
	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15, 
	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v16, 
	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17, 
	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v18  
} ;

const int dimuonHLTBitIndex[21] = {
	HLT_Mu17_Mu8_v2, 
	HLT_Mu17_Mu8_v3,  
	HLT_Mu17_Mu8_v4, 
	HLT_Mu17_Mu8_v6, 
	HLT_Mu17_Mu8_v7, 
	HLT_Mu17_Mu8_v10, 
	HLT_Mu17_Mu8_v11, 
	HLT_Mu17_Mu8_v16, 
	HLT_Mu17_Mu8_v17,  
	HLT_Mu17_Mu8_v18, 
	HLT_Mu17_Mu8_v19, 
	HLT_Mu17_Mu8_v21, 
	HLT_Mu17_Mu8_v22,  
	HLT_Mu17_TkMu8_v3, 
	HLT_Mu17_TkMu8_v4,
	HLT_Mu17_TkMu8_v9, 
	HLT_Mu17_TkMu8_v10, 
	HLT_Mu17_TkMu8_v12, 
	HLT_Mu17_TkMu8_v11, 
	HLT_Mu17_TkMu8_v13, 
	HLT_Mu17_TkMu8_v14 
} ;

const int MuEGHLTBitIndex[12] = { 
	HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4,
	HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5,
	HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6,
	HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7,
	HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8,
	HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9,
	HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4,
	HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5,
	HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6,
	HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7,
	HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8,
	HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9 
} ; 

static const int nmaxZ = 15 ; 

bool isData = 1 ;	
TChain* thisChain = new TChain("os2lTree"); 

//bool isGoodElectron (EvtInfoBranches*, LepInfoBranches*, int&, int&) ; 
//bool isGoodMuon (LepInfoBranches*,int& ,int& ) ; 
//bool isGoodTau (LepInfoBranches* ,int& ,int& ) ; 
void evtLoop () ; 

TFile* outfile ; 
TH1D* h_ptHiggs; 

TH1D* hscEta ; 
TH1D* hpfEta ; 

void readBprimeKitNtuple() {

	outfile = new TFile("MyFile.root","RECREATE") ; 
	outfile->cd() ; 
	h_ptHiggs = new TH1D("h_ptHiggs","Higgs pT at Gen level",1000,0.,1000.) ; 
	hscEta = new TH1D("hscEta","electron SC eta",200,-4.,4.) ; 
	hpfEta = new TH1D("hpfEta","electron PF eta",200,-4.,4.) ; 

	char inputFile[1000] = "/afs/cern.ch/work/d/devdatta/Analysis/FourthGen/CMSSW_5_3_6_patch1/src/MultileptonBprime/OS2LSelector/test/NTuples/SkimmedData/NTuple_OSOF_DoubleElectron_Run2012A-13Jul2012-v1_190456-193686.root" ; 
	ifstream infile ; 

	thisChain->Add(inputFile);
	std::cout << " chain name =" << thisChain->GetName() << std::endl ; 
	std::cout << " chain entries =" << thisChain->GetEntries() << std::endl ; 
	evtLoop();

	thisChain->Reset() ;

	outfile->Write() ; 

	return ; 

}

void evtLoop () {

	GenInfoBranches    GenInfo ; 
	EvtInfoBranches    EvtInfo ; 
	VertexInfoBranches VtxInfo ; 
	LepInfoBranches    LepInfo ; 
	JetInfoBranches    JetInfo ; 
	JetInfoBranches    WJetInfo ; 
	SelectedJet    selectedjet ;
	SelectedJet    jet_jesHigh ; 
	SelectedJet    jet_jesLow  ; 
	SelectedJet    jet_jerHigh ; 
	SelectedJet    jet_jerLow  ; 

	GenInfo.Register(thisChain) ; 
	EvtInfo.Register(thisChain) ; 
	VtxInfo.Register(thisChain) ; 
	LepInfo.Register(thisChain,"PFLepInfo") ; 
	JetInfo.Register(thisChain,"PFJetInfo") ; 
	WJetInfo.Register(thisChain,"WJetInfo") ;
	selectedjet.Register(thisChain) ; 
	jet_jesHigh.Register(thisChain,"Jets_JESHigh") ; 
	jet_jesLow.Register(thisChain,"Jets_JESLow") ; 
	jet_jerHigh.Register(thisChain,"Jets_JERHigh") ; 
	jet_jerLow.Register(thisChain,"Jets_JERLow") ; 

	//// counters for signal: Bprime decay modes:
	int nbZbZ(0) ; 
	int nbZtW(0) ; 
	int ntWtW(0) ; 
	int nbHbH(0) ; 
	int nbHbZ(0) ; 
	int nbHtW(0) ; 

	/**\ Event selection flags */
	int  nGoodVtxs(0) ; 

	/**\ Looping over events */ 
	//for (int iEntry=0;iEntry<thisChain->GetEntries();++iEntry) { 
	for (int iEntry = 0; iEntry <  20000; ++iEntry) { 

		thisChain->GetEntry(iEntry) ; 

		//// Pick Higgs bosons 
		double pth0(0); 
		int indexh0(0); 
		for (int ii=0; ii < GenInfo.Size; ++ii) {
			if (GenInfo.Status[ii] == 3 && TMath::Abs(GenInfo.PdgID[ii]) == 25 ) {
				if (GenInfo.Pt[ii] > pth0) { pth0 = GenInfo.Pt[ii] ; indexh0 = ii ; }
			}
		}
		//// Fill leading Higgs boson pT
		h_ptHiggs->Fill(GenInfo.Pt[indexh0]) ; 


		//// Get Bprime decay modes:
		if (!isData) {
			if ( (EvtInfo.McbprimeMode[0] == 1 && EvtInfo.McbprimeMode[1] == 3) || (EvtInfo.McbprimeMode[0] == 3 && EvtInfo.McbprimeMode[1] == 1) ) ++nbZtW ; 
			if (EvtInfo.McbprimeMode[0] == 1 && EvtInfo.McbprimeMode[1] == 1) ++ntWtW ; 
			if (EvtInfo.McbprimeMode[0] == 3 && EvtInfo.McbprimeMode[1] == 3) ++nbZbZ ; 
			if (EvtInfo.McbprimeMode[0] == 4 && EvtInfo.McbprimeMode[1] == 4) ++nbHbH ; 
			if ( (EvtInfo.McbprimeMode[0] == 3 && EvtInfo.McbprimeMode[1] == 4) || (EvtInfo.McbprimeMode[0] == 4 && EvtInfo.McbprimeMode[1] == 3) ) ++nbHbZ ; 
			if ( (EvtInfo.McbprimeMode[0] == 1 && EvtInfo.McbprimeMode[1] == 4) || (EvtInfo.McbprimeMode[0] == 4 && EvtInfo.McbprimeMode[1] == 1) ) ++nbHtW ; 
		}

		/**\ Appy HLT requirement */
		bool isPassHLT(false) ;
		for (int iHLT = 0; iHLT < 21; ++iHLT) {
			if (EvtInfo.TrgBook[dimuonHLTBitIndex[iHLT]] == 1 ) {
				isPassHLT = true ;
				break ; 
			}
			else { isPassHLT = false ; } 
		}

		if (isPassHLT == true) std::cout << " Event " << iEntry << " passes HLT\n" ; 

		nGoodVtxs = 0 ; 
		/**\ Select good vertices */
		for (int iVtx=0; iVtx < VtxInfo.Size; ++iVtx) { 
			if (   VtxInfo.Type[iVtx]==1 
					&& VtxInfo.isFake[iVtx]==false 
					&& VtxInfo.Ndof[iVtx]>4 
					&& VtxInfo.Rho[iVtx]<2. 
					&& VtxInfo.z[iVtx]<24.) { ++nGoodVtxs ; }  
		}
		if (nGoodVtxs < 1)  { continue ; } 

		EvtInfoBranches* evtinfo = &EvtInfo ; 
		LepInfoBranches* lepInfo = &LepInfo ; 

		for (int iLep = 0; iLep < LepInfo.Size; ++iLep) {
	                if (LepInfo.EgammaCutBasedEleIdMEDIUM[iLep] == true) { 
				TLorentzVector ele14v (LepInfo.Px[iLep],LepInfo.Py[iLep],LepInfo.Pz[iLep],LepInfo.Energy[iLep]) ; 
				hscEta->Fill(LepInfo.Eta[iLep]) ; 
				hpfEta->Fill(ele14v.Eta()) ; 
			}
		}

		/**\ Select fat jets reconstructed using the CA8 algorithm */ 

		for (int iWJet = 0; iWJet < WJetInfo.Size; ++iWJet) {
			//if (WJetInfo.CombinedSVBJetTags[iWJet] > .244) std::cout << " This jet passes CSVL\n" ;
			//if (WJetInfo.CombinedSVBJetTags[iWJet] > .679) std::cout << " This jet passes CSVM\n" ;
			//if (WJetInfo.CombinedSVBJetTags[iWJet] > .898) std::cout << " This jet passes CSVT\n" ;
		}

		/**\ Select jets */

		int nbtags(0) ; 

		for (int iJet=0;iJet < JetInfo.Size;++iJet) {

			if (JetInfo.Pt[iJet] < 30. ||
					TMath::Abs(JetInfo.Eta[iJet]) > 2.4 || 
					JetInfo.CHF[iJet] <= 0.00 ||
					JetInfo.NHF[iJet] >= 0.99 ||
					JetInfo.CEF[iJet] >= 0.99 ||
					JetInfo.NEF[iJet] >= 0.99 ||
					JetInfo.NCH[iJet] <= 0.00  
			   ) {
				continue; 
			}  // apply jet ID 

			//// This jet passes jet ID. Can be used for analysis 

		}

	} // End of event loop 

	std::cout << " Events processed = " << thisChain->GetEntries() << std::endl ; 

} // evtLoop () 

//#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h" 
//#include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h" 
//
//bool isGoodElectron (EvtInfoBranches* EvtInfo, LepInfoBranches* LepInfo, int& entry, int& ele) { 
//
//	if (LepInfo->EgammaCutBasedEleIdMEDIUM[ele] == true) { return true ; } 
//	else { return false ; }
//
//	if ( TMath::Abs(LepInfo->LeptonType[ele]) != 11 ) return false ;
//	if ( LepInfo->Pt[ele] < 30. ) return false ;
//
//	int idx = -1;
//	if((fabs(LepInfo->Eta[ele])<=1.4442)){
//		idx = 0;
//	} else if(((fabs(LepInfo->Eta[ele])<2.4)&&(fabs(LepInfo->Eta[ele])>=1.566))){ 
//		idx = 1;
//	}
//
//	if( idx == -1 ) return false ;
//
//	float        cut_dEtaIn[2]         = {999.9, 999.9};
//	float        cut_dPhiIn[2]         = {999.9, 999.9};
//	float        cut_sigmaIEtaIEta[2]  = {999.9, 999.9};
//	float        cut_hoe[2]            = {999.9, 999.9};
//	float        cut_ooemoop[2]        = {999.9, 999.9};
//	float        cut_d0vtx[2]          = {999.9, 999.9};
//	float        cut_dzvtx[2]          = {999.9, 999.9};
//	float        cut_iso[2]            = {999.9, 999.9};
//	bool         cut_vtxFit[2]         = {false, false};
//	unsigned int cut_mHits[2]          = {999,   999  };
//
//	cut_dEtaIn[0]        = 0.004; cut_dEtaIn[1]        = 0.007;
//	cut_dPhiIn[0]        = 0.060; cut_dPhiIn[1]        = 0.030;
//	cut_sigmaIEtaIEta[0] = 0.010; cut_sigmaIEtaIEta[1] = 0.030;
//	cut_hoe[0]           = 0.120; cut_hoe[1]           = 0.100;
//	cut_ooemoop[0]       = 0.050; cut_ooemoop[1]       = 0.050;
//	cut_d0vtx[0]         = 0.020; cut_d0vtx[1]         = 0.020;
//	cut_dzvtx[0]         = 0.100; cut_dzvtx[1]         = 0.100;
//	cut_iso[0]           = 0.150; cut_iso[1]           = 0.150;
//	cut_vtxFit[0]        = true ; cut_vtxFit[1]        = true;
//	cut_mHits[0]         = 1    ; cut_mHits[1]         = 1;
//
//	const float        pt                = LepInfo->Pt[ele] ;
//	const float        dEtaIn            = LepInfo->EldeltaEta[ele] ; 
//	const float        dPhiIn            = LepInfo->EldeltaPhi[ele] ; 
//	const float        sigmaIEtaIEta     = LepInfo->ElsigmaIetaIeta[ele] ; 
//	const float        hoe               = LepInfo->ElHadoverEm[ele] ; 
//	const float        ooemoop           = ((1.0/LepInfo->ElEcalE[ele]) - (LepInfo->ElEoverP[ele]/LepInfo->ElEcalE[ele])) ; 
//	const float        d0vtx             = LepInfo->ElTrackDxy_PV[ele] ; 
//	const float        dzvtx             = LepInfo->ElTrackDz[ele] ; 
//	const float        iso_ch            = LepInfo->ChargedHadronIso[ele] ; 
//	const float        iso_em            = LepInfo->PhotonIso[ele] ;
//	const float        iso_nh            = LepInfo->NeutralHadronIso[ele] ; 
//	const float        isorhocorr03      = LepInfo->IsoRhoCorrR03[ele] ; 
//	const bool         vtxFitConversion  = LepInfo->ElhasConv[ele] ; 
//	const unsigned int mHits             = LepInfo->NumberOfExpectedInnerHits[ele] ; 
//	const double       rho               = std::max(double(EvtInfo->RhoPU[0]),0.0); ; 
//
//	// effective area for isolation
//	float AEff = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, LepInfo->Eta[ele], ElectronEffectiveArea::kEleEAData2012);
//
//	// apply to neutrals 
//	double iso_n = std::max((iso_nh + iso_em) - (rho * AEff), 0.0); 
//
//	// compute final isolation 
//	double iso = (iso_n + iso_ch) / pt;
//
//	bool isGoodElectronCutbased(false) ; 
//
//	// test cuts
//	if ( (fabs(dEtaIn) < cut_dEtaIn[idx])             
//			&& (fabs(dPhiIn) < cut_dPhiIn[idx])             
//			&& (sigmaIEtaIEta < cut_sigmaIEtaIEta[idx])     
//			&& (hoe < cut_hoe[idx])                         
//			&& (fabs(ooemoop) < cut_ooemoop[idx])           
//			&& (fabs(d0vtx) < cut_d0vtx[idx])               
//			&& (fabs(dzvtx) < cut_dzvtx[idx])               
//			&& (!cut_vtxFit[idx] || !vtxFitConversion)      
//			&& (mHits <= cut_mHits[idx])                    
//			&& (iso < cut_iso[idx])                         
//			&& ( (LepInfo->Charge[ele] == 1 && LepInfo->ChargeCtf[ele] == 1 && LepInfo->ChargeGsf[ele] == 1) || 
//				(LepInfo->Charge[ele] == -1 && LepInfo->ChargeCtf[ele] == -1 && LepInfo->ChargeGsf[ele] == -1) )  
//	   ) { isGoodElectronCutbased = true ; } 
//
//	return isGoodElectronCutbased ; 
//
//}
//
//bool  isGoodMuon (LepInfoBranches* LepInfo,int& entry,int& mu) { 
//
//	bool isGoodMuon(false) ; 
//
//	if(   TMath::Abs(LepInfo->LeptonType[mu]) == 13 
//			&& LepInfo->isPFMuon[mu] == true 
//			&& (LepInfo->MuType[mu] & (1<<5)) != 0 
//			&& LepInfo->innerTracknormalizedChi2[mu] < 10 
//			&& LepInfo->MuNMuonhits[mu] > 0
//			&& LepInfo->MuNMatchedStations[mu] > 1
//			&& TMath::Abs(LepInfo->MuInnerTrackDxy_PV[mu]) < 0.2
//			&& TMath::Abs(LepInfo->MuInnerTrackDz[mu]) < 0.5
//			&& LepInfo->MuNPixelLayers[mu] > 0
//			&& LepInfo->MuNTrackLayersWMeasurement[mu] > 5  
//			&& LepInfo->Pt[mu] > 30.
//			&& TMath::Abs(LepInfo->Eta[mu]) < 2.4 
//			&& (( LepInfo->ChargedHadronIsoR04[mu] + max(LepInfo->NeutralHadronIsoR04[mu] + LepInfo->PhotonIsoR04[mu] - 0.5*LepInfo->sumPUPtR04[mu], 0.0))/LepInfo->Pt[mu]<0.15) 
//	  ) {
//		isGoodMuon = true ; 	   
//	} else isGoodMuon = false ; 
//
//	return isGoodMuon ; 	
//
//}
//
//bool  isGoodTau (LepInfoBranches* LepInfo,int& entry,int& tau) {
//	bool isGoodTau(false) ;
//	if ( LepInfo->isPFTau[tau] 
//			//&& LepInfo->againstElectronLoose[tau] > .5 
//			&& LepInfo->againstElectronMVA[tau] > .5 
//			&& LepInfo->againstMuonLoose[tau] > .5 
//			&& LepInfo->Pt[tau] > 20.
//			&& LepInfo->byLooseCombinedIsolationDeltaBetaCorr[tau] == 1) isGoodTau = true ; 
//	else isGoodTau = false ; 
//	return isGoodTau ; 
//}

