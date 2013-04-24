#include "MultileptonBprime/OS2LSelector/interface/OS2LSelector.h"
#include "BprimeKit/Headers/interface/format.h" 
#include "BprimeKit/Headers/interface/TriggerBooking.h"
#include "BprimeKit/Utils/interface/checkEvt.h"
#include "BprimeKit/Utils/interface/getPUWeight.h"
#include "BprimeKit/Utils/interface/Histograms.h" 
//#include "BprimeKit/Utils/interface/mathutils.h"

#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <assert.h>
#include <TLorentzVector.h> 

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

const double massZboson = 91.1876 ; // GeV/c2 

/**\ Constructor */
OS2LSelector::OS2LSelector (const std::string&infiles,const bool&isData,const std::string&jsonfile,const std::string&outfile) :
  inputfiles_ (infiles),
  outfileName_((char*)outfile.c_str()), 
  outfile_ ((char*)outfile.c_str(),"RECREATE"),
  isData_ (isData),
  jsonfile_(jsonfile), 
  thisChain_ (0),   	
  evtProc_ (0), 
  evtAccept_ (0), 
  ngoodElectrons_ (0), 
  ngoodMuons_ (0) 
{
  std::cout << " Creating this instance of OS2LSelector\n " ; 

  /**\ Initialize histograms */
  outfile_.cd() ; 
  hist_ = new Histograms() ;
  hist_->setSumw2() ; 
  m_hcutFlow_           = new TH1D("hCutFlow","Cut flow",25,.5,25.5) ; 
  m_hnvtx_PUWt_         = new TH1D("hnvtx_PUWt","N(vtx) with PU reweighting",100,0.,100.) ; 
  m_hnvtx_noPUWt_       = new TH1D("hnvtx_noPUWt","N(vtx) without PU reweighting",100,0.,100.) ; 
  m_hnvtx_PUWt_ZWin_    = new TH1D("hnvtx_PUWt_ZWin","N(vtx) with PU reweighting",100,0.,100.) ; 
  m_hnvtx_noPUWt_ZWin_  = new TH1D("hnvtx_noPUWt_ZWin","N(vtx) without PU reweighting",100,0.,100.) ; 
  m_hnvtx_PUWt_ZVeto_   = new TH1D("hnvtx_PUWt_ZVeto","N(vtx) with PU reweighting",100,0.,100.) ; 
  m_hnvtx_noPUWt_ZVeto_ = new TH1D("hnvtx_noPUWt_ZVeto","N(vtx) without PU reweighting",100,0.,100.) ; 
  m_hpte_               = new TH1D("hpte","electron pt",1000,0.,1000.) ; 
  m_hetae_              = new TH1D("hetae","electron eta",1000,-4.,4.) ; 
  m_hptee_              = new TH1D("hptee","p_{T}(e^{+}e^{-})",1000,0.,1000.) ; 
  m_hetaee_             = new TH1D("hetaee","#eta(e^{+}e^{-}))",1000,-4.,4.) ; 
  m_hnjets_             = new TH1D("hnjets","Njets",21,-.5,20.5) ; 
  m_hnbtags_            = new TH1D("hnbtags","Nbtags",21,-.5,20.5) ; 

}

 /**\ Destructor */ 
OS2LSelector::~OS2LSelector () {
	std::cout << "OS2LSelector::~OS2LSelector\n" ; 
}

void OS2LSelector::process () {

  char inputFile[1000] ; 
  ifstream infile ; 
  infile.open ((char*)inputfiles_.c_str(),ifstream::in) ; 
  assert(!infile.fail()) ;

  while (infile.good()) { 

	  infile >> inputFile ; 
	  if(strncmp(inputFile,"#",1)==0) continue; 
	  if(strncmp(inputFile,"%",1)==0) break ; 

	  if ( strstr(inputFile,"Reduced")!=0 ) {
		  thisChain_ = new TChain("root") ;   
	  }
	  else if (strstr(inputFile,"filter2L")!=0) {
		  thisChain_ = new TChain("bprimeKit/root") ; 	    
	  }

	  thisChain_->Add(inputFile); 
	  std::cout << "\ninfile " << inputFile << std::endl ; 
  }
  infile.close() ;   
  evtProc_ = thisChain_->GetEntries(); 
  evtLoop();
  thisChain_->Reset() ;

  outfile_.Write() ; 
  outfile_.Close() ; 

  std::cout << " outfile written\n" ; 

  return ; 

}

#include "BprimeKit/Headers/interface/JetCollection.h" 

void::OS2LSelector::evtLoop () {

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

	int    nossf_el(0)      ; 
	int    nossf_mu(0)      ; 
	int    nosof(0)         ; 

	double weightPU0(1)     ; // Event flag 
	double weightPUP(1)     ; // Event flag 
	double weightPUM(1)     ; // Event flag 
	int    nGoodVtxs(0)     ; // Event flag  
	int    nGoodElp(0)      ; // Event flag  
	int    nGoodElm(0)      ; // Event flag  
	int    nGoodMup(0)      ; // Event flag  
	int    nGoodMum(0)      ; // Event flag  
	int    elPIndex[nmaxZ]  ; // Event flag  
	int    elMIndex[nmaxZ]  ; // Event flag  
	int    muPIndex[nmaxZ]  ; // Event flag  
	int    muMIndex[nmaxZ]  ; // Event flag  
	bool   isOSSF_el(0)     ; // Event flag  
	bool   isOSSF_mu(0)     ; // Event flag  
	bool   isOSOF(0)        ; // Event flag  
	bool   isZVeto(0)       ; // Event flag 
	bool   isZVetoLow(0)    ; // Event flag 
	bool   isZVetoHigh(0)   ; // Event flag 
	SelectedJet selectedjet ; // Event flag 
	SelectedJet jet_jesHigh ; // Event flag  
	SelectedJet jet_jesLow  ; // Event flag  
	SelectedJet jet_jerHigh ; // Event flag  
	SelectedJet jet_jerLow  ; // Event flag  

	TString name = "NTuple_" ;
	name += outfileName_ ; 
	TFile* fntuple = new TFile(name,"RECREATE") ; 
	fntuple->cd() ; 
	TTree* os2lTree = thisChain_->CloneTree(0); 
	os2lTree->SetName("os2lTree") ; 
	os2lTree->Branch("weightPU0",   &weightPU0,   "weightPU0/D"); 
	os2lTree->Branch("weightPUP",   &weightPUP,   "weightPUP/D");
	os2lTree->Branch("weightPUM",   &weightPUM,   "weightPUM/D"); 
	os2lTree->Branch("nGoodVtxs",   &nGoodVtxs,   "nGoodVtxs/D"); 
	//os2lTree->Branch("nGoodElp",    &nGoodElp,    "nGoodElp/I") ; 
	//os2lTree->Branch("nGoodElm",    &nGoodElm,    "nGoodElm/I") ; 
	//os2lTree->Branch("nGoodMup",    &nGoodMup,    "nGoodMup/I") ; 
	//os2lTree->Branch("nGoodMum",    &nGoodMum,    "nGoodMum/I") ; 
	os2lTree->Branch("elPIndex",    &elPIndex[0], "elPIndex/I") ;
	os2lTree->Branch("elMIndex",    &elMIndex[0], "elMIndex/I") ; 
	os2lTree->Branch("muPIndex",    &muPIndex[0], "muPIndex/I") ;
	os2lTree->Branch("muMIndex",    &muMIndex[0], "muMIndex/I") ; 
	os2lTree->Branch("isOSSF_el",   &isOSSF_el,   "isOSSF_el/O"); 
	os2lTree->Branch("isOSSF_mu",   &isOSSF_mu,   "isOSSF_mu/O"); 
	os2lTree->Branch("isOSOF",      &isOSOF,      "isOSOF/O"); 
	os2lTree->Branch("isZVeto",     &isZVeto,     "isZVeto/O"); 
	os2lTree->Branch("isZVetoLow",  &isZVetoLow,  "isZVetoLow/O"); 
	os2lTree->Branch("isZVetoHigh", &isZVetoHigh, "isZVetoHight/O"); 
	selectedjet.RegisterTree(os2lTree) ; 
	jet_jesHigh.RegisterTree(os2lTree,"Jets_JESHigh") ;  
	jet_jesLow.RegisterTree(os2lTree,"Jets_JESLow") ;  
	jet_jerHigh.RegisterTree(os2lTree,"Jets_JERHigh") ;  
	jet_jerLow.RegisterTree(os2lTree,"Jets_JERLow") ;  

	m_hcutFlow_->SetBinContent(1,evtProc_) ; 

	if (isData_) {
		MakeJsonMap(jsonfile_);
	}

	std::cout << " Events processed = " << thisChain_->GetEntries() << std::endl ; 

	/**\ Looping over events */ 
	for (int iEntry = 0; iEntry < thisChain_->GetEntries(); ++iEntry) { 

		weightPU0 = 1. ;
		weightPUP = 1. ;
		weightPUM = 1. ;

		thisChain_->GetEntry(iEntry) ; 

		/**\ Get PU weight of event for MC */
		if (!isData_) {
			int n0 = -1; 
			for(int i = 0; i < EvtInfo.nBX; i++) { 
				if (EvtInfo.BXPU[i] ==  0) { n0  = EvtInfo.nPU[i] ; }
			}
			std::string mode = "electron"; 
			double weights[3] = {1,1,1} ; 
			getPUWeight(mode,n0,weights); 
			weightPUP = weights[0] ; 
			weightPU0 = weights[1] ; 
			weightPUM = weights[2] ; 
		} 
		else { weightPUP = 1; weightPU0 = 1; weightPUM = 1 ;} 

		/**\ Appy HLT requirement */
		bool passHLT_DoubleE(false) ;
		for (int iHLT = 0; iHLT < 10; ++iHLT) {
			if (EvtInfo.TrgBook[doubleElectronHLTBitIndex[iHLT]] == 1 ) {
				passHLT_DoubleE = true ;
				break ; 
			}
			else { passHLT_DoubleE = false ; } 
		}

		bool passHLT_DoubleMu(false) ;
		for (int iHLT = 0; iHLT < 21; ++iHLT) {
			if (EvtInfo.TrgBook[dimuonHLTBitIndex[iHLT]] == 1 ) {
				passHLT_DoubleMu = true ;
				break ; 
			}
			else { passHLT_DoubleMu = false ; } 
		}

		bool passHLT_MuEG(false) ;
		for (int iHLT = 0; iHLT < 12; ++iHLT) {
			if (EvtInfo.TrgBook[MuEGHLTBitIndex[iHLT]] == 1 ) {
				passHLT_MuEG = true ;
				break ; 
			}
			else { passHLT_MuEG = false ; } 
		}

		//if (isData_) 
		if (passHLT_DoubleE == false && passHLT_DoubleMu == false && passHLT_MuEG == false) continue ; 

		nGoodVtxs    = 0     ; 
	        nGoodElp     = 0 ;  
	        nGoodElm     = 0 ;  
	        nGoodMup     = 0 ;  
	        nGoodMum     = 0 ;  
		for (int ii = 0; ii < nmaxZ; ++ii) {
			elPIndex[ii] = -1 ; 
			elMIndex[ii] = -1 ; 
			muPIndex[ii] = -1 ; 
			muMIndex[ii] = -1 ; 
		}
		isOSSF_el    = false ; 
		isOSSF_mu    = false ; 
		isOSOF       = false ; 
		isZVeto      = true  ; 
		isZVetoLow   = true  ; 
		isZVetoHigh  = true  ; 

		int ntaus(0) ;

		if(isData_ && !isGoodEvt(EvtInfo.RunNo,EvtInfo.LumiNo)) continue;    

		/**\ Select good vertices */
		for (int iVtx=0; iVtx < VtxInfo.Size; ++iVtx) { 
			if (   VtxInfo.Type[iVtx]==1 
					&& VtxInfo.isFake[iVtx]==false 
					&& VtxInfo.Ndof[iVtx]>4 
					&& VtxInfo.Rho[iVtx]<2. 
					&& VtxInfo.z[iVtx]<24.) { ++nGoodVtxs ; }  
		}


		hist_->hnGoodVtxsNoPUWt_->Fill(nGoodVtxs, 1.) ;  
		hist_->hnGoodVtxs_->Fill(nGoodVtxs, weightPU0) ; 

		if (nGoodVtxs < 1)  { continue ; } 

		EvtInfoBranches* evtinfo = &EvtInfo ; 
		LepInfoBranches* lepInfo = &LepInfo ; 


		for (int iLep=0; iLep < LepInfo.Size; ++iLep) {

			//if (isGoodElectron(evtinfo,lepInfo,iEntry,iLep)) {
			//	if (LepInfo.ChargeGsf[iLep] > 0) { elPIndex[nGoodElp] = iLep ;++nGoodElp; } 
			//	else if (LepInfo.ChargeGsf[iLep] < 0) { elMIndex[nGoodElm] = iLep ;  ++nGoodElm; } 
			//}

			//if (isGoodMuon(lepInfo,iEntry,iLep)) {
			//	if (LepInfo.Charge[iLep] > 0) { muPIndex[nGoodMup] = iLep ;++nGoodMup; }
			//	else if (LepInfo.Charge[iLep] < 0) { muMIndex[nGoodMum] = iLep ;++nGoodMum; }
			//}

			if (isGoodTau(lepInfo,iEntry,iLep)) { 
				++ntaus ; // Good tau counting 
				continue ; 
			}

			if (isGoodElectron(evtinfo,lepInfo,iEntry,iLep)) {
				if (LepInfo.ChargeGsf[iLep] > 0) { elPIndex[nGoodElp] = iLep ;++nGoodElp; } 
				else if (LepInfo.ChargeGsf[iLep] < 0) { elMIndex[nGoodElm] = iLep ;  ++nGoodElm; } 
				for (int jLep=iLep+1; jLep < LepInfo.Size; ++jLep) { 
					if(    LepInfo.ChargeGsf[iLep]*LepInfo.ChargeGsf[jLep] < 0
							&& isGoodElectron(evtinfo,lepInfo,iEntry,jLep) ) { 
						if (LepInfo.ChargeGsf[jLep] > 0) { elPIndex[nGoodElp] = jLep ; } 
						else if (LepInfo.ChargeGsf[jLep] < 0) { elMIndex[nGoodElm] = jLep ;  } 
						TLorentzVector ele14v (LepInfo.Px[iLep],LepInfo.Py[iLep],LepInfo.Pz[iLep],LepInfo.Energy[iLep]) ; 
						TLorentzVector ele24v (LepInfo.Px[jLep],LepInfo.Py[jLep],LepInfo.Pz[jLep],LepInfo.Energy[jLep]) ; 
						TLorentzVector zcand4v (ele14v+ele24v) ;
						hist_->hmee_->Fill(zcand4v.Mag(),weightPU0) ; 
						m_hptee_->Fill(zcand4v.Pt(),weightPU0) ; 
						hist_->hZPt_->Fill(zcand4v.Pt(),weightPU0) ; 
						m_hetaee_->Fill(zcand4v.PseudoRapidity(),weightPU0) ; 
						if (zcand4v.Mag() > 60. && zcand4v.Mag() < 120.) {
							isZVeto = false ; 	
							isZVetoLow = false; 	    
							isZVetoHigh = false ; 	    
						} 
						else {
							isZVeto = true ; 
							if (zcand4v.Mag() < 60.) {
								isZVetoLow = true ; 	    
								isZVetoHigh = false ; 
							} 
							else if (zcand4v.Mag() > 120.) { 
								isZVetoLow = false ; 	    
								isZVetoHigh = true ; 	    
							}
							else {
								std::cout << ">>>>>Warning: M(l+l-) = " << zcand4v.Mag() << " not in range\n" ; 
							}
						}
					} // Select 2nd good electron 
				}
				continue ; 
			} // Select 1st good electron 

			if (isGoodMuon(lepInfo,iEntry,iLep)) {
				if (LepInfo.Charge[iLep] > 0) { muPIndex[nGoodMup] = iLep ;++nGoodMup; }
				else if (LepInfo.Charge[iLep] < 0) { muMIndex[nGoodMum] = iLep ;++nGoodMum; }
				for (int jLep=iLep+1; jLep < LepInfo.Size; ++jLep) { 
					if(    LepInfo.Charge[iLep]*LepInfo.Charge[jLep] < 0
							&& isGoodMuon(lepInfo,iEntry,jLep) ) { 
						if (LepInfo.Charge[jLep] > 0) { muPIndex[nGoodMup] = jLep ; }
						else if (LepInfo.Charge[jLep] < 0) { muMIndex[nGoodMum] = jLep ; }
						TLorentzVector mu14v (LepInfo.Px[iLep],LepInfo.Py[iLep],LepInfo.Pz[iLep],LepInfo.Energy[iLep]) ; 
						TLorentzVector mu24v (LepInfo.Px[jLep],LepInfo.Py[jLep],LepInfo.Pz[jLep],LepInfo.Energy[jLep]) ; 
						TLorentzVector zcand4v (mu14v+mu24v) ;
						if (zcand4v.Mag() > 60. && zcand4v.Mag() < 120.) {
							isZVeto = false ; 
							isZVetoLow = false; 	    
							isZVetoHigh = false ; 	    
						} 
						else {
							isZVeto = true ; 
							if (zcand4v.Mag() < 60.) {
								isZVetoLow = true ; 	    
								isZVetoHigh = false ; 	    
							} 
							else if (zcand4v.Mag() > 120.) { 
								isZVetoLow = false; 	    
								isZVetoHigh = true ; 	    
							}
							else {
								std::cout << ">>>>>Warning: M(l+l-) = " << zcand4v.Mag() << " not in range\n" ; 
							}
						}
					} // Select 2nd good muon 
				} 
				continue ; 
			} // Select 1st good muon 

		} // Looping over leptons 

		if (ntaus > 0) continue ; // tau-veto applied 

		if ( (passHLT_DoubleMu && !passHLT_DoubleE  && !passHLT_MuEG   ) && ( nGoodMup == 1 && nGoodMum == 1 && nGoodElp == 0 && nGoodElm == 0 ) )  { isOSSF_mu = true ; ++nossf_mu ; }
		if ( (passHLT_DoubleE  && !passHLT_DoubleMu && !passHLT_MuEG   ) && ( nGoodMup == 0 && nGoodMum == 0 && nGoodElp == 1 && nGoodElm == 1 ) )  { isOSSF_el = true ; ++nossf_el ; } 
		if ( (passHLT_MuEG     && !passHLT_DoubleMu && !passHLT_DoubleE) && ( (nGoodMup == 1 && nGoodElm == 1 && nGoodElp == 0 && nGoodMum == 0) || (nGoodMup == 0 && nGoodElm == 0 && nGoodElp == 1 && nGoodMum == 1) ) ) { isOSOF = true ; ++nosof ; } 

		//// Make ZeeCands 
		//TLorentzVector eleP4v  ;
		//TLorentzVector eleM4v  ;
		//TLorentzVector zee4v ; 
		//double mdiff(massZboson) ; 
		//int theElPIndex, theElMIndex ; 
		//for (int ii = 0; ii < nmaxZ; ++ii) {
		//	eleP4v.SetPxPyPzE(LepInfo.Px[elPIndex[ii]],LepInfo.Py[elPIndex[ii]],LepInfo.Pz[elPIndex[ii]],LepInfo.Energy[elPIndex[ii]]) ;
		//	for (int jj = 0; jj < nmaxZ; ++jj) {
		//		eleM4v.SetPxPyPzE(LepInfo.Px[elMIndex[jj]],LepInfo.Py[elMIndex[jj]],LepInfo.Pz[elMIndex[jj]],LepInfo.Energy[elMIndex[jj]]) ;
		//		zee4v = eleP4v  + eleM4v ; 
		//		if ( TMath::Abs(zee4v.Mag() - massZboson) < mdiff  ) {
		//			mdiff = TMath::Abs(zee4v.Mag() - massZboson);  
		//			theElPIndex = elPIndex[ii] ; 
		//			theElMIndex = elMIndex[jj] ; 
		//		}
		//	}
		//}
		//eleP4v.SetPxPyPzE(LepInfo.Px[theElPIndex],LepInfo.Py[theElPIndex],LepInfo.Pz[theElPIndex],LepInfo.Energy[theElPIndex]) ;
		//eleM4v.SetPxPyPzE(LepInfo.Px[theElMIndex],LepInfo.Py[theElMIndex],LepInfo.Pz[theElMIndex],LepInfo.Energy[theElMIndex]) ;
		//zee4v = eleP4v  + eleM4v ; 
		//if (zee4v.Mag() > 76 && zee4v.Mag() < 106) { 
		//	hasZeeCand = true ; 
		//}
		//else {
		//	isZVeto = true ; 
		//	if (zee4v.Mag() < 76 ) {
		//		isZVetoLow = true ; 
		//	}
		//	else if (zee4v.Mag() > 106) {
		//		isZVetoHigh = true ; 
		//	}
		//}
		//mdiff = massZboson ; 

		//// Make ZmumuCands 
		//TLorentzVector muP4v  ;
		//TLorentzVector muM4v  ;
		//TLorentzVector zmumu4v ; 
		//int theMuPIndex, theMuMIndex ; 
		//for (int ii = 0; ii < nmaxZ; ++ii) {
		//	muP4v.SetPxPyPzE(LepInfo.Px[muPIndex[ii]],LepInfo.Py[muPIndex[ii]],LepInfo.Pz[muPIndex[ii]],LepInfo.Energy[muPIndex[ii]]) ;
		//	for (int jj = 0; jj < nmaxZ; ++jj) {
		//		muM4v.SetPxPyPzE(LepInfo.Px[muMIndex[jj]],LepInfo.Py[muMIndex[jj]],LepInfo.Pz[muMIndex[jj]],LepInfo.Energy[muMIndex[jj]]) ;
		//		zmumu4v = muP4v  + muM4v ; 
		//		if ( TMath::Abs(zmumu4v.Mag() - massZboson) < mdiff  ) {
		//			mdiff = TMath::Abs(zmumu4v.Mag() - massZboson);  
		//			theMuPIndex = muPIndex[ii] ; 
		//			theMuMIndex = muMIndex[jj] ; 
		//		}
		//	}
		//}
		//muP4v.SetPxPyPzE(LepInfo.Px[theMuPIndex],LepInfo.Py[theMuPIndex],LepInfo.Pz[theMuPIndex],LepInfo.Energy[theMuPIndex]) ;
		//muM4v.SetPxPyPzE(LepInfo.Px[theMuMIndex],LepInfo.Py[theMuMIndex],LepInfo.Pz[theMuMIndex],LepInfo.Energy[theMuMIndex]) ;
		//zmumu4v = muP4v  + muM4v ; 
		//if (zmumu4v.Mag() > 76 && zmumu4v.Mag() < 106) { 
		//	hasZeeCand = true ; 
		//}
		//else {
		//	isZVeto = true ; 
		//	if (zmumu4v.Mag() < 76 ) {
		//		isZVetoLow = true ; 
		//	}
		//	else if (zmumu4v.Mag() > 106) {
		//		isZVetoHigh = true ; 
		//	}
		//}

		if (isOSSF_mu == false && isOSSF_el == false && isOSOF == false) continue ; //// selecting OS2L events 

		m_hnvtx_PUWt_->Fill(nGoodVtxs, weightPU0) ; 
		m_hnvtx_noPUWt_->Fill(nGoodVtxs, 1) ; 

		if (isZVeto) {
			m_hnvtx_PUWt_ZVeto_->Fill(nGoodVtxs, weightPU0) ; 
			m_hnvtx_noPUWt_ZVeto_->Fill(nGoodVtxs, 1) ; 
		}
		else {
			m_hnvtx_PUWt_ZWin_->Fill(nGoodVtxs, weightPU0) ; 
			m_hnvtx_noPUWt_ZWin_->Fill(nGoodVtxs, 1) ; 
		}

		/**\ Select jets */

		selectedjet.Size = 0 ; 
		int nbtags(0) ; 

		for (int iJet=0;iJet < JetInfo.Size;++iJet) {

			bool jetLepOverlap(false) ; 
			for (int iLep=0; iLep < LepInfo.Size; ++iLep) {
				if(isGoodMuon(lepInfo,iEntry,iLep) || isGoodElectron(evtinfo,lepInfo,iEntry,iLep)) { 
					TLorentzVector lep4v (LepInfo.Px[iLep],LepInfo.Py[iLep],LepInfo.Pz[iLep],LepInfo.Energy[iLep]) ;
					TLorentzVector jet4v (JetInfo.Px[iJet],JetInfo.Py[iJet],JetInfo.Pz[iJet],JetInfo.Energy[iJet]) ; 
					if (lep4v.DeltaR(jet4v)<.5) { jetLepOverlap = true; break ; }
				} 
			}
			if (jetLepOverlap) continue ;  // Reject jets overlapping with leptons 

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

			selectedjet.index[selectedjet.Size] = iJet ; 
			selectedjet.flavour[selectedjet.Size] = TMath::Abs(JetInfo.GenFlavor[iJet]) ; 
			JetInfo.CombinedSVBJetTags[iJet] < .244 ? selectedjet.btagUntag[selectedjet.Size]  = 1 : selectedjet.btagUntag[selectedjet.Size]  = 0 ; 
			JetInfo.CombinedSVBJetTags[iJet] > .244 ? selectedjet.btagLoose[selectedjet.Size]  = 1 : selectedjet.btagLoose[selectedjet.Size]  = 0 ; 
			JetInfo.CombinedSVBJetTags[iJet] > .679 ? selectedjet.btagMedium[selectedjet.Size] = 1 : selectedjet.btagMedium[selectedjet.Size] = 0 ; 
			JetInfo.CombinedSVBJetTags[iJet] > .898 ? selectedjet.btagTight[selectedjet.Size]  = 1 : selectedjet.btagTight[selectedjet.Size]  = 0 ; 
			selectedjet.mass  [selectedjet.Size] = JetInfo.Mass[iJet]   ; 
			selectedjet.energy[selectedjet.Size] = JetInfo.Energy[iJet] ; 
			selectedjet.pt    [selectedjet.Size] = JetInfo.Pt[iJet]     ; 
			selectedjet.eta   [selectedjet.Size] = JetInfo.Eta[iJet]    ; 
			selectedjet.phi   [selectedjet.Size] = JetInfo.Phi[iJet]    ; 
			selectedjet.Size++ ; 
			if (selectedjet.btagMedium[selectedjet.Size] == 1) ++nbtags ;

		} // selectedjet 
		m_hnjets_->Fill(selectedjet.Size,weightPU0) ; 
		m_hnbtags_->Fill(nbtags,weightPU0) ;
		m_hcutFlow_->AddBinContent(2,1) ; 

		/**\ Select jets with JES +1sigma */

		jet_jesHigh.Size = 0 ; 
		nbtags = 0 ;  

		for (int iJet=0;iJet < JetInfo.Size;++iJet) {

			bool jetLepOverlap(false) ; 

			double jeshigh = (1. + JetInfo.Unc[iJet]) ; 
			TLorentzVector jet4vhigh(JetInfo.Px[iJet]*jeshigh,JetInfo.Py[iJet]*jeshigh,JetInfo.Pz[iJet]*jeshigh,JetInfo.Energy[iJet]*jeshigh) ;

			for (int iLep=0; iLep < LepInfo.Size; ++iLep) {
				if(isGoodMuon(lepInfo,iEntry,iLep) || isGoodElectron(evtinfo,lepInfo,iEntry,iLep)) { 
					TLorentzVector lep4v (LepInfo.Px[iLep],LepInfo.Py[iLep],LepInfo.Pz[iLep],LepInfo.Energy[iLep]) ;
					if (lep4v.DeltaR(jet4vhigh) < .5) { jetLepOverlap = true; break ; } 
				} 
			}
			if (jetLepOverlap) continue ;  // Reject jets overlapping with leptons 

			if (jet4vhigh.Pt() < 30. ||
					TMath::Abs(jet4vhigh.PseudoRapidity()) > 2.4 || 
					JetInfo.CHF[iJet] <= 0.00 ||
					JetInfo.NHF[iJet] >= 0.99 ||
					JetInfo.CEF[iJet] >= 0.99 ||
					JetInfo.NEF[iJet] >= 0.99 ||
					JetInfo.NCH[iJet] <= 0.00  
			   ) {
				continue; 
			}  // apply jet ID 

			jet_jesHigh.index[jet_jesHigh.Size] = iJet ; 
			jet_jesHigh.flavour[jet_jesHigh.Size] = TMath::Abs(JetInfo.GenFlavor[iJet]) ; 
			JetInfo.CombinedSVBJetTags[iJet] < .244 ? jet_jesHigh.btagUntag[jet_jesHigh.Size]  = 1 : jet_jesHigh.btagUntag[jet_jesHigh.Size]  = 0 ; 
			JetInfo.CombinedSVBJetTags[iJet] > .244 ? jet_jesHigh.btagLoose[jet_jesHigh.Size]  = 1 : jet_jesHigh.btagLoose[jet_jesHigh.Size]  = 0 ; 
			JetInfo.CombinedSVBJetTags[iJet] > .679 ? jet_jesHigh.btagMedium[jet_jesHigh.Size] = 1 : jet_jesHigh.btagMedium[jet_jesHigh.Size] = 0 ; 
			JetInfo.CombinedSVBJetTags[iJet] > .898 ? jet_jesHigh.btagTight[jet_jesHigh.Size]  = 1 : jet_jesHigh.btagTight[jet_jesHigh.Size]  = 0 ; 
			jet_jesHigh.mass  [jet_jesHigh.Size] = jet4vhigh.Mag()    ; 
			jet_jesHigh.energy[jet_jesHigh.Size] = jet4vhigh.Energy() ; 
			jet_jesHigh.pt    [jet_jesHigh.Size] = jet4vhigh.Pt()     ; 
			jet_jesHigh.eta   [jet_jesHigh.Size] = jet4vhigh.PseudoRapidity()    ; 
			jet_jesHigh.phi   [jet_jesHigh.Size] = jet4vhigh.Phi()    ; 
			jet_jesHigh.Size++ ; 
			if (jet_jesHigh.btagMedium[jet_jesHigh.Size] == 1) ++nbtags ;

		} // jet_jesHigh 

		/**\ Select jets with JES -1sigma */

		jet_jesLow.Size = 0 ; 
		nbtags = 0 ;  

		for (int iJet=0;iJet < JetInfo.Size;++iJet) {

			bool jetLepOverlap(false) ; 

			double jeslow = (1. - JetInfo.Unc[iJet]) ; 
			TLorentzVector jet4vlow(JetInfo.Px[iJet]*jeslow,JetInfo.Py[iJet]*jeslow,JetInfo.Pz[iJet]*jeslow,JetInfo.Energy[iJet]*jeslow) ;

			for (int iLep=0; iLep < LepInfo.Size; ++iLep) {
				if(isGoodMuon(lepInfo,iEntry,iLep) || isGoodElectron(evtinfo,lepInfo,iEntry,iLep)) { 
					TLorentzVector lep4v (LepInfo.Px[iLep],LepInfo.Py[iLep],LepInfo.Pz[iLep],LepInfo.Energy[iLep]) ;
					if (lep4v.DeltaR(jet4vlow) < .5) { jetLepOverlap = true; break ; } 
				} 
			}
			if (jetLepOverlap) continue ;  // Reject jets overlapping with leptons 

			if (jet4vlow.Pt() < 30. ||
					TMath::Abs(jet4vlow.PseudoRapidity()) > 2.4 || 
					JetInfo.CHF[iJet] <= 0.00 ||
					JetInfo.NHF[iJet] >= 0.99 ||
					JetInfo.CEF[iJet] >= 0.99 ||
					JetInfo.NEF[iJet] >= 0.99 ||
					JetInfo.NCH[iJet] <= 0.00  
			   ) {
				continue; 
			}  // apply jet ID 

			jet_jesLow.index[jet_jesLow.Size] = iJet ; 
			jet_jesLow.flavour[jet_jesLow.Size] = TMath::Abs(JetInfo.GenFlavor[iJet]) ; 
			JetInfo.CombinedSVBJetTags[iJet] < .244 ? jet_jesLow.btagUntag[jet_jesLow.Size]  = 1 : jet_jesLow.btagUntag[jet_jesLow.Size]  = 0 ; 
			JetInfo.CombinedSVBJetTags[iJet] > .244 ? jet_jesLow.btagLoose[jet_jesLow.Size]  = 1 : jet_jesLow.btagLoose[jet_jesLow.Size]  = 0 ; 
			JetInfo.CombinedSVBJetTags[iJet] > .679 ? jet_jesLow.btagMedium[jet_jesLow.Size] = 1 : jet_jesLow.btagMedium[jet_jesLow.Size] = 0 ; 
			JetInfo.CombinedSVBJetTags[iJet] > .898 ? jet_jesLow.btagTight[jet_jesLow.Size]  = 1 : jet_jesLow.btagTight[jet_jesLow.Size]  = 0 ; 
			jet_jesLow.mass  [jet_jesLow.Size] = jet4vlow.Mag()    ; 
			jet_jesLow.energy[jet_jesLow.Size] = jet4vlow.Energy() ; 
			jet_jesLow.pt    [jet_jesLow.Size] = jet4vlow.Pt()     ; 
			jet_jesLow.eta   [jet_jesLow.Size] = jet4vlow.PseudoRapidity()       ; 
			jet_jesLow.phi   [jet_jesLow.Size] = jet4vlow.Phi()       ; 
			jet_jesLow.Size++ ; 
			if (jet_jesLow.btagMedium[jet_jesLow.Size] == 1) ++nbtags ;

		} // jet_jesLow 

		if (!isData_) { 
			/**\ Select jets with JER +1sigma */

			jet_jerHigh.Size = 0 ; 
			nbtags = 0 ; 

			for (int iJet=0;iJet < JetInfo.Size;++iJet) {

				bool jetLepOverlap(false) ; 

				double deltapt = 0.1*TMath::Abs(JetInfo.Pt[iJet] - JetInfo.GenJetPt[iJet]) ;
				double jerhigh = max(0.,double(1. + deltapt/JetInfo.Pt[iJet])) ;
				TLorentzVector jer4vhigh(JetInfo.Px[iJet]*jerhigh,JetInfo.Py[iJet]*jerhigh,JetInfo.Pz[iJet]*jerhigh,JetInfo.Energy[iJet]*jerhigh) ; 

				for (int iLep=0; iLep < LepInfo.Size; ++iLep) {
					if(isGoodMuon(lepInfo,iEntry,iLep) || isGoodElectron(evtinfo,lepInfo,iEntry,iLep)) { 
						TLorentzVector lep4v (LepInfo.Px[iLep],LepInfo.Py[iLep],LepInfo.Pz[iLep],LepInfo.Energy[iLep]) ;
						if (lep4v.DeltaR(jer4vhigh) < .5) { jetLepOverlap = true; break ; } 
					} 
				}
				if (jetLepOverlap) continue ;  // Reject jets overlapping with leptons 

				if ( jer4vhigh.Pt() < 0.001 ) std::cout << " >>>>>>>>>>>>>>> JERHigh PT = " << jer4vhigh.Pt() << std::endl ; 

				if (jer4vhigh.Pt() < 30. ||
						TMath::Abs(jer4vhigh.PseudoRapidity()) > 2.4 || 
						JetInfo.CHF[iJet] <= 0.00 ||
						JetInfo.NHF[iJet] >= 0.99 ||
						JetInfo.CEF[iJet] >= 0.99 ||
						JetInfo.NEF[iJet] >= 0.99 ||
						JetInfo.NCH[iJet] <= 0.00  
				   ) {
					continue; 
				}  // apply jet ID 

				jet_jerHigh.index[jet_jerHigh.Size] = iJet ; 
				jet_jerHigh.flavour[jet_jerHigh.Size] = TMath::Abs(JetInfo.GenFlavor[iJet]) ; 
				JetInfo.CombinedSVBJetTags[iJet] < .244 ? jet_jerHigh.btagUntag[jet_jerHigh.Size]  = 1 : jet_jerHigh.btagUntag[jet_jerHigh.Size]  = 0 ; 
				JetInfo.CombinedSVBJetTags[iJet] > .244 ? jet_jerHigh.btagLoose[jet_jerHigh.Size]  = 1 : jet_jerHigh.btagLoose[jet_jerHigh.Size]  = 0 ; 
				JetInfo.CombinedSVBJetTags[iJet] > .679 ? jet_jerHigh.btagMedium[jet_jerHigh.Size] = 1 : jet_jerHigh.btagMedium[jet_jerHigh.Size] = 0 ; 
				JetInfo.CombinedSVBJetTags[iJet] > .898 ? jet_jerHigh.btagTight[jet_jerHigh.Size]  = 1 : jet_jerHigh.btagTight[jet_jerHigh.Size]  = 0 ; 
				jet_jerHigh.mass  [jet_jerHigh.Size] = jer4vhigh.Mag()    ; 
				jet_jerHigh.energy[jet_jerHigh.Size] = jer4vhigh.Energy() ; 
				jet_jerHigh.pt    [jet_jerHigh.Size] = jer4vhigh.Pt()     ; 
				jet_jerHigh.eta   [jet_jerHigh.Size] = jer4vhigh.Eta()    ; 
				jet_jerHigh.phi   [jet_jerHigh.Size] = jer4vhigh.Phi()    ; 
				jet_jerHigh.Size++ ; 
				if (jet_jerHigh.btagMedium[jet_jerHigh.Size] == 1) ++nbtags ;

			} // jet_jerHigh 

			/**\ Select jets with JER -1sigma */

			jet_jerLow.Size = 0 ; 
			nbtags = 0 ; 

			for (int iJet=0;iJet < JetInfo.Size;++iJet) {

				bool jetLepOverlap(false) ; 

				double deltapt = 0.1*TMath::Abs(JetInfo.Pt[iJet] - JetInfo.GenJetPt[iJet]) ;
				if (deltapt/JetInfo.Pt[iJet] > 1.) deltapt = 0.1*JetInfo.Pt[iJet] ; 
				double jerlow  = max(0.,double(1. - deltapt/JetInfo.Pt[iJet])) ; 
				TLorentzVector jer4vlow(JetInfo.Px[iJet]*jerlow,JetInfo.Py[iJet]*jerlow,JetInfo.Pz[iJet]*jerlow,JetInfo.Energy[iJet]*jerlow) ;

				for (int iLep=0; iLep < LepInfo.Size; ++iLep) {
					if(isGoodMuon(lepInfo,iEntry,iLep) || isGoodElectron(evtinfo,lepInfo,iEntry,iLep)) { 
						TLorentzVector lep4v (LepInfo.Px[iLep],LepInfo.Py[iLep],LepInfo.Pz[iLep],LepInfo.Energy[iLep]) ;
						if (lep4v.DeltaR(jer4vlow) < .5) { jetLepOverlap = true; break ; } 
					} 
				}
				if (jetLepOverlap) continue ;  // Reject jets overlapping with leptons 

				if (jer4vlow.Pt() < 30. ||
						TMath::Abs(jer4vlow.PseudoRapidity()) > 2.4 || 
						JetInfo.CHF[iJet] <= 0.00 ||
						JetInfo.NHF[iJet] >= 0.99 ||
						JetInfo.CEF[iJet] >= 0.99 ||
						JetInfo.NEF[iJet] >= 0.99 ||
						JetInfo.NCH[iJet] <= 0.00  
				   ) {
					continue; 
				}  // apply jet ID 

				jet_jerLow.index[jet_jerLow.Size] = iJet ; 
				jet_jerLow.flavour[jet_jerLow.Size] = TMath::Abs(JetInfo.GenFlavor[iJet]) ; 
				JetInfo.CombinedSVBJetTags[iJet] < .244 ? jet_jerLow.btagUntag[jet_jerLow.Size]  = 1 : jet_jerLow.btagUntag[jet_jerLow.Size]  = 0 ; 
				JetInfo.CombinedSVBJetTags[iJet] > .244 ? jet_jerLow.btagLoose[jet_jerLow.Size]  = 1 : jet_jerLow.btagLoose[jet_jerLow.Size]  = 0 ; 
				JetInfo.CombinedSVBJetTags[iJet] > .679 ? jet_jerLow.btagMedium[jet_jerLow.Size] = 1 : jet_jerLow.btagMedium[jet_jerLow.Size] = 0 ; 
				JetInfo.CombinedSVBJetTags[iJet] > .898 ? jet_jerLow.btagTight[jet_jerLow.Size]  = 1 : jet_jerLow.btagTight[jet_jerLow.Size]  = 0 ; 
				jet_jerLow.mass  [jet_jerLow.Size] = jer4vlow.Mag()    ; 
				jet_jerLow.energy[jet_jerLow.Size] = jer4vlow.Energy() ; 
				jet_jerLow.pt    [jet_jerLow.Size] = jer4vlow.Pt()     ; 
				jet_jerLow.eta   [jet_jerLow.Size] = jer4vlow.Eta()    ; 
				jet_jerLow.phi   [jet_jerLow.Size] = jer4vlow.Phi()    ; 
				jet_jerLow.Size++ ; 
				if (jet_jerLow.btagMedium[jet_jerLow.Size] == 1) ++nbtags ;

			} // jet_jerLow 
		} 

		os2lTree->Fill() ; 

	} // End of event loop 

	std::cout << " \nOSSF_el = " << nossf_el << "\n" ; 
	std::cout << " \nOSSF_mu = " << nossf_mu << "\n" ; 
	std::cout << " \nOSOF = " << nosof << "\n" ; 

	fntuple->Write() ; 
	fntuple->Close() ; 

} // OS2LSelector::evtLoop () 

#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h" 
#include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h" 

bool OS2LSelector::isGoodElectron (EvtInfoBranches* EvtInfo, LepInfoBranches* LepInfo, int& entry, int& ele) { 

	if ( TMath::Abs(LepInfo->LeptonType[ele]) != 11 ) return false ;
	if ( LepInfo->Pt[ele] < 30. ) return false ;

	if (LepInfo->EgammaCutBasedEleIdMEDIUM[ele] == true) { return true ; }
	else { return false ; }

	int idx = -1;
	if((fabs(LepInfo->Eta[ele])<=1.4442)){
		idx = 0;
	} else if(((fabs(LepInfo->Eta[ele])<2.4)&&(fabs(LepInfo->Eta[ele])>=1.566))){ 
		idx = 1;
	}

	if( idx == -1 ) return false ;

	float cut_dEtaIn[2]         = {999.9, 999.9};
	float cut_dPhiIn[2]         = {999.9, 999.9};
	float cut_sigmaIEtaIEta[2]  = {999.9, 999.9};
	float cut_hoe[2]            = {999.9, 999.9};
	float cut_ooemoop[2]        = {999.9, 999.9};
	float cut_d0vtx[2]          = {999.9, 999.9};
	float cut_dzvtx[2]          = {999.9, 999.9};
	float cut_iso[2]            = {999.9, 999.9};
	bool cut_vtxFit[2]          = {false, false};
	unsigned int cut_mHits[2]   = {999, 999};

	cut_dEtaIn[0]        = 0.004; cut_dEtaIn[1]        = 0.007;
	cut_dPhiIn[0]        = 0.060; cut_dPhiIn[1]        = 0.030;
	cut_sigmaIEtaIEta[0] = 0.010; cut_sigmaIEtaIEta[1] = 0.030;
	cut_hoe[0]           = 0.120; cut_hoe[1]           = 0.100;
	cut_ooemoop[0]       = 0.050; cut_ooemoop[1]       = 0.050;
	cut_d0vtx[0]         = 0.020; cut_d0vtx[1]         = 0.020;
	cut_dzvtx[0]         = 0.100; cut_dzvtx[1]         = 0.100;
	cut_iso[0]           = 0.150; cut_iso[1]           = 0.150;
	cut_vtxFit[0]        = true ; cut_vtxFit[1]        = true;
	cut_mHits[0]         = 1    ; cut_mHits[1]         = 1;

	const float        pt                = LepInfo->Pt[ele] ;
	const float        dEtaIn            = LepInfo->EldeltaEta[ele] ; 
	const float        dPhiIn            = LepInfo->EldeltaPhi[ele] ; 
	const float        sigmaIEtaIEta     = LepInfo->ElsigmaIetaIeta[ele] ; 
	const float        hoe               = LepInfo->ElHadoverEm[ele] ; 
	const float        ooemoop           = ((1.0/LepInfo->ElEcalE[ele]) - (LepInfo->ElEoverP[ele]/LepInfo->ElEcalE[ele])) ; 
	const float        d0vtx             = LepInfo->ElTrackDxy_PV[ele] ; 
	const float        dzvtx             = LepInfo->ElTrackDz[ele] ; 
	const float        iso_ch            = LepInfo->ChargedHadronIso[ele] ; 
	const float        iso_em            = LepInfo->PhotonIso[ele] ;
	const float        iso_nh            = LepInfo->NeutralHadronIso[ele] ; 
	const bool         vtxFitConversion  = LepInfo->ElhasConv[ele] ; 
	const unsigned int mHits             = LepInfo->NumberOfExpectedInnerHits[ele] ; 
	const double       rho               = std::max(double(EvtInfo->RhoPU[0]),0.0); ; 

	// effective area for isolation
	float AEff = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, LepInfo->Eta[ele], ElectronEffectiveArea::kEleEAData2012);

	// apply to neutrals 
	double iso_n = std::max((iso_nh + iso_em) - (rho * AEff), 0.0); 

	// compute final isolation 
	double iso = (iso_n + iso_ch) / pt;

	bool isGoodElectronCutbased(false) ; 

	// test cuts
	if (fabs(dEtaIn) < cut_dEtaIn[idx])             
		if (fabs(dPhiIn) < cut_dPhiIn[idx])             
			if (sigmaIEtaIEta < cut_sigmaIEtaIEta[idx])     
				if (hoe < cut_hoe[idx])                         
					if (fabs(ooemoop) < cut_ooemoop[idx])           
						if (fabs(d0vtx) < cut_d0vtx[idx])               
							if (fabs(dzvtx) < cut_dzvtx[idx])               
								if (!cut_vtxFit[idx] || !vtxFitConversion)      
									if (mHits <= cut_mHits[idx])                    
										if (iso < cut_iso[idx])                         
											if ( (LepInfo->Charge[ele] == 1 && LepInfo->ChargeCtf[ele] == 1 && LepInfo->ChargeGsf[ele] == 1) || (LepInfo->Charge[ele] == -1 && LepInfo->ChargeCtf[ele] == -1 && LepInfo->ChargeGsf[ele] == -1) )  
												isGoodElectronCutbased = true ;  

	return isGoodElectronCutbased ; 

}

bool OS2LSelector:: isGoodMuon (LepInfoBranches* LepInfo,int& entry,int& mu) { 

	bool isGoodMuon(false) ; 

	if(   TMath::Abs(LepInfo->LeptonType[mu]) == 13 
			&& LepInfo->isPFMuon[mu] == true 
			&& LepInfo->innerTracknormalizedChi2[mu] < 10 
			&& LepInfo->MuNMuonhits[mu] > 0
			&& LepInfo->MuNMatchedStations[mu] > 1
			&& TMath::Abs(LepInfo->MuInnerTrackDxy_PV[mu]) < 0.2
			&& TMath::Abs(LepInfo->MuInnerTrackDz[mu]) < 0.5
			&& LepInfo->MuNPixelLayers[mu] > 0
			&& LepInfo->MuNTrackLayersWMeasurement[mu] > 5  
			&& LepInfo->Pt[mu] > 30.
			&& TMath::Abs(LepInfo->Eta[mu]) < 2.1 
			&& (( LepInfo->ChargedHadronIsoR04[mu] + max(LepInfo->NeutralHadronIsoR04[mu] + LepInfo->PhotonIsoR04[mu] - 0.5*LepInfo->sumPUPtR04[mu], 0.0))/LepInfo->Pt[mu]<0.15) 
	  ) {
		isGoodMuon = true ; 	   
	} else isGoodMuon = false ; 

	return isGoodMuon ; 	

}

bool OS2LSelector:: isGoodTau (LepInfoBranches* LepInfo,int& entry,int& tau) {
	bool isGoodTau(false) ;
	if ( LepInfo->isPFTau[tau] 
			//&& LepInfo->againstElectronLoose[tau] > .5 
			&& LepInfo->againstElectronMVA[tau] > .5 
			&& LepInfo->againstMuonLoose[tau] > .5 
			&& LepInfo->Pt[tau] > 20.
			&& LepInfo->byLooseCombinedIsolationDeltaBetaCorr[tau] == 1) isGoodTau = true ; 
	else isGoodTau = false ; 
	return isGoodTau ; 
}

