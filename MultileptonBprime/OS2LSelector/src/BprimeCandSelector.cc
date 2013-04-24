#include "MultileptonBprime/OS2LSelector/interface/BprimeCandSelector.h"
#include "BprimeKit/Utils/interface/mathutils.h"
#include "BprimeKit/Headers/interface/format.h" 
#include "BprimeKit/Utils/interface/BTagWeight.h"
#include "MultileptonBprime/OS2LSelector/interface/BookHistograms.h" 

#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <sstream> 
#include <assert.h>
#include <TLorentzVector.h> 

using namespace std; 

double BprimeCandSelector::STBins[NSTBins] = {0.,250.,500.,750.,1000.,1500.,15000.} ;    

/**\ Constructor */
BprimeCandSelector::BprimeCandSelector (const std::string&infiles,const bool&isData,const std::string&jsonfile,const std::string&outfile) :
	OS2LSelector(infiles,isData,jsonfile,outfile), 	
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
	std::cout << " Creating this instance of BprimeCandSelector\n " ; 

	/**\ Initialize histograms */
	outfile_.cd() ; 
	hist_ = new BookHistograms() ;

	TDirectory* dyee = outfile_.mkdir("DYee") ; 
	dyee->cd() ; 
	TDirectory* dyee_DY = dyee->mkdir("DY") ; 
	TDirectory* dyee_DY4Jets = dyee->mkdir("DY4Jets") ; 
	TDirectory* dyee_DYSTLt500 = dyee->mkdir("DYSTLt500") ; 
	TDirectory* dyee_DYSTGt500 = dyee->mkdir("DYSTGt500") ; 
	TDirectory* dyee_DYSTGt500NoBTag = dyee->mkdir("DYSTGt500NoBTag") ; 
	dyee_DY->cd() ; 
	hist_DYee_DY_ = new BookHistograms() ; 
	dyee_DY4Jets->cd() ; 
	hist_DYee_DY4Jets_ = new BookHistograms() ; 
	dyee_DYSTLt500->cd() ; 
	hist_DYee_DYSTLt500_ = new BookHistograms() ;
	dyee_DYSTGt500->cd() ; 
	hist_DYee_DYSTGt500_ = new BookHistograms() ;
	dyee_DYSTGt500NoBTag->cd() ; 
	hist_DYee_DYSTGt500NoBTag_ = new BookHistograms() ; 
	outfile_.cd() ; 

	TDirectory* dymumu = outfile_.mkdir("DYmumu") ; 
	dymumu->cd() ; 
	TDirectory* dymumu_DY = dymumu->mkdir("DY") ; 
	TDirectory* dymumu_DY4Jets = dymumu->mkdir("DY4Jets") ; 
	TDirectory* dymumu_DYSTLt500 = dymumu->mkdir("DYSTLt500") ; 
	TDirectory* dymumu_DYSTGt500 = dymumu->mkdir("DYSTGt500") ; 
	TDirectory* dymumu_DYSTGt500NoBTag = dymumu->mkdir("DYSTGt500NoBTag") ; 
	dymumu_DY->cd() ; 
	hist_DYmumu_DY_ = new BookHistograms() ; 
	dymumu_DY4Jets->cd() ; 
	hist_DYmumu_DY4Jets_ = new BookHistograms() ; 
	dymumu_DYSTLt500->cd() ; 
	hist_DYmumu_DYSTLt500_ = new BookHistograms() ;
	dymumu_DYSTGt500->cd() ; 
	hist_DYmumu_DYSTGt500_ = new BookHistograms() ;
	dymumu_DYSTGt500NoBTag->cd() ; 
	hist_DYmumu_DYSTGt500NoBTag_ = new BookHistograms() ; 
	outfile_.cd() ; 

	TDirectory* emu = outfile_.mkdir("EMu") ; 
	emu->cd() ; 
	TDirectory* emu_DY = emu->mkdir("DY") ; 
	TDirectory* emu_DY4Jets = emu->mkdir("DY4Jets") ; 
	TDirectory* emu_DYSTLt500 = emu->mkdir("DYSTLt500") ; 
	TDirectory* emu_DYSTGt500 = emu->mkdir("DYSTGt500") ; 
	TDirectory* emu_DYSTGt500NoBTag = emu->mkdir("DYSTGt500NoBTag") ; 
	emu_DY->cd() ; 
	hist_EMu_DY_ = new BookHistograms() ; 
	emu_DY4Jets->cd() ; 
	hist_EMu_DY4Jets_ = new BookHistograms() ; 
	emu_DYSTLt500->cd() ; 
	hist_EMu_DYSTLt500_ = new BookHistograms() ;
	emu_DYSTGt500->cd() ; 
	hist_EMu_DYSTGt500_ = new BookHistograms() ;
	emu_DYSTGt500NoBTag->cd() ; 
	hist_EMu_DYSTGt500NoBTag_ = new BookHistograms() ; 
	outfile_.cd() ; 

	m_hst[0][0] = new TH1D("m_hst_ee_isZVetoLow", "ST ee isZVetoLow", NSTBins-1,STBins) ; 
	m_hst[0][1] = new TH1D("m_hst_ee_isZVetoHigh","ST ee isZVetoHigh",NSTBins-1,STBins) ; 
	m_hst[0][2] = new TH1D("m_hst_ee_isZWin",     "ST ee isZWin",     NSTBins-1,STBins) ; 
	m_hst[1][0] = new TH1D("m_hst_mm_isZVetoLow", "ST mm isZVetoLow", NSTBins-1,STBins) ; 
	m_hst[1][1] = new TH1D("m_hst_mm_isZVetoHigh","ST mm isZVetoHigh",NSTBins-1,STBins) ; 
	m_hst[1][2] = new TH1D("m_hst_mm_isZWin",     "ST mm isZWin",     NSTBins-1,STBins) ; 
	m_hst[2][0] = new TH1D("m_hst_em_isZVetoLow", "ST em isZVetoLow", NSTBins-1,STBins) ; 
	m_hst[2][1] = new TH1D("m_hst_em_isZVetoHigh","ST em isZVetoHigh",NSTBins-1,STBins) ; 
	m_hst[2][2] = new TH1D("m_hst_em_isZWin",     "ST em isZWin",     NSTBins-1,STBins) ; 

	m_hevtsel_      = new TH1D("hevtsel","Selected events",21,.5,21.5) ; 
	m_hevtsel_ee_   = new TH1D("hevtsel_ee","Selected events ee ch.",21,.5,21.5) ; 
	m_hevtsel_emu_  = new TH1D("hevtsel_emu","Selected events emu ch.",21,.5,21.5) ; 
	m_hevtsel_mumu_ = new TH1D("hevtsel_mumu","Selected events mumu ch.",21,.5,21.5) ; 

	if (!isData_ && (inputfiles_.find("BprimeBprimeTo") != std::string::npos) ) { 

		int nbr(0) ; 
		int sumbr(0) ; 
		for (int itw = 0; itw <= 10; ++ itw ) {
			for (int ibz = 0; ibz <= 10; ++ibz) {
				for (int ibh = 0; ibh <= 10; ++ibh) {  
					sumbr = 0 ; 
					sumbr += itw ; 
					sumbr += ibz ; 
					sumbr += ibh ; 
					if ( sumbr == 10 ) {  
						std::cout << " BR(tW) : BR(bZ) : BR(bH) :: " << itw << " : " << ibz << " : " << ibh <<  " sum br = " << sumbr/10. << std::endl ;  

						std::stringstream hname_ee_isZVetoLow  ; 
						std::stringstream hname_ee_isZVetoHigh ; 
						std::stringstream hname_ee_isZWin      ; 
						std::stringstream hname_mm_isZVetoLow  ; 
						std::stringstream hname_mm_isZVetoHigh ; 
						std::stringstream hname_mm_isZWin      ; 
						std::stringstream hname_em_isZVetoLow  ; 
						std::stringstream hname_em_isZVetoHigh ; 
						std::stringstream hname_em_isZWin      ; 
						hname_ee_isZVetoLow  << "hst_ee_isZVetoLow"  << "_" << itw << "_" << ibz << "_" << ibh ; 
						hname_ee_isZVetoHigh << "hst_ee_isZVetoHigh" << "_" << itw << "_" << ibz << "_" << ibh ; 
						hname_ee_isZWin      << "hst_ee_isZWin"      << "_" << itw << "_" << ibz << "_" << ibh ; 
						hname_mm_isZVetoLow  << "hst_mm_isZVetoLow"  << "_" << itw << "_" << ibz << "_" << ibh ; 
						hname_mm_isZVetoHigh << "hst_mm_isZVetoHigh" << "_" << itw << "_" << ibz << "_" << ibh ; 
						hname_mm_isZWin      << "hst_mm_isZWin"      << "_" << itw << "_" << ibz << "_" << ibh ; 
						hname_em_isZVetoLow  << "hst_em_isZVetoLow"  << "_" << itw << "_" << ibz << "_" << ibh ; 
						hname_em_isZVetoHigh << "hst_em_isZVetoHigh" << "_" << itw << "_" << ibz << "_" << ibh ; 
						hname_em_isZWin      << "hst_em_isZWin"      << "_" << itw << "_" << ibz << "_" << ibh ; 

						m_hst_sig[0][0][nbr] = new TH1D((hname_ee_isZVetoLow .str()).c_str(),"",NSTBins-1,STBins) ; 
						m_hst_sig[0][1][nbr] = new TH1D((hname_ee_isZVetoHigh.str()).c_str(),"",NSTBins-1,STBins) ; 
						m_hst_sig[0][2][nbr] = new TH1D((hname_ee_isZWin     .str()).c_str(),"",NSTBins-1,STBins) ; 
						m_hst_sig[1][0][nbr] = new TH1D((hname_mm_isZVetoLow .str()).c_str(),"",NSTBins-1,STBins) ; 
						m_hst_sig[1][1][nbr] = new TH1D((hname_mm_isZVetoHigh.str()).c_str(),"",NSTBins-1,STBins) ; 
						m_hst_sig[1][2][nbr] = new TH1D((hname_mm_isZWin     .str()).c_str(),"",NSTBins-1,STBins) ; 
						m_hst_sig[2][0][nbr] = new TH1D((hname_em_isZVetoLow .str()).c_str(),"",NSTBins-1,STBins) ; 
						m_hst_sig[2][1][nbr] = new TH1D((hname_em_isZVetoHigh.str()).c_str(),"",NSTBins-1,STBins) ; 
						m_hst_sig[2][2][nbr] = new TH1D((hname_em_isZWin     .str()).c_str(),"",NSTBins-1,STBins) ; 

						std::stringstream hname_ee ; 
						std::stringstream hname_mm ; 
						std::stringstream hname_em ; 
						hname_ee << "hevtselSig_ee_" << itw << "_"  << ibz << "_" << ibh ; 
						hname_mm << "hevtselSig_mm_" << itw << "_"  << ibz << "_" << ibh ; 
						hname_em << "hevtselSig_em_" << itw << "_"  << ibz << "_" << ibh ; 
						m_hevtselSig_ee_   [nbr] = new TH1D((hname_ee.str()).c_str(),((hname_ee).str()).c_str(),51,.5,51.5) ; 

						m_hevtselSig_emu_  [nbr] = new TH1D((hname_mm.str()).c_str(),((hname_mm).str()).c_str(),51,.5,51.5) ; 
						m_hevtselSig_mumu_ [nbr] = new TH1D((hname_em.str()).c_str(),((hname_em).str()).c_str(),51,.5,51.5) ; 
						++nbr ; 
					}
				}
			}
		}
		std::cout << " nbr = " << nbr << std::endl ; 
	}

}

/**\ Destructor */ 
BprimeCandSelector::~BprimeCandSelector () {
}

void BprimeCandSelector::process () {

	char inputFile[2000] ; 
	ifstream infile ; 
	infile.open ((char*)inputfiles_.c_str(),ifstream::in) ; 
	assert(!infile.fail()) ;
	thisChain_ = new TChain("os2lTree") ; 

	while (infile.good()) { 
		infile >> inputFile ; 
		if(strncmp(inputFile,"#",1)==0) continue; 
		if(strncmp(inputFile,"%",1)==0) break ; 
		thisChain_->Add(inputFile); 
		std::cout << "\ninfile = " << inputFile << std::endl ; 
	}

	infile.close() ;   
	evtProc_ = thisChain_->GetEntries(); 
	evtLoop();
	thisChain_->Reset() ;

	outfile_.Write() ; 
	outfile_.Close() ; 

	return ; 

}

#include "BprimeKit/Headers/interface/SelectedJet.h" 
#include "BprimeKit/Headers/interface/JetCollection.h" 

void::BprimeCandSelector::evtLoop () {

	GenInfoBranches    GenInfo ; 
	EvtInfoBranches    EvtInfo ; 
	VertexInfoBranches VtxInfo ; 
	LepInfoBranches    LepInfo ; 
	JetInfoBranches    JetInfo ; 
	SelectedJet    selectedjet ;
	SelectedJet    jet_jesHigh ; 
	SelectedJet    jet_jesLow  ; 
	SelectedJet    jet_jerHigh ; 
	SelectedJet    jet_jerLow  ; 

	double weightPU0(1)   ; // Event flag 
	double weightPUP(1)   ; // Event flag 
	double weightPUM(1)   ; // Event flag 
	int    nGoodVtxs(0)   ; // Event flag  
	int    elPIndex       ; // Event flag  
	int    elMIndex       ; // Event flag  
	int    muPIndex       ; // Event flag  
	int    muMIndex       ; // Event flag  
	bool   isZVeto(0)     ; // Event flag 
	bool   isZVetoLow(0)  ; // Event flag 
	bool   isZVetoHigh(0) ; // Event flag 
	bool   isOSSF_el(0)   ; // Event flag  
	bool   isOSSF_mu(0)   ; // Event flag  
	bool   isOSOF(0)      ; // Event flag  
	double ST(0.) ; 
	double ST_JESHigh(0.) ; 
	double ST_JESLow(0.) ; 
	double ST_JERHigh(0.) ; 
	double ST_JERLow(0.) ; 

	GenInfo.Register(thisChain_) ; 
	EvtInfo.Register(thisChain_) ; 
	VtxInfo.Register(thisChain_) ; 
	LepInfo.Register(thisChain_,"PFLepInfo") ; 
	JetInfo.Register(thisChain_,"PFJetInfo") ; 
	selectedjet.Register(thisChain_) ; 
	jet_jesHigh.Register(thisChain_,"Jets_JESHigh") ; 
	jet_jesLow.Register(thisChain_,"Jets_JESLow") ; 
	jet_jerHigh.Register(thisChain_,"Jets_JERHigh") ; 
	jet_jerLow.Register(thisChain_,"Jets_JERLow") ; 

	thisChain_->SetBranchAddress("weightPU0",   &weightPU0)    ; 
	thisChain_->SetBranchAddress("weightPUP",   &weightPUP)    ; 
	thisChain_->SetBranchAddress("weightPUM",   &weightPUM)    ; 
	thisChain_->SetBranchAddress("nGoodVtxs",   &nGoodVtxs)    ; 
	thisChain_->SetBranchAddress("elPIndex",    &elPIndex)     ;
	thisChain_->SetBranchAddress("elMIndex",    &elMIndex)     ; 
	thisChain_->SetBranchAddress("muPIndex",    &muPIndex)     ;
	thisChain_->SetBranchAddress("muMIndex",    &muMIndex)     ; 
	thisChain_->SetBranchAddress("isOSSF_el",   &isOSSF_el)    ; 
	thisChain_->SetBranchAddress("isOSSF_mu",   &isOSSF_mu)    ; 
	thisChain_->SetBranchAddress("isOSOF",      &isOSOF)       ; 
	thisChain_->SetBranchAddress("isZVeto",     &isZVeto)      ; 
	thisChain_->SetBranchAddress("isZVetoLow",  &isZVetoLow)   ; 
	thisChain_->SetBranchAddress("isZVetoHigh", &isZVetoHigh)  ; 

	TString name = "NTuple_" ;
	name += outfileName_ ; 
	TFile* fntuple = new TFile(name,"RECREATE") ; 
	fntuple->cd() ; 
	TTree* bprimeTree = thisChain_->CloneTree(0); 
	bprimeTree->Branch("ST", &ST, "ST/D") ; 
	bprimeTree->Branch("ST_JESHigh",&ST_JESHigh,"ST_JESHigh/D") ; 
	bprimeTree->Branch("ST_JESLow", &ST_JESLow, "ST_JESLow/D") ; 
	bprimeTree->Branch("ST_JERHigh",&ST_JERHigh,"ST_JERHigh/D") ; 
	bprimeTree->Branch("ST_JERLow", &ST_JERLow, "ST_JERLow/D") ; 
	bprimeTree->SetName("bprimeTree") ; 

	std::cout << " Events processed = " << thisChain_->GetEntries() << std::endl ; 

	/**\ Looping over events */ 
	for (int iEntry = 0; iEntry < thisChain_->GetEntries(); ++iEntry) { 

		m_hevtsel_->Fill(1) ; 

		ST = 0. ; 
		ST_JESHigh = 0. ; 
		ST_JESLow = 0. ; 
		ST_JERHigh = 0. ; 
		ST_JERLow = 0. ; 

		isZVeto = true ; 
		isZVetoLow = true ; 
		isZVetoHigh = true ; 

		bool btagm(false) ; 
		bool btagm_jesHigh(false) ; 
		bool btagm_jesLow(false) ; 
		bool btagm_jerHigh(false) ; 
		bool btagm_jerLow(false) ; 

		int nbjets(0) ; 
		int nbjets_jesHigh(0) ; 
		int nbjets_jesLow(0) ; 
		int nbjets_jerHigh(0) ; 
		int nbjets_jerLow(0) ; 

		std::vector<JetVariables> JetCollection; //// Needed for b-tag weight 

		thisChain_->GetEntry(iEntry) ; 

		EvtInfoBranches* evtinfo = &EvtInfo ; 
		LepInfoBranches* lepInfo = &LepInfo ; 

		if (!isOSSF_el && !isOSSF_mu && !isOSOF) continue ; 

		//// Identify leading OS leptons 

		int lep0(-1) ;
		int lep1(-1) ;
		//if (isOSSF_el) {
		//	lep0 = elPIndex ; lep1 = elMIndex ; 
		//}
		//else if (isOSSF_mu)  {
		//	lep0 = muPIndex ; lep1 = muMIndex ; 
		//}
		//else if (isOSOF)  {
		//	if (elPIndex >= 0 && muMIndex >= 0) { lep0 = elPIndex ; lep1 = muMIndex ; } 
		//	else if (elMIndex >= 0 && muPIndex >= 0) { lep0 = elMIndex ; lep1 = muPIndex ; } 
		//}

		for (int iLep=0; iLep < LepInfo.Size; ++iLep) {
			if (isGoodElectron(evtinfo,lepInfo,iEntry,iLep) || isGoodMuon(lepInfo,iEntry,iLep)) {
				lep0 = iLep ; 
				break ;
			}
		}

		for (int jLep=lep0+1; jLep < LepInfo.Size; ++jLep) { 
			if (isGoodElectron(evtinfo,lepInfo,iEntry,jLep) || isGoodMuon(lepInfo,iEntry,jLep)) {
				if (LepInfo.Charge[lep0]*LepInfo.Charge[jLep] < 0) { 
					lep1 = jLep ;  
					break ; 
				}
			}
		}

		TLorentzVector lep0_4v(LepInfo.Px[lep0],LepInfo.Py[lep0],LepInfo.Pz[lep0],LepInfo.Energy[lep0]) ; 
		TLorentzVector lep1_4v(LepInfo.Px[lep1],LepInfo.Py[lep1],LepInfo.Pz[lep1],LepInfo.Energy[lep1]) ; 
		if (isOSSF_el) {
			hist_DYee_DY_->fill(hist_DYee_DY_->m_hmee_,(lep0_4v+lep1_4v).Mag(),weightPU0) ; 
			if ((lep0_4v+lep1_4v).Mag() > 76. && (lep0_4v+lep1_4v).Mag() < 106.) isZVeto = false ; 
			else {
				isZVeto = true ; 
				if ((lep0_4v+lep1_4v).Mag() < 76.) { isZVetoHigh = false ; isZVetoLow = true ; } 
				else if ((lep0_4v+lep1_4v).Mag() > 106.) { isZVetoHigh = true ; isZVetoLow = false ; } 
			}
		}
		if (isOSSF_mu) {
			hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hmmumu_,(lep0_4v+lep1_4v).Mag(),weightPU0) ; 
			if ((lep0_4v+lep1_4v).Mag() > 76. && (lep0_4v+lep1_4v).Mag() < 106.) isZVeto = false ; 
			else {
				isZVeto = true ; 
				if ((lep0_4v+lep1_4v).Mag() < 76.) { isZVetoHigh = false ; isZVetoLow = true ; } 
				else if ((lep0_4v+lep1_4v).Mag() > 106.) { isZVetoHigh = true ; isZVetoLow = false ; } 
			}
		}
		if (isOSOF) {
			hist_EMu_DY_->fill(hist_EMu_DY_->m_hmemu_,(lep0_4v+lep1_4v).Mag(),weightPU0) ; 
		}

		int el0(-1) ;
		int el1(-1) ;
		//if (LepInfo.LeptonType[lep0] == 13) el0 = lep0 ; 
		//if (LepInfo.LeptonType[lep1] == 13) el1 = lep1 ; 
		for (int iLep=0; iLep < LepInfo.Size; ++iLep) {
			if (isGoodElectron(evtinfo,lepInfo,iEntry,iLep)) {
				el0 = iLep ; break ; 
			}
		}
		for (int jLep=el0+1; jLep < LepInfo.Size; ++jLep) { 
			if (isGoodElectron(evtinfo,lepInfo,iEntry,jLep)) {
				if (LepInfo.Charge[el0]*LepInfo.Charge[jLep] < 0) { 
					el1 = jLep ;  
					break ; 
				}
			}
		}

		if (el0 > -1) {
			if (isOSSF_el) {
				hist_DYee_DY_->fill(hist_DYee_DY_->m_hpte0_,LepInfo.Pt[el0],weightPU0) ; 
				hist_DYee_DY_->fill(hist_DYee_DY_->m_hetae0_,LepInfo.Eta[el0],weightPU0) ; 
				if (isZVeto) {
					hist_DYee_DY_->fill(hist_DYee_DY_->m_hpte0_notdy_,LepInfo.Pt[el0],weightPU0) ; 
					hist_DYee_DY_->fill(hist_DYee_DY_->m_hetae0_notdy_,LepInfo.Eta[el0],weightPU0) ; 
				}
				else {
					hist_DYee_DY_->fill(hist_DYee_DY_->m_hpte0_dy_,LepInfo.Pt[el0],weightPU0) ; 
					hist_DYee_DY_->fill(hist_DYee_DY_->m_hetae0_dy_,LepInfo.Eta[el0],weightPU0) ; 
				}
			}
			if (isOSOF) {
				hist_EMu_DY_->fill(hist_EMu_DY_->m_hpte0_,LepInfo.Pt[el0],weightPU0) ; 
				hist_EMu_DY_->fill(hist_EMu_DY_->m_hetae0_,LepInfo.Eta[el0],weightPU0) ; 
				if (isZVeto) {
					hist_EMu_DY_->fill(hist_EMu_DY_->m_hpte0_notdy_,LepInfo.Pt[el0],weightPU0) ; 
					hist_EMu_DY_->fill(hist_EMu_DY_->m_hetae0_notdy_,LepInfo.Eta[el0],weightPU0) ; 
				}
				else {
					hist_EMu_DY_->fill(hist_EMu_DY_->m_hpte0_dy_,LepInfo.Pt[el0],weightPU0) ; 
					hist_EMu_DY_->fill(hist_EMu_DY_->m_hetae0_dy_,LepInfo.Eta[el0],weightPU0) ; 
				}
			}
		}
		if (el1 > -1) {
			if (isOSSF_el) {
				hist_DYee_DY_->fill(hist_DYee_DY_->m_hpte1_,LepInfo.Pt[el0],weightPU0) ; 
				hist_DYee_DY_->fill(hist_DYee_DY_->m_hetae1_,LepInfo.Eta[el0],weightPU0) ; 
				if (isZVeto) {
					hist_DYee_DY_->fill(hist_DYee_DY_->m_hpte1_notdy_,LepInfo.Pt[el0],weightPU0) ; 
					hist_DYee_DY_->fill(hist_DYee_DY_->m_hetae1_notdy_,LepInfo.Eta[el0],weightPU0) ; 
				}
				else {
					hist_DYee_DY_->fill(hist_DYee_DY_->m_hpte1_dy_,LepInfo.Pt[el0],weightPU0) ; 
					hist_DYee_DY_->fill(hist_DYee_DY_->m_hetae1_dy_,LepInfo.Eta[el0],weightPU0) ; 
				}
			}
			if (isOSOF) {
				hist_EMu_DY_->fill(hist_EMu_DY_->m_hpte1_,LepInfo.Pt[el0],weightPU0) ; 
				hist_EMu_DY_->fill(hist_EMu_DY_->m_hetae1_,LepInfo.Eta[el0],weightPU0) ; 
				if (isZVeto) {
					hist_EMu_DY_->fill(hist_EMu_DY_->m_hpte1_notdy_,LepInfo.Pt[el0],weightPU0) ; 
					hist_EMu_DY_->fill(hist_EMu_DY_->m_hetae1_notdy_,LepInfo.Eta[el0],weightPU0) ; 
				}
				else {
					hist_EMu_DY_->fill(hist_EMu_DY_->m_hpte1_dy_,LepInfo.Pt[el0],weightPU0) ; 
					hist_EMu_DY_->fill(hist_EMu_DY_->m_hetae1_dy_,LepInfo.Eta[el0],weightPU0) ; 
				}
			}
		}

		int mu0(-1) ;
		int mu1(-1) ;
		//if (LepInfo.LeptonType[lep0] == 13) mu0 = lep0 ; 
		//if (LepInfo.LeptonType[lep1] == 13) mu1 = lep1 ; 
		for (int iLep=0; iLep < LepInfo.Size; ++iLep) {
			if (isGoodMuon(lepInfo,iEntry,iLep)) {
				mu0 = iLep ; break ; 
			}
		}
		for (int jLep=mu0+1; jLep < LepInfo.Size; ++jLep) { 
			if (isGoodMuon(lepInfo,iEntry,jLep) && (LepInfo.Charge[mu0]*LepInfo.Charge[jLep] < 0)) { 
				mu1 = jLep ;  break ; 
			}
		}

		if (mu0 > -1) {
			if (isOSSF_mu) {
				hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hptm0_,LepInfo.Pt[mu0],weightPU0) ; 
				hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hetam0_,LepInfo.Eta[mu0],weightPU0) ; 
				if (isZVeto) {
					hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hptm0_notdy_,LepInfo.Pt[mu0],weightPU0) ; 
					hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hetam0_notdy_,LepInfo.Eta[mu0],weightPU0) ; 
				}
				else {
					hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hptm0_dy_,LepInfo.Pt[mu0],weightPU0) ; 
					hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hetam0_dy_,LepInfo.Eta[mu0],weightPU0) ; 
				}
			}
			if (isOSOF) {
				hist_EMu_DY_->fill(hist_EMu_DY_->m_hptm0_,LepInfo.Pt[mu0],weightPU0) ; 
				hist_EMu_DY_->fill(hist_EMu_DY_->m_hetam0_,LepInfo.Eta[mu0],weightPU0) ; 
				if (isZVeto) {
					hist_EMu_DY_->fill(hist_EMu_DY_->m_hptm0_notdy_,LepInfo.Pt[mu0],weightPU0) ; 
					hist_EMu_DY_->fill(hist_EMu_DY_->m_hetam0_notdy_,LepInfo.Eta[mu0],weightPU0) ; 
				}
				else {
					hist_EMu_DY_->fill(hist_EMu_DY_->m_hptm0_dy_,LepInfo.Pt[mu0],weightPU0) ; 
					hist_EMu_DY_->fill(hist_EMu_DY_->m_hetam0_dy_,LepInfo.Eta[mu0],weightPU0) ; 
				}
			}
		}
		if (mu1 > -1) {
		}
		if (isOSSF_mu) {
			hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hptm1_,LepInfo.Pt[mu1],weightPU0) ; 
			hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hetam1_,LepInfo.Eta[mu1],weightPU0) ; 
			if (isZVeto) {
				hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hptm1_notdy_,LepInfo.Pt[mu1],weightPU0) ; 
				hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hetam1_notdy_,LepInfo.Eta[mu1],weightPU0) ; 
			}
			else {
				hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hptm1_dy_,LepInfo.Pt[mu1],weightPU0) ; 
				hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hetam1_dy_,LepInfo.Eta[mu1],weightPU0) ; 
			}
		}
		if (isOSOF) {
			hist_EMu_DY_->fill(hist_EMu_DY_->m_hptm1_,LepInfo.Pt[mu1],weightPU0) ; 
			hist_EMu_DY_->fill(hist_EMu_DY_->m_hetam1_,LepInfo.Eta[mu1],weightPU0) ; 
			if (isZVeto) {
				hist_EMu_DY_->fill(hist_EMu_DY_->m_hptm1_notdy_,LepInfo.Pt[mu1],weightPU0) ; 
				hist_EMu_DY_->fill(hist_EMu_DY_->m_hetam1_notdy_,LepInfo.Eta[mu1],weightPU0) ; 
			}
			else {
				hist_EMu_DY_->fill(hist_EMu_DY_->m_hptm1_dy_,LepInfo.Pt[mu1],weightPU0) ; 
				hist_EMu_DY_->fill(hist_EMu_DY_->m_hetam1_dy_,LepInfo.Eta[mu1],weightPU0) ; 
			}
		}

		//// Calculate ST 
		ST         = LepInfo.Pt[lep0] + LepInfo.Pt[lep1] ; 
		ST_JESHigh = LepInfo.Pt[lep0] + LepInfo.Pt[lep1] ; 
		ST_JESLow  = LepInfo.Pt[lep0] + LepInfo.Pt[lep1] ; 
		ST_JERHigh = LepInfo.Pt[lep0] + LepInfo.Pt[lep1] ; 
		ST_JERLow  = LepInfo.Pt[lep0] + LepInfo.Pt[lep1] ; 

		ST         += EvtInfo.PFMET ; 
		ST_JESHigh += EvtInfo.PFMET ; 
		ST_JESLow  += EvtInfo.PFMET ; 
		ST_JERHigh += EvtInfo.PFMET ; 
		ST_JERLow  += EvtInfo.PFMET ; 

		for (int iJet = 0; iJet < selectedjet.Size; ++iJet) { 
			ST += selectedjet.pt[iJet] ; 
			if ( selectedjet.btagMedium[iJet] == true ) { btagm = true ; ++nbjets ; } 
			int flavour = selectedjet.flavour[iJet] ; 
			double ptj = selectedjet.pt[iJet] ;
			double etaj = selectedjet.eta[iJet] ;
			bool btag = selectedjet.btagMedium[iJet] ; 
			JetVariables Jet(iJet,btag,flavour,ptj,etaj) ;
			JetCollection.push_back(Jet) ;
		}

		for (int iJet = 0; iJet < jet_jesHigh.Size; ++iJet) { 
			ST_JESHigh += jet_jesHigh.pt[iJet] ; 
			if ( jet_jesHigh.btagMedium[iJet] == true ) { btagm_jesHigh = true ;  ++nbjets_jesHigh ; }
		}

		for (int iJet = 0; iJet < jet_jesLow.Size; ++iJet) { 
			ST_JERLow += jet_jesLow.pt[iJet] ; 
			if ( jet_jesLow.btagMedium[iJet] == true ) { btagm_jesLow = true ; ++nbjets_jesLow ; } 
		}

		for (int iJet = 0; iJet < jet_jerHigh.Size; ++iJet) { 
			ST_JERHigh += jet_jerHigh.pt[iJet] ; 
			if ( jet_jerHigh.btagMedium[iJet] == true ) { btagm_jerHigh = true ; ++nbjets_jerHigh ; } 
		}

		for (int iJet = 0; iJet < jet_jerLow.Size; ++iJet) { 
			ST_JESLow += jet_jerLow.pt[iJet] ; 
			if ( jet_jerLow.btagMedium[iJet] == true ) { btagm_jerLow = true ; ++nbjets_jerLow ; }  
		}

		if (isOSSF_el) {
			hist_DYee_DY_->fill(hist_DYee_DY_->m_hnjets_,selectedjet.Size,weightPU0) ; 
			hist_DYee_DY_->fill(hist_DYee_DY_->m_hst_,ST,weightPU0) ; 
			hist_DYee_DY_->fill(hist_DYee_DY_->m_hmet_,EvtInfo.PFMET,weightPU0) ; 
			if (isZVeto) {
				hist_DYee_DY_->fill(hist_DYee_DY_->m_hnjets_notdy_,selectedjet.Size,weightPU0) ; 
				hist_DYee_DY_->fill(hist_DYee_DY_->m_hst_notdy_,ST,weightPU0) ; 
				hist_DYee_DY_->fill(hist_DYee_DY_->m_hmet_notdy_,EvtInfo.PFMET,weightPU0) ; 
			}
			else {
				hist_DYee_DY_->fill(hist_DYee_DY_->m_hnjets_dy_,selectedjet.Size,weightPU0) ; 
				hist_DYee_DY_->fill(hist_DYee_DY_->m_hst_dy_,ST,weightPU0) ; 
				hist_DYee_DY_->fill(hist_DYee_DY_->m_hmet_dy_,EvtInfo.PFMET,weightPU0) ; 
			}
		}

		if (isOSSF_mu) {
			hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hst_,ST,weightPU0) ; 
			hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hmet_,EvtInfo.PFMET,weightPU0) ; 
			hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hnjets_,selectedjet.Size,weightPU0) ; 
			if (isZVeto) {
				hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hnjets_notdy_,selectedjet.Size,weightPU0) ; 
				hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hst_notdy_,ST,weightPU0) ; 
				hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hmet_notdy_,EvtInfo.PFMET,weightPU0) ; 
			}
			else {
				hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hnjets_dy_,selectedjet.Size,weightPU0) ; 
				hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hst_dy_,ST,weightPU0) ; 
				hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hmet_dy_,EvtInfo.PFMET,weightPU0) ; 
			}
		}

		if (isOSOF) {
			hist_EMu_DY_->fill(hist_EMu_DY_->m_hnjets_,selectedjet.Size,weightPU0) ; 
			hist_EMu_DY_->fill(hist_EMu_DY_->m_hst_,ST,weightPU0) ; 
			hist_EMu_DY_->fill(hist_EMu_DY_->m_hmet_,EvtInfo.PFMET,weightPU0) ; 
			if (isZVeto) {
				hist_EMu_DY_->fill(hist_EMu_DY_->m_hnjets_notdy_,selectedjet.Size,weightPU0) ; 
				hist_EMu_DY_->fill(hist_EMu_DY_->m_hst_notdy_,ST,weightPU0) ; 
				hist_EMu_DY_->fill(hist_EMu_DY_->m_hmet_notdy_,EvtInfo.PFMET,weightPU0) ; 
			}
			else {
				hist_EMu_DY_->fill(hist_EMu_DY_->m_hnjets_dy_,selectedjet.Size,weightPU0) ; 
				hist_EMu_DY_->fill(hist_EMu_DY_->m_hst_dy_,ST,weightPU0) ; 
				hist_EMu_DY_->fill(hist_EMu_DY_->m_hmet_dy_,EvtInfo.PFMET,weightPU0) ; 
			}
		}

		//// Get b-tagging weight

		double btagWtMean = 1. ;
		double btagWtHigh_bc = 1. ;
		double btagWtLow_bc = 1. ;
		double btagWtHigh_l = 1. ;
		double btagWtLow_l = 1. ;

		BTagWeight bTagWeight(JetCollection,isData_,.679,0,1,nbjets) ;
		TString tag="mean" ;
		btagWtMean = bTagWeight.weight(tag) ;
		tag="errorP_bc" ;
		btagWtHigh_bc = bTagWeight.weight(tag) ;
		tag="errorM_bc" ;
		btagWtLow_bc = bTagWeight.weight(tag) ;
		tag="errorP_l" ;
		btagWtHigh_l = bTagWeight.weight(tag) ;
		tag="errorM_l" ;
		btagWtLow_l = bTagWeight.weight(tag) ;

		double evtwt              = 1. ; 
		double evtwt_btag_bcHigh  = 1. ; 
		double evtwt_btag_bcLow   = 1. ; 
		double evtwt_btag_lHigh   = 1. ; 
		double evtwt_btag_lLow    = 1. ; 
		if (!isData_) {
			evtwt              *= weightPU0*btagWtMean ; 
			evtwt_btag_bcHigh  *= weightPU0*btagWtHigh_bc ; 
			evtwt_btag_bcLow   *= weightPU0*btagWtLow_bc ; 
			evtwt_btag_lHigh   *= weightPU0*btagWtHigh_l ; 
			evtwt_btag_lLow    *= weightPU0*btagWtLow_l ; 
		}

		if (isOSSF_el) {
			hist_DYee_DY_->fill(hist_DYee_DY_->m_hnbjets_,nbjets,evtwt) ; 
			if (isZVeto) {
				hist_DYee_DY_->fill(hist_DYee_DY_->m_hnbjets_notdy_,nbjets,evtwt) ; 
			}
			else {
				hist_DYee_DY_->fill(hist_DYee_DY_->m_hnbjets_dy_,nbjets,evtwt) ; 
			}
		}

		if (isOSSF_mu) {
			hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hnbjets_,nbjets,evtwt) ; 
			if (isZVeto) {
				hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hnbjets_notdy_,nbjets,evtwt) ; 
			}
			else {
				hist_DYmumu_DY_->fill(hist_DYmumu_DY_->m_hnbjets_dy_,nbjets,evtwt) ; 
			}
		}

		if (isOSOF) {
			hist_EMu_DY_->fill(hist_EMu_DY_->m_hnbjets_,nbjets,evtwt) ; 
			if (isZVeto) {
				hist_EMu_DY_->fill(hist_EMu_DY_->m_hnbjets_notdy_,nbjets,evtwt) ; 
			}
			else {
				hist_EMu_DY_->fill(hist_EMu_DY_->m_hnbjets_dy_,nbjets,evtwt) ; 
			}
		}

		if (selectedjet.Size >= 4) {

			if (btagm == true) {
				if (isOSSF_el && isZVetoLow)  m_hst[0][0] -> Fill(ST,evtwt) ; 
				if (isOSSF_el && isZVetoHigh) m_hst[0][1] -> Fill(ST,evtwt) ; 
				if (isOSSF_el && !isZVeto)    m_hst[0][2] -> Fill(ST,evtwt) ; 
				if (isOSSF_mu && isZVetoLow)  m_hst[1][0] -> Fill(ST,evtwt) ; 
				if (isOSSF_mu && isZVetoHigh) m_hst[1][1] -> Fill(ST,evtwt) ; 
				if (isOSSF_mu && !isZVeto)    m_hst[1][2] -> Fill(ST,evtwt) ; 
				if (isOSOF && isZVetoLow)     m_hst[2][0] -> Fill(ST,evtwt) ; 
				if (isOSOF && isZVetoHigh)    m_hst[2][1] -> Fill(ST,evtwt) ; 
				if (isOSOF && !isZVeto)       m_hst[2][2] -> Fill(ST,evtwt) ; 
			}

			if (isOSSF_el) {
				if (!isZVeto) hist_DYee_DY4Jets_->fill(hist_DYee_DY4Jets_->m_hst_dy_,ST,weightPU0) ; 
				else hist_DYee_DY4Jets_->fill(hist_DYee_DY4Jets_->m_hst_notdy_,ST,weightPU0) ; 
			}
			if (isOSSF_mu) {
				if (!isZVeto) hist_DYmumu_DY4Jets_->fill(hist_DYmumu_DY4Jets_->m_hst_dy_,ST,weightPU0) ; 
				else hist_DYmumu_DY4Jets_->fill(hist_DYmumu_DY4Jets_->m_hst_notdy_,ST,weightPU0) ; 
			}
			if (isOSOF) {
				if (!isZVeto) hist_EMu_DY4Jets_->fill(hist_EMu_DY4Jets_->m_hst_dy_,ST,weightPU0) ; 
				else hist_EMu_DY4Jets_->fill(hist_EMu_DY4Jets_->m_hst_notdy_,ST,weightPU0) ; 
			}

			if (ST > 500.) {
				if (btagm == true) {
					m_hevtsel_->Fill(2,evtwt) ; 
					if (isOSSF_el) {
						hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hmee_,(lep0_4v+lep1_4v).Mag(),weightPU0) ;
						hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hpte0_,LepInfo.Pt[el0],weightPU0) ; 
						hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hetae0_,LepInfo.Eta[el0],weightPU0) ; 
						if (isZVeto) {
							hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hpte0_notdy_,LepInfo.Pt[el0],weightPU0) ; 
							hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hetae0_notdy_,LepInfo.Eta[el0],weightPU0) ; 
						}
						else {
							hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hpte0_dy_,LepInfo.Pt[el0],weightPU0) ; 
							hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hetae0_dy_,LepInfo.Eta[el0],weightPU0) ; 
						}
						hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hpte1_,LepInfo.Pt[el0],weightPU0) ; 
						hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hetae1_,LepInfo.Eta[el0],weightPU0) ; 
						if (isZVeto) {
							hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hpte1_notdy_,LepInfo.Pt[el0],weightPU0) ; 
							hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hetae1_notdy_,LepInfo.Eta[el0],weightPU0) ; 
						}
						else {
							hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hpte1_dy_,LepInfo.Pt[el0],weightPU0) ; 
							hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hetae1_dy_,LepInfo.Eta[el0],weightPU0) ; 
						}
						hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hst_,ST,weightPU0) ; 
						hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hmet_,EvtInfo.PFMET,weightPU0) ; 
						if (isZVeto) {
							hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hst_notdy_,ST,weightPU0) ; 
							hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hmet_notdy_,EvtInfo.PFMET,weightPU0) ; 
						}
						else {
							hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hst_dy_,ST,weightPU0) ; 
							hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hmet_dy_,EvtInfo.PFMET,weightPU0) ; 
						}

						hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hnjets_,selectedjet.Size,weightPU0) ; 
						if (isZVeto) {
							hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hnjets_notdy_,selectedjet.Size,weightPU0) ; 
						}
						else {
							hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hnjets_dy_,selectedjet.Size,weightPU0) ; 
						}
						hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hnbjets_,nbjets,evtwt) ; 
						if (isZVeto) {
							hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hnbjets_notdy_,nbjets,evtwt) ; 
						}
						else {
							hist_DYee_DYSTGt500_->fill(hist_DYee_DYSTGt500_->m_hnbjets_dy_,nbjets,evtwt) ; 
						}
					} // isOSSF_el 
					if (isOSSF_mu) {
						hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hmmumu_,(lep0_4v+lep1_4v).Mag(),weightPU0)  ; 
						hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hptm0_,LepInfo.Pt[mu0],weightPU0) ; 
						hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hetam0_,LepInfo.Eta[mu0],weightPU0) ; 
						if (isZVeto) {
							hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hptm0_notdy_,LepInfo.Pt[mu0],weightPU0) ; 
							hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hetam0_notdy_,LepInfo.Eta[mu0],weightPU0) ; 
						}
						else {
							hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hptm0_dy_,LepInfo.Pt[mu0],weightPU0) ; 
							hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hetam0_dy_,LepInfo.Eta[mu0],weightPU0) ; 
						}
						hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hptm1_,LepInfo.Pt[mu1],weightPU0) ; 
						hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hetam1_,LepInfo.Eta[mu1],weightPU0) ; 
						if (isZVeto) {
							hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hptm1_notdy_,LepInfo.Pt[mu1],weightPU0) ; 
							hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hetam1_notdy_,LepInfo.Eta[mu1],weightPU0) ; 
						}
						else {
							hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hptm1_dy_,LepInfo.Pt[mu1],weightPU0) ; 
							hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hetam1_dy_,LepInfo.Eta[mu1],weightPU0) ; 
						}
						hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hst_,ST,weightPU0) ; 
						hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hmet_,EvtInfo.PFMET,weightPU0) ; 
						if (isZVeto) {
							hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hst_notdy_,ST,weightPU0) ; 
							hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hmet_notdy_,EvtInfo.PFMET,weightPU0) ; 
						}
						else {
							hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hst_dy_,ST,weightPU0) ; 
							hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hmet_dy_,EvtInfo.PFMET,weightPU0) ; 
						}

						hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hnjets_,selectedjet.Size,weightPU0) ; 
						if (isZVeto) {
							hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hnjets_notdy_,selectedjet.Size,weightPU0) ; 
						}
						else {
							hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hnjets_dy_,selectedjet.Size,weightPU0) ; 
						}
						hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hnbjets_,nbjets,evtwt) ; 
						if (isZVeto) {
							hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hnbjets_notdy_,nbjets,evtwt) ; 
						}
						else {
							hist_DYmumu_DYSTGt500_->fill(hist_DYmumu_DYSTGt500_->m_hnbjets_dy_,nbjets,evtwt) ; 
						}
					} // isOSSF_mu 
					if (isOSOF) {
						hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hmemu_,(lep0_4v+lep1_4v).Mag(),weightPU0) ; 
						hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hpte0_,LepInfo.Pt[el0],weightPU0) ; 
						hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hetae0_,LepInfo.Eta[el0],weightPU0) ; 
						if (isZVeto) {
							hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hpte0_notdy_,LepInfo.Pt[el0],weightPU0) ; 
							hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hetae0_notdy_,LepInfo.Eta[el0],weightPU0) ; 
						}
						else {
							hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hpte0_dy_,LepInfo.Pt[el0],weightPU0) ; 
							hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hetae0_dy_,LepInfo.Eta[el0],weightPU0) ; 
						}
						hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hpte1_,LepInfo.Pt[el0],weightPU0) ; 
						hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hetae1_,LepInfo.Eta[el0],weightPU0) ; 
						if (isZVeto) {
							hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hpte1_notdy_,LepInfo.Pt[el0],weightPU0) ; 
							hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hetae1_notdy_,LepInfo.Eta[el0],weightPU0) ; 
						}
						else {
							hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hpte1_dy_,LepInfo.Pt[el0],weightPU0) ; 
							hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hetae1_dy_,LepInfo.Eta[el0],weightPU0) ; 
						}
						hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hptm0_,LepInfo.Pt[mu0],weightPU0) ; 
						hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hetam0_,LepInfo.Eta[mu0],weightPU0) ; 
						if (isZVeto) {
							hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hptm0_notdy_,LepInfo.Pt[mu0],weightPU0) ; 
							hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hetam0_notdy_,LepInfo.Eta[mu0],weightPU0) ; 
						}
						else {
							hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hptm0_dy_,LepInfo.Pt[mu0],weightPU0) ; 
							hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hetam0_dy_,LepInfo.Eta[mu0],weightPU0) ; 
						}
						hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hptm1_,LepInfo.Pt[mu1],weightPU0) ; 
						hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hetam1_,LepInfo.Eta[mu1],weightPU0) ; 
						if (isZVeto) {
							hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hptm1_notdy_,LepInfo.Pt[mu1],weightPU0) ; 
							hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hetam1_notdy_,LepInfo.Eta[mu1],weightPU0) ; 
						}
						else {
							hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hptm1_dy_,LepInfo.Pt[mu1],weightPU0) ; 
							hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hetam1_dy_,LepInfo.Eta[mu1],weightPU0) ; 
						}
						hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hst_,ST,weightPU0) ; 
						hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hmet_,EvtInfo.PFMET,weightPU0) ; 
						if (isZVeto) {
							hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hst_notdy_,ST,weightPU0) ; 
							hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hmet_notdy_,EvtInfo.PFMET,weightPU0) ; 
						}
						else {
							hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hst_dy_,ST,weightPU0) ; 
							hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hmet_dy_,EvtInfo.PFMET,weightPU0) ; 
						}

						hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hnjets_,selectedjet.Size,weightPU0) ; 
						if (isZVeto) {
							hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hnjets_notdy_,selectedjet.Size,weightPU0) ; 
						}
						else {
							hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hnjets_dy_,selectedjet.Size,weightPU0) ; 
						}
						hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hnbjets_,nbjets,evtwt) ; 
						if (isZVeto) {
							hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hnbjets_notdy_,nbjets,evtwt) ; 
						}
						else {
							hist_EMu_DYSTGt500_->fill(hist_EMu_DYSTGt500_->m_hnbjets_dy_,nbjets,evtwt) ; 
						}
					} // isOSOF 
				} // has b tagged jet selection 
			}  // ST > 500 GeV selection 
		} // Njets > = 4 selection 

		if (isOSSF_el || isOSSF_mu || isOSOF) {
			if (jet_jesHigh.Size >= 4)  
				if (ST_JESHigh > 500.) 
					if (btagm_jesHigh == true) 
						m_hevtsel_->Fill(3,evtwt) ; 

			if (jet_jesLow.Size >= 4)  
				if (ST_JESLow > 500.) 
					if (btagm_jesLow == true) 
						m_hevtsel_->Fill(4,evtwt) ; 

			if (jet_jerHigh.Size >= 4)  
				if (ST_JESHigh > 500.) 
					if (btagm_jerHigh == true) 
						m_hevtsel_->Fill(5,evtwt) ; 

			if (jet_jerLow.Size >= 4)  
				if (ST_JESLow > 500.) 
					if (btagm_jerLow == true) 
						m_hevtsel_->Fill(6,evtwt) ; 

			if (selectedjet.Size >= 4)  
				if (ST_JESLow > 500.) 
					if (btagm == true) 
						m_hevtsel_->Fill(7,evtwt_btag_bcHigh) ; 

			if (selectedjet.Size >= 4)  
				if (ST_JESLow > 500.) 
					if (btagm == true) 
						m_hevtsel_->Fill(8,evtwt_btag_bcLow) ; 

			if (selectedjet.Size >= 4)  
				if (ST_JESLow > 500.) 
					if (btagm == true) 
						m_hevtsel_->Fill(9,evtwt_btag_lHigh) ; 

			if (selectedjet.Size >= 4)  
				if (ST_JESLow > 500.) 
					if (btagm == true) 
						m_hevtsel_->Fill(10,evtwt_btag_lLow) ; 
		}

		//// Filling dielectron channel 
		if (isOSSF_el && !isOSSF_mu && !isOSOF) {
			if (selectedjet.Size >= 4)  
				if (ST > 500.) 
					if (btagm == true)			  
						m_hevtsel_ee_->Fill(2,evtwt) ; 

			if (jet_jesHigh.Size >= 4)  
				if (ST_JESHigh > 500.) 
					if (btagm_jesHigh == true) 
						m_hevtsel_ee_->Fill(3,evtwt) ; 

			if (jet_jesLow.Size >= 4)  
				if (ST_JESLow > 500.) 
					if (btagm_jesLow == true) 
						m_hevtsel_ee_->Fill(4,evtwt) ; 

			if (jet_jerHigh.Size >= 4)  
				if (ST_JESHigh > 500.) 
					if (btagm_jerHigh == true) 
						m_hevtsel_ee_->Fill(5,evtwt) ; 

			if (jet_jerLow.Size >= 4)  
				if (ST_JESLow > 500.) 
					if (btagm_jerLow == true) 
						m_hevtsel_ee_->Fill(6,evtwt) ; 

			if (selectedjet.Size >= 4)  
				if (ST_JESLow > 500.) 
					if (btagm == true) 
						m_hevtsel_ee_->Fill(7,evtwt_btag_bcHigh) ; 

			if (selectedjet.Size >= 4)  
				if (ST_JESLow > 500.) 
					if (btagm == true) 
						m_hevtsel_ee_->Fill(8,evtwt_btag_bcLow) ; 

			if (selectedjet.Size >= 4)  
				if (ST_JESLow > 500.) 
					if (btagm == true) 
						m_hevtsel_ee_->Fill(9,evtwt_btag_lHigh) ; 

			if (selectedjet.Size >= 4)  
				if (ST_JESLow > 500.) 
					if (btagm == true) 
						m_hevtsel_ee_->Fill(10,evtwt_btag_lLow) ; 

		}

		//// Filling dimuon channel 
		if (isOSSF_mu && !isOSSF_el && !isOSOF) {
			if (selectedjet.Size >= 4)  
				if (ST > 500.) 
					if (btagm == true)			  
						m_hevtsel_mumu_->Fill(2,evtwt) ; 

			if (jet_jesHigh.Size >= 4)  
				if (ST_JESHigh > 500.) 
					if (btagm_jesHigh == true) 
						m_hevtsel_mumu_->Fill(3,evtwt) ; 

			if (jet_jesLow.Size >= 4)  
				if (ST_JESLow > 500.) 
					if (btagm_jesLow == true) 
						m_hevtsel_mumu_->Fill(4,evtwt) ; 

			if (jet_jerHigh.Size >= 4)  
				if (ST_JESHigh > 500.) 
					if (btagm_jerHigh == true) 
						m_hevtsel_mumu_->Fill(5,evtwt) ; 

			if (jet_jerLow.Size >= 4)  
				if (ST_JESLow > 500.) 
					if (btagm_jerLow == true) 
						m_hevtsel_mumu_->Fill(6,evtwt) ; 

			if (selectedjet.Size >= 4)  
				if (ST_JESLow > 500.) 
					if (btagm == true) 
						m_hevtsel_mumu_->Fill(7,evtwt_btag_bcHigh) ; 

			if (selectedjet.Size >= 4)  
				if (ST_JESLow > 500.) 
					if (btagm == true) 
						m_hevtsel_mumu_->Fill(8,evtwt_btag_bcLow) ; 

			if (selectedjet.Size >= 4)  
				if (ST_JESLow > 500.) 
					if (btagm == true) 
						m_hevtsel_mumu_->Fill(9,evtwt_btag_lHigh) ; 

			if (selectedjet.Size >= 4)  
				if (ST_JESLow > 500.) 
					if (btagm == true) 
						m_hevtsel_mumu_->Fill(10,evtwt_btag_lLow) ; 

		}

		//// Filling emu channel 
		if (isOSOF && !isOSSF_el && !isOSSF_mu) {
			if (selectedjet.Size >= 4)  
				if (ST > 500.) 
					if (btagm == true)			  
						m_hevtsel_emu_->Fill(2,evtwt) ; 

			if (jet_jesHigh.Size >= 4)  
				if (ST_JESHigh > 500.) 
					if (btagm_jesHigh == true) 
						m_hevtsel_emu_->Fill(3,evtwt) ; 

			if (jet_jesLow.Size >= 4)  
				if (ST_JESLow > 500.) 
					if (btagm_jesLow == true) 
						m_hevtsel_emu_->Fill(4,evtwt) ; 

			if (jet_jerHigh.Size >= 4)  
				if (ST_JESHigh > 500.) 
					if (btagm_jerHigh == true) 
						m_hevtsel_emu_->Fill(5,evtwt) ; 

			if (jet_jerLow.Size >= 4)  
				if (ST_JESLow > 500.) 
					if (btagm_jerLow == true) 
						m_hevtsel_emu_->Fill(6,evtwt) ; 

			if (selectedjet.Size >= 4)  
				if (ST_JESLow > 500.) 
					if (btagm == true) 
						m_hevtsel_emu_->Fill(7,evtwt_btag_bcHigh) ; 

			if (selectedjet.Size >= 4)  
				if (ST_JESLow > 500.) 
					if (btagm == true) 
						m_hevtsel_emu_->Fill(8,evtwt_btag_bcLow) ; 

			if (selectedjet.Size >= 4)  
				if (ST_JESLow > 500.) 
					if (btagm == true) 
						m_hevtsel_emu_->Fill(9,evtwt_btag_lHigh) ; 

			if (selectedjet.Size >= 4)  
				if (ST_JESLow > 500.) 
					if (btagm == true) 
						m_hevtsel_emu_->Fill(10,evtwt_btag_lLow) ; 

		}

		if (!isData_ && (inputfiles_.find("BprimeBprimeTo") != std::string::npos) ) { 

			int nbr(0) ; 
			int sumbr(0) ; 
			for (int itw = 0; itw <= 10; ++ itw ) {
				for (int ibz = 0; ibz <= 10; ++ibz) {
					for (int ibh = 0; ibh <= 10; ++ibh) {  
						sumbr = 0 ; 
						sumbr += itw ; 
						sumbr += ibz ; 
						sumbr += ibh ; 
						if ( sumbr == 10 ) {  
							int bprimemode0(EvtInfo.McbprimeMode[0]) ; 
							int bprimemode1(EvtInfo.McbprimeMode[1]) ; 
							double brevtwt             =  evtwt             ; 
							double brevtwt_btag_bcHigh =  evtwt_btag_bcHigh ; 
							double brevtwt_btag_bcLow  =  evtwt_btag_bcLow  ; 
							double brevtwt_btag_lHigh  =  evtwt_btag_lHigh  ; 
							double brevtwt_btag_lLow   =  evtwt_btag_lLow   ; 

							if (bprimemode0 == 1 && bprimemode1 == 1) {
								brevtwt             *= double(itw*itw)*4./100. ; 
								brevtwt_btag_bcHigh *= double(itw*itw)*4./100. ; 
								brevtwt_btag_bcLow  *= double(itw*itw)*4./100. ; 
								brevtwt_btag_lHigh  *= double(itw*itw)*4./100. ; 
								brevtwt_btag_lLow   *= double(itw*itw)*4./100. ; 
							} 
							else if (bprimemode0 == 1 && bprimemode1 == 3) {
								brevtwt             *= double(itw*ibz)*4./100. ; 
								brevtwt_btag_bcHigh *= double(itw*ibz)*4./100. ; 
								brevtwt_btag_bcLow  *= double(itw*ibz)*4./100. ; 
								brevtwt_btag_lHigh  *= double(itw*ibz)*4./100. ; 
								brevtwt_btag_lLow   *= double(itw*ibz)*4./100. ; 
							}
							else if (bprimemode0 == 1 && bprimemode1 == 4) {
								brevtwt             *= double(itw*ibh)*4./100. ; 
								brevtwt_btag_bcHigh *= double(itw*ibh)*4./100. ; 
								brevtwt_btag_bcLow  *= double(itw*ibh)*4./100. ; 
								brevtwt_btag_lHigh  *= double(itw*ibh)*4./100. ; 
								brevtwt_btag_lLow   *= double(itw*ibh)*4./100. ; 
							}
							else if (bprimemode0 == 3 && bprimemode1 == 1) {
								brevtwt             *= double(ibz*itw)*4./100. ; 
								brevtwt_btag_bcHigh *= double(ibz*itw)*4./100. ; 
								brevtwt_btag_bcLow  *= double(ibz*itw)*4./100. ; 
								brevtwt_btag_lHigh  *= double(ibz*itw)*4./100. ; 
								brevtwt_btag_lLow   *= double(ibz*itw)*4./100. ; 
							}
							else if (bprimemode0 == 3 && bprimemode1 == 3) {
								brevtwt             *= double(ibz*ibz)*4./100. ; 
								brevtwt_btag_bcHigh *= double(ibz*ibz)*4./100. ; 
								brevtwt_btag_bcLow  *= double(ibz*ibz)*4./100. ; 
								brevtwt_btag_lHigh  *= double(ibz*ibz)*4./100. ; 
								brevtwt_btag_lLow   *= double(ibz*ibz)*4./100. ; 
							}
							else if (bprimemode0 == 3 && bprimemode1 == 4) {
								brevtwt             *= double(ibz*ibh)*4./100. ; 
								brevtwt_btag_bcHigh *= double(ibz*ibh)*4./100. ; 
								brevtwt_btag_bcLow  *= double(ibz*ibh)*4./100. ; 
								brevtwt_btag_lHigh  *= double(ibz*ibh)*4./100. ; 
								brevtwt_btag_lLow   *= double(ibz*ibh)*4./100. ; 
							}
							else if (bprimemode0 == 4 && bprimemode1 == 1) {
								brevtwt             *= double(ibh*itw)*4./100. ; 
								brevtwt_btag_bcHigh *= double(ibh*itw)*4./100. ; 
								brevtwt_btag_bcLow  *= double(ibh*itw)*4./100. ; 
								brevtwt_btag_lHigh  *= double(ibh*itw)*4./100. ; 
								brevtwt_btag_lLow   *= double(ibh*itw)*4./100. ; 
							}
							else if (bprimemode0 == 4 && bprimemode1 == 3) {
								brevtwt             *= double(ibh*ibz)*4./100. ; 
								brevtwt_btag_bcHigh *= double(ibh*ibz)*4./100. ; 
								brevtwt_btag_bcLow  *= double(ibh*ibz)*4./100. ; 
								brevtwt_btag_lHigh  *= double(ibh*ibz)*4./100. ; 
								brevtwt_btag_lLow   *= double(ibh*ibz)*4./100. ; 
							}
							else if (bprimemode0 == 4 && bprimemode1 == 4) {
								brevtwt             *= double(ibh*ibh)*4./100. ; 
								brevtwt_btag_bcHigh *= double(ibh*ibh)*4./100. ; 
								brevtwt_btag_bcLow  *= double(ibh*ibh)*4./100. ; 
								brevtwt_btag_lHigh  *= double(ibh*ibh)*4./100. ; 
								brevtwt_btag_lLow   *= double(ibh*ibh)*4./100. ; 
							}
							else {
								brevtwt             =  0. ; 
								brevtwt_btag_bcHigh =  0. ; 
								brevtwt_btag_bcLow  =  0. ; 
								brevtwt_btag_lHigh  =  0. ; 
								brevtwt_btag_lLow   =  0. ; 
							}

							if (selectedjet.Size >= 4 && btagm == true) {
								if (isOSSF_el && isZVetoLow)  m_hst_sig[0][0][nbr] -> Fill(ST,brevtwt) ; 
								if (isOSSF_el && isZVetoHigh) m_hst_sig[0][1][nbr] -> Fill(ST,brevtwt) ; 
								if (isOSSF_el && !isZVeto)    m_hst_sig[0][2][nbr] -> Fill(ST,brevtwt) ; 
								if (isOSSF_mu && isZVetoLow)  m_hst_sig[1][0][nbr] -> Fill(ST,brevtwt) ; 
								if (isOSSF_mu && isZVetoHigh) m_hst_sig[1][1][nbr] -> Fill(ST,brevtwt) ; 
								if (isOSSF_mu && !isZVeto)    m_hst_sig[1][2][nbr] -> Fill(ST,brevtwt) ; 
								if (isOSOF && isZVetoLow)     m_hst_sig[2][0][nbr] -> Fill(ST,brevtwt) ; 
								if (isOSOF && isZVetoHigh)    m_hst_sig[2][1][nbr] -> Fill(ST,brevtwt) ; 
								if (isOSOF && !isZVeto)       m_hst_sig[2][2][nbr] -> Fill(ST,brevtwt) ; 
							}

							//// Filling dielectron channel 
							if (isOSSF_el && !isOSSF_mu && !isOSOF) {
								if (selectedjet.Size >= 4)  
									if (ST > 500.) 
										if (btagm == true)			  
											m_hevtselSig_ee_[nbr]->Fill(2,brevtwt) ; 

								if (jet_jesHigh.Size >= 4)  
									if (ST_JESHigh > 500.) 
										if (btagm_jesHigh == true) 
											m_hevtselSig_ee_[nbr]->Fill(3,brevtwt) ; 

								if (jet_jesLow.Size >= 4)  
									if (ST_JESLow > 500.) 
										if (btagm_jesLow == true) 
											m_hevtselSig_ee_[nbr]->Fill(4,brevtwt) ; 

								if (jet_jerHigh.Size >= 4)  
									if (ST_JESHigh > 500.) 
										if (btagm_jerHigh == true) 
											m_hevtselSig_ee_[nbr]->Fill(5,brevtwt) ; 

								if (jet_jerLow.Size >= 4)  
									if (ST_JESLow > 500.) 
										if (btagm_jerLow == true) 
											m_hevtselSig_ee_[nbr]->Fill(6,brevtwt) ; 

								if (selectedjet.Size >= 4)  
									if (ST_JESLow > 500.) 
										if (btagm == true) 
											m_hevtselSig_ee_[nbr]->Fill(7,brevtwt_btag_bcHigh) ; 

								if (selectedjet.Size >= 4)  
									if (ST_JESLow > 500.) 
										if (btagm == true) 
											m_hevtselSig_ee_[nbr]->Fill(8,brevtwt_btag_bcLow) ; 

								if (selectedjet.Size >= 4)  
									if (ST_JESLow > 500.) 
										if (btagm == true) 
											m_hevtselSig_ee_[nbr]->Fill(9,brevtwt_btag_lHigh) ; 

								if (selectedjet.Size >= 4)  
									if (ST_JESLow > 500.) 
										if (btagm == true) 
											m_hevtselSig_ee_[nbr]->Fill(10,brevtwt_btag_lLow) ; 

							}

							//// Filling dimuon channel 
							if (isOSSF_mu && !isOSSF_el && !isOSOF) {
								if (selectedjet.Size >= 4)  
									if (ST > 500.) 
										if (btagm == true)			  
											m_hevtselSig_mumu_[nbr]->Fill(2,brevtwt) ; 

								if (jet_jesHigh.Size >= 4)  
									if (ST_JESHigh > 500.) 
										if (btagm_jesHigh == true) 
											m_hevtselSig_mumu_[nbr]->Fill(3,brevtwt) ; 

								if (jet_jesLow.Size >= 4)  
									if (ST_JESLow > 500.) 
										if (btagm_jesLow == true) 
											m_hevtselSig_mumu_[nbr]->Fill(4,brevtwt) ; 

								if (jet_jerHigh.Size >= 4)  
									if (ST_JESHigh > 500.) 
										if (btagm_jerHigh == true) 
											m_hevtselSig_mumu_[nbr]->Fill(5,brevtwt) ; 

								if (jet_jerLow.Size >= 4)  
									if (ST_JESLow > 500.) 
										if (btagm_jerLow == true) 
											m_hevtselSig_mumu_[nbr]->Fill(6,brevtwt) ; 

								if (selectedjet.Size >= 4)  
									if (ST_JESLow > 500.) 
										if (btagm == true) 
											m_hevtselSig_mumu_[nbr]->Fill(7,brevtwt_btag_bcHigh) ; 

								if (selectedjet.Size >= 4)  
									if (ST_JESLow > 500.) 
										if (btagm == true) 
											m_hevtselSig_mumu_[nbr]->Fill(8,brevtwt_btag_bcLow) ; 

								if (selectedjet.Size >= 4)  
									if (ST_JESLow > 500.) 
										if (btagm == true) 
											m_hevtselSig_mumu_[nbr]->Fill(9,brevtwt_btag_lHigh) ; 

								if (selectedjet.Size >= 4)  
									if (ST_JESLow > 500.) 
										if (btagm == true) 
											m_hevtselSig_mumu_[nbr]->Fill(10,brevtwt_btag_lLow) ; 

							}

							//// Filling emu channel 
							if (isOSOF && !isOSSF_el && !isOSSF_mu) {
								if (selectedjet.Size >= 4)  
									if (ST > 500.) 
										if (btagm == true)			  
											m_hevtselSig_emu_[nbr]->Fill(2,brevtwt) ; 

								if (jet_jesHigh.Size >= 4)  
									if (ST_JESHigh > 500.) 
										if (btagm_jesHigh == true) 
											m_hevtselSig_emu_[nbr]->Fill(3,brevtwt) ; 

								if (jet_jesLow.Size >= 4)  
									if (ST_JESLow > 500.) 
										if (btagm_jesLow == true) 
											m_hevtselSig_emu_[nbr]->Fill(4,brevtwt) ; 

								if (jet_jerHigh.Size >= 4)  
									if (ST_JESHigh > 500.) 
										if (btagm_jerHigh == true) 
											m_hevtselSig_emu_[nbr]->Fill(5,brevtwt) ; 

								if (jet_jerLow.Size >= 4)  
									if (ST_JESLow > 500.) 
										if (btagm_jerLow == true) 
											m_hevtselSig_emu_[nbr]->Fill(6,brevtwt) ; 

								if (selectedjet.Size >= 4)  
									if (ST_JESLow > 500.) 
										if (btagm == true) 
											m_hevtselSig_emu_[nbr]->Fill(7,brevtwt_btag_bcHigh) ; 

								if (selectedjet.Size >= 4)  
									if (ST_JESLow > 500.) 
										if (btagm == true) 
											m_hevtselSig_emu_[nbr]->Fill(8,brevtwt_btag_bcLow) ; 

								if (selectedjet.Size >= 4)  
									if (ST_JESLow > 500.) 
										if (btagm == true) 
											m_hevtselSig_emu_[nbr]->Fill(9,brevtwt_btag_lHigh) ; 

								if (selectedjet.Size >= 4)  
									if (ST_JESLow > 500.) 
										if (btagm == true) 
											m_hevtselSig_emu_[nbr]->Fill(10,brevtwt_btag_lLow) ; 

							}

							++nbr ; 

						} //// Correct BRs: itw+ibz+ibh = 1 
					} //// ibh loop 
				} //// ibz loop 
			} //// itw loop  
		} //// if is signal

	} //// End of event loop 

	std::cout << " End of event loop\n" ; 

	fntuple->Write() ; 
	fntuple->Close() ; 

} //// BprimeCandSelector::evtLoop () 

