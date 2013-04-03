// -*- C++ -*-
//
// Package:    HiggsTagValidation
// Class:      HiggsTagValidation
// 
/**\class HiggsTagValidation HiggsTagValidation.cc HiggsTagging/HiggsTagValidation/src/HiggsTagValidation.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Devdatta Majumder,13 2-054,+41227671675,
//         Created:  Wed Mar 20 10:59:16 CET 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/Track.h" 
#include "DataFormats/TrackReco/interface/TrackFwd.h" 
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TH1D.h"
#include "TH2D.h"

#include "HiggsTagging/HiggsTagValidation/src/CommHistProducer.cc"

//
// class declaration
//

class HiggsTagValidation : public edm::EDAnalyzer {
	public:
		explicit HiggsTagValidation(const edm::ParameterSet&);
		~HiggsTagValidation();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

		// typedefs
		typedef std::vector<pat::Jet> PatJetCollection;

	private:
		virtual void beginJob() ;
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob() ;

		virtual void beginRun(edm::Run const&, edm::EventSetup const&);
		virtual void endRun(edm::Run const&, edm::EventSetup const&);
		virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
		virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

		int TaggedJet(pat::Jet , edm::Handle<reco::JetTagCollection> ); 

		// ----------member data ---------------------------

		edm::Service<TFileService> fs;

		// ----------Configurables--------------------------
		const bool              m_isData; 
		const bool              m_useEventWeight;
		const int               m_bosonPdgId;
		const edm::InputTag     m_ipTrackTag;  
		const edm::InputTag     m_genParticleTag ; 
		const edm::InputTag     m_jetsTag ; 
		const edm::InputTag     m_groomedJetsTag ; 
		const edm::InputTag     m_groomedBasicJetsTag ;  
		const edm::InputTag     m_subJetsTag ; 
		const edm::InputTag     m_pvTag ; 
		const double            m_jetRadius;
		const std::vector<int>  m_jetFlavourPdgId ; 
		const double            m_jetPtMin;
		const double            m_jetAbsEtaMax;
		const double            m_jetMassMin;
		const double            m_jetMassMax; 
		const double            m_dRSubjetMatch;  

		TH1D *h1_BosonPt         ; 

		TH1D* h1_All_NJets           ; 
		TH1D* h1_All_NJetsMatched    ;
		TH1D* h1_All_NSubjets        ;
		TH1D* h1_All_NJetsBtagged    ;
		TH1D* h1_All_NSubjetsBtagged ;

		TH1D* h1_All_JetPt       ;
		TH1D* h1_All_JetEta      ;
		TH1D* h1_All_JetMass     ;
		TH1D* h1_All_NTracks     ;
		TH1D *h1_All_csvBDisc    ; 

		std::map<int,TH1D*> h1_Flav_JetPt     ;
		std::map<int,TH1D*> h1_Flav_JetEta    ;
		std::map<int,TH1D*> h1_Flav_JetMass   ;
		std::map<int,TH1D*> h1_Flav_NTracks   ;
		std::map<int,TH1D*> h1_Flav_csvBDisc  ; 

		TH1D* h1_JetPt           ;
		TH1D* h1_JetEta          ;
		TH1D* h1_JetMass         ;
		TH1D* h1_JetMassDrop     ;
		TH1D* h1_JetMassDropOrg  ;
		TH1D* h1_NTracks         ;
		TH1D *h1_csvBDisc        ; 

		TH1D* h1_Tagger          ;
		TH1D* h1_Tagger_TCHE     ;
		TH1D* h1_Tagger_TCHP     ;
		TH1D* h1_Tagger_JP	 ;
		TH1D* h1_Tagger_SSVHE    ;
		TH1D* h1_Tagger_SSVHP    ;
		TH1D* h1_Tagger_CSV      ;
		TH1D* h1_Tagger_MU       ;

		TH1D* hAllFlav_Flavour      ;
		TH1D* hAllFlav_Tagger       ;
		TH1D* hAllFlav_Tagger_Gam   ;
		TH1D* hAllFlav_Tagger_K0s   ;
		TH1D* hAllFlav_Tagger_Lam   ;
		TH1D* hAllFlav_Tagger_Bwd   ;
		TH1D* hAllFlav_Tagger_Cwd   ;
		TH1D* hAllFlav_Tagger_Tau   ;
		TH1D* hAllFlav_Tagger_Int   ;
		TH1D* hAllFlav_Tagger_Fak   ;
		TH1D* hAllFlav_Tagger_Bad   ;
		TH1D* hAllFlav_Tagger_Oth   ;

		TH1D* hLightFlav_Tagger     ;
		TH1D* hGluonFlav_Tagger     ;
		TH1D* hUDSFlav_Tagger       ;
		TH1D* hCFlav_Tagger         ;
		TH1D* hBFlav_Tagger         ;

		TH1D*  IPSign_cat0   ; 
		TH1D*  IPSign_cat1   ;
		TH1D*  IPSign_cat2   ;
		TH1D*  IPSign_cat3   ;
		TH1D*  IPSign_cat4   ;
		TH1D*  IPSign_cat5   ;
		TH1D*  IPSign_cat6   ;
		TH1D*  IPSign_cat7   ;
		TH1D*  IPSign_cat8   ;
		TH1D*  IPSign_cat9   ;
		TH1D*  TrackProbaNeg      ; 
		TH1D*  TrackProbaNeg_Cat0 ;
		TH1D*  TrackProbaNeg_Cat1 ;
		TH1D*  TrackProbaNeg_Cat2 ;
		TH1D*  TrackProbaNeg_Cat3 ;
		TH1D*  TrackProbaNeg_Cat4 ;
		TH1D*  TrackProbaNeg_Cat5 ;
		TH1D*  TrackProbaNeg_Cat6 ;
		TH1D*  TrackProbaNeg_Cat7 ;
		TH1D*  TrackProbaNeg_Cat8 ;
		TH1D*  TrackProbaNeg_Cat9 ;
		TH1D*  TrackProbJet80      ; 
		TH1D*  TrackProbJet80_Cat0 ;
		TH1D*  TrackProbJet80_Cat1 ;
		TH1D*  TrackProbJet80_Cat2 ;
		TH1D*  TrackProbJet80_Cat3 ;
		TH1D*  TrackProbJet80_Cat4 ;
		TH1D*  TrackProbJet80_Cat5 ;
		TH1D*  TrackProbJet80_Cat6 ;
		TH1D*  TrackProbJet80_Cat7 ;
		TH1D*  TrackProbJet80_Cat8 ;
		TH1D*  TrackProbJet80_Cat9 ;

		CommHistProducer histProducer ; 

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HiggsTagValidation::HiggsTagValidation(const edm::ParameterSet& iConfig) : 

	m_isData(iConfig.getParameter<bool>("IsData")), 
	m_useEventWeight(iConfig.getParameter<bool>("UseEventWeight")),
	m_bosonPdgId(iConfig.getParameter<int>("BosonPdgId")),
	m_ipTrackTag(iConfig.getParameter<edm::InputTag>("IpTrackTag")), 
	m_genParticleTag(iConfig.getParameter<edm::InputTag>("GenParticleTag")),
	m_jetsTag(iConfig.getParameter<edm::InputTag>("JetsTag")),
	m_groomedJetsTag(iConfig.getParameter<edm::InputTag>("GroomedJetsTag")),
	//  m_groomedBasicJetsTag(iConfig.getParameter<edm::InputTag>("GroomedBasicJetsTag")),
	m_subJetsTag(iConfig.getParameter<edm::InputTag>("SubJetsTag")),
	m_pvTag(iConfig.getParameter<edm::InputTag>("PvTag")), 
	m_jetRadius(iConfig.getParameter<double>("JetRadius")),
	m_jetFlavourPdgId(iConfig.getParameter<std::vector<int> >("JetFlavourPdgId")), 
	m_jetPtMin(iConfig.getParameter<double>("JetPtMin")),
	m_jetAbsEtaMax(iConfig.getParameter<double>("JetAbsEtaMax")),
	m_jetMassMin(iConfig.getParameter<double>("JetMassMin")),
	m_jetMassMax(iConfig.getParameter<double>("JetMassMax")), 
	m_dRSubjetMatch(iConfig.getParameter<double>("DRSubjetMatch")),
	histProducer(fs,m_isData)

{
	//now do what ever initialization is needed
	
	int ptBins=250;
	double ptMin=0., ptMax=1000. ;
	float AbsEtaMax = m_jetAbsEtaMax ; 

	histProducer.AddHisto("jet_multi"    ,"number of jets",                 20,0,20    );
	histProducer.AddHisto("jet_pt_all"   ,"pT of all jets",                 ptMax/10,0,ptMax);
	histProducer.AddHisto("jet_pt_sv"    ,"pT of jets containing a SV",     ptMax/10,0,ptMax);
	histProducer.AddHisto("jet_eta_all"  ,"eta of all jets",                50,-2.5,2.5);
	histProducer.AddHisto("jet_phi_all"  ,"phi of all jets",                20,-5,5    );
	histProducer.AddHisto("jet_mass_all"  ,"mass of all jets",              ptMax/10,0,ptMax);

	histProducer.AddHisto("jet_track_multi"  ,      "number of tracks in the jets",                40,-0.5,39.5  );
	histProducer.AddHisto("jet_track_multi_sel"  ,  "number of selected tracks in the jets",       40,-0.5,39.5  );
	histProducer.AddHisto("jet_track_charge"   ,    "charge of the tracks",                      2,-2.,2.);
	histProducer.AddHisto("jet_track_chi2"   ,      "normalized chi2 of the tracks",               100,0.,30.    );
	histProducer.AddHisto("jet_track_nHit" ,        "number of hits ",               35,-0.5, 34.5 );
	histProducer.AddHisto("jet_track_HPix"   ,      "number of hits in the Pixel",                 10,-0.5, 9.5  );

	histProducer.AddHisto("jet_ipTagTrack_dxy",      "Dxy of impactParameterTag tracks",     100,-1.,1.); 
	histProducer.AddHisto("jet_ipTagTrack_dz",       "Dz  of impactParameterTag tracks",     100,-10.,10.); 
	histProducer.AddHisto("jet_ipTagTrack_zIP",       "z IP  of impactParameterTag tracks",     100,-10.,10.); 

	histProducer.AddHisto("jet_ipTagNSV"    ,"No. of subjets",                                    20,0,20    ); 

	histProducer.AddHisto("jet_csvBDisc", "Jet CSV b discriminator", 200, 0.,1.) ; 

	histProducer.AddHisto("subjet_multi"    ,"number of subjets",                                    20,0,20    );
	histProducer.AddHisto("subjet_pt_all"   ,"pT of all subjets",                                    ptMax/10,0,ptMax );
	histProducer.AddHisto("subjet_pt_sv"    ,"pT of subjets containing a SV",                        ptMax/10,0,ptMax);
	histProducer.AddHisto("subjet_eta_all"  ,"eta of all subjets",                                   50,-2.5,2.5);
	histProducer.AddHisto("subjet_phi_all"  ,"phi of all subjets",                                   20,-5,5    );
	histProducer.AddHisto("subjet_mass_all"  ,"mass of all subjets",                                 ptMax/10,0,ptMax);
	histProducer.AddHisto("subjet_massdrop_all"  ,"mass drop of all subjets",                        40,0.,2.);
	histProducer.AddHisto("subjet_massdropOrg_all"  ,"mass drop w.r.t original jet of all subjets",  20,0.,2.);

	histProducer.AddHisto("subjet_track_multi"  ,      "number of tracks in the jets",       40,-0.5,39.5  );
	histProducer.AddHisto("subjet_track_multi_sel"  ,  "number of selected tracks in the jets",       40,-0.5,39.5  );
	histProducer.AddHisto("subjet_track_charge"   ,    "charge of the tracks",                      2,-2.,2.);
	histProducer.AddHisto("subjet_track_chi2"   ,      "normalized chi2 of the tracks",               100,0.,30.    );
	histProducer.AddHisto("subjet_track_nHit" ,        "number of hits ",               35,-0.5, 34.5 );
	histProducer.AddHisto("subjet_track_HPix"   ,      "number of hits in the Pixel",                 10,-0.5, 9.5  );

	histProducer.AddHisto("subjet_csvBDisc", "Subjet CSV b discriminator", 200, 0.,1.) ; 

	h1_BosonPt         = fs->make<TH1D>("h1_BosonPt",";p_{T} [GeV];",ptBins,ptMin,ptMax);

	h1_All_NJets           = fs->make<TH1D>("h1_All_NJets","nb. of jets"                                          ,21, -0.5, 20.5);
	h1_All_NJetsMatched    = fs->make<TH1D>("h1_All_NJetsMatched","nb. of jets matched to groomed jets"           ,21, -0.5, 20.5);
	h1_All_NSubjets        = fs->make<TH1D>("h1_All_NSubjets","nb. of subjets"                                    ,21, -0.5, 20.5);
	h1_All_NJetsBtagged    = fs->make<TH1D>("h1_All_NJetsBtagged","nb. of b-tagged jets matched to groomed jets"  ,21, -0.5, 20.5);
	h1_All_NSubjetsBtagged = fs->make<TH1D>("h1_All_NSubjetsBtagged","nb. of b-tagged subjets"                    ,21, -0.5, 20.5);

	h1_All_JetPt       = fs->make<TH1D>("h1_All_JetPt","pt(jet)"                                 ,50,  20., 520.);
	h1_All_JetMass     = fs->make<TH1D>("h1_All_JetMass","M(jet)"                                ,100,  20.,1020.);
	h1_All_JetEta      = fs->make<TH1D>("h1_All_JetEta","|#eta(jet)|"                            ,50, -2.5, 2.5 );
	h1_All_NTracks     = fs->make<TH1D>("h1_All_NTracks","nb. of tracks"                         ,101,-0.5, 100.5);
	h1_All_csvBDisc    = fs->make<TH1D>("h1_All_csvBDisc",";CSV b tag discriminator;"            ,100, 0.,  1.) ; 

	h1_JetPt           = fs->make<TH1D>("h1_JetPt","pt(jet)"                                     ,50,  20., 520.);
	h1_JetEta          = fs->make<TH1D>("h1_JetEta","|#eta(jet)|"                                ,50,-2.5 , 2.5 );
	h1_JetMass         = fs->make<TH1D>("h1_JetMass","M(jet)"                                    ,100,  20.,1020.);
	h1_JetMassDrop     = fs->make<TH1D>("h1_JetMassDrop","M(subjet)/M(jet)"                      ,50,    0.,   2.);
	h1_JetMassDropOrg  = fs->make<TH1D>("h1_JetMassDropOrg","M(subjet)/M(org jet)"               ,50,    0.,   2.);
	h1_NTracks         = fs->make<TH1D>("h1_NTracks","nb. of tracks"                   ,101,-0.5,100.5);
	h1_csvBDisc        = fs->make<TH1D>("h1_csvBDisc",";CSV b tag discriminator;",100,0.,1.) ; 

	h1_Tagger          = fs->make<TH1D>("h1_Tagger","Tagger"             ,100,-25.,25.);
	h1_Tagger_TCHE	   = fs->make<TH1D>("h1_Tagger_TCHE","Tagger_TCHE"   ,100,-25.,25.);
	h1_Tagger_TCHP	   = fs->make<TH1D>("h1_Tagger_TCHP","Tagger_TCHP"   ,100,-25.,25.);
	h1_Tagger_JP       = fs->make<TH1D>("h1_Tagger_JP","Tagger_JP"       ,100,-25.,25.);
	h1_Tagger_SSVHE	   = fs->make<TH1D>("h1_Tagger_SSVHE ","Tagger_SSVHE",100,-25.,25.);
	h1_Tagger_SSVHP	   = fs->make<TH1D>("h1_Tagger_SSVHP","Tagger_SSVHP" ,100,-25.,25.);
	h1_Tagger_CSV	   = fs->make<TH1D>("h1_Tagger_CSV","Tagger_CSV"     ,100,-25.,25.);
	h1_Tagger_MU       = fs->make<TH1D>("h1_Tagger_MU","Tagger_MU"       ,100,-25.,25.);

	hAllFlav_Flavour        = fs->make<TH1D>("hAllFlav_Flavour","Flavour"   ,22,-0.5,21.5);
	hAllFlav_Tagger	        = fs->make<TH1D>("hAllFlav_Tagger","Tagger"     ,100,-25.,25.);
	hAllFlav_Tagger_Gam	= fs->make<TH1D>("hAllFlav_Tagger_Gam","Tagger" ,100,-25.,25.);
	hAllFlav_Tagger_K0s	= fs->make<TH1D>("hAllFlav_Tagger_K0s","Tagger" ,100,-25.,25.);
	hAllFlav_Tagger_Lam	= fs->make<TH1D>("hAllFlav_Tagger_Lam","Tagger" ,100,-25.,25.);
	hAllFlav_Tagger_Bwd	= fs->make<TH1D>("hAllFlav_Tagger_Bwd","Tagger" ,100,-25.,25.);
	hAllFlav_Tagger_Cwd	= fs->make<TH1D>("hAllFlav_Tagger_Cwd","Tagger" ,100,-25.,25.);
	hAllFlav_Tagger_Tau	= fs->make<TH1D>("hAllFlav_Tagger_Tau","Tagger" ,100,-25.,25.);
	hAllFlav_Tagger_Int	= fs->make<TH1D>("hAllFlav_Tagger_Int","Tagger" ,100,-25.,25.);
	hAllFlav_Tagger_Fak	= fs->make<TH1D>("hAllFlav_Tagger_Fak","Tagger" ,100,-25.,25.);
	hAllFlav_Tagger_Bad	= fs->make<TH1D>("hAllFlav_Tagger_Bad","Tagger" ,100,-25.,25.);
	hAllFlav_Tagger_Oth	= fs->make<TH1D>("hAllFlav_Tagger_Oth","Tagger" ,100,-25.,25.);

	hLightFlav_Tagger	= fs->make<TH1D>("hLightFlav_Tagger","Tagger",100,-25.,25.);
	hGluonFlav_Tagger	= fs->make<TH1D>("hGluonFlav_Tagger","Tagger",100,-25.,25.);
	hUDSFlav_Tagger	        = fs->make<TH1D>("hUDSFlav_Tagger","Tagger"  ,100,-25.,25.);
	hCFlav_Tagger 	        = fs->make<TH1D>("hCFlav_Tagger","Tagger"    ,100,-25.,25.);
	hBFlav_Tagger 	        = fs->make<TH1D>("hBFlav_Tagger","Tagger"    ,100,-25.,25.);

	IPSign_cat0		= fs->make<TH1D>("IPSign_cat0","-IP/#sigma",50, 0., 25.);
	IPSign_cat1		= fs->make<TH1D>("IPSign_cat1","-IP/#sigma",50, 0., 25.);
	IPSign_cat2		= fs->make<TH1D>("IPSign_cat2","-IP/#sigma",50, 0., 25.);
	IPSign_cat3		= fs->make<TH1D>("IPSign_cat3","-IP/#sigma",50, 0., 25.);
	IPSign_cat4		= fs->make<TH1D>("IPSign_cat4","-IP/#sigma",50, 0., 25.);
	IPSign_cat5		= fs->make<TH1D>("IPSign_cat5","-IP/#sigma",50, 0., 25.);
	IPSign_cat6		= fs->make<TH1D>("IPSign_cat6","-IP/#sigma",50, 0., 25.);
	IPSign_cat7		= fs->make<TH1D>("IPSign_cat7","-IP/#sigma",50, 0., 25.);
	IPSign_cat8		= fs->make<TH1D>("IPSign_cat8","-IP/#sigma",50, 0., 25.);
	IPSign_cat9		= fs->make<TH1D>("IPSign_cat9","-IP/#sigma",50, 0., 25.);
	TrackProbaNeg 	        = fs->make<TH1D>("TrackProbaNEG","TrackProbaNEG",51, 0., 1.02);
	TrackProbaNeg_Cat0	= fs->make<TH1D>("TrackProbaNEG Cat0","TrackProbaNEG Cat0",51, 0., 1.02);
	TrackProbaNeg_Cat1	= fs->make<TH1D>("TrackProbaNEG Cat1","TrackProbaNEG Cat1",51, 0., 1.02);
	TrackProbaNeg_Cat2	= fs->make<TH1D>("TrackProbaNEG Cat2","TrackProbaNEG Cat2",51, 0., 1.02);
	TrackProbaNeg_Cat3	= fs->make<TH1D>("TrackProbaNEG Cat3","TrackProbaNEG Cat3",51, 0., 1.02);
	TrackProbaNeg_Cat4	= fs->make<TH1D>("TrackProbaNEG Cat4","TrackProbaNEG Cat4",51, 0., 1.02);
	TrackProbaNeg_Cat5	= fs->make<TH1D>("TrackProbaNEG Cat5","TrackProbaNEG Cat5",51, 0., 1.02);
	TrackProbaNeg_Cat6	= fs->make<TH1D>("TrackProbaNEG Cat6","TrackProbaNEG Cat6",51, 0., 1.02);
	TrackProbaNeg_Cat7	= fs->make<TH1D>("TrackProbaNEG Cat7","TrackProbaNEG Cat7",51, 0., 1.02);
	TrackProbaNeg_Cat8	= fs->make<TH1D>("TrackProbaNEG Cat8","TrackProbaNEG Cat8",51, 0., 1.02);
	TrackProbaNeg_Cat9	= fs->make<TH1D>("TrackProbaNEG Cat9","TrackProbaNEG Cat9",51, 0., 1.02);
	TrackProbJet80      	= fs->make<TH1D>("TrackProbJet80","TrackProbJet80",51, 0., 1.02);
	TrackProbJet80_Cat0	= fs->make<TH1D>("TrackProbJet80 Cat0","TrackProbJet80 Cat0",51, 0., 1.02);
	TrackProbJet80_Cat1	= fs->make<TH1D>("TrackProbJet80 Cat1","TrackProbJet80 Cat1",51, 0., 1.02);
	TrackProbJet80_Cat2	= fs->make<TH1D>("TrackProbJet80 Cat2","TrackProbJet80 Cat2",51, 0., 1.02);
	TrackProbJet80_Cat3	= fs->make<TH1D>("TrackProbJet80 Cat3","TrackProbJet80 Cat3",51, 0., 1.02);
	TrackProbJet80_Cat4	= fs->make<TH1D>("TrackProbJet80 Cat4","TrackProbJet80 Cat4",51, 0., 1.02);
	TrackProbJet80_Cat5	= fs->make<TH1D>("TrackProbJet80 Cat5","TrackProbJet80 Cat5",51, 0., 1.02);
	TrackProbJet80_Cat6	= fs->make<TH1D>("TrackProbJet80 Cat6","TrackProbJet80 Cat6",51, 0., 1.02);
	TrackProbJet80_Cat7	= fs->make<TH1D>("TrackProbJet80 Cat7","TrackProbJet80 Cat7",51, 0., 1.02);
	TrackProbJet80_Cat8	= fs->make<TH1D>("TrackProbJet80 Cat8","TrackProbJet80 Cat8",51, 0., 1.02);
	TrackProbJet80_Cat9	= fs->make<TH1D>("TrackProbJet80 Cat9","TrackProbJet80 Cat9",51, 0., 1.02);

}


HiggsTagValidation::~HiggsTagValidation()   
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
	void
HiggsTagValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

	typedef reco::LeafCandidate::LorentzVector LorentzVector;

	double eventWeight = 1.;
	if( !iEvent.isRealData() && m_useEventWeight ) {
		edm::Handle<GenEventInfoProduct> genEvtInfoProduct;
		iEvent.getByLabel("generator", genEvtInfoProduct);

		eventWeight = genEvtInfoProduct->weight();
	}

	edm::LogInfo("EventWeight") << "Event weight = " << eventWeight ;

	edm::Handle<reco::GenParticleCollection> genParticles;

	std::vector<LorentzVector> bFromGSplit4v ; 
	if (!m_isData) {
		iEvent.getByLabel(m_genParticleTag,genParticles);
		//// Identify b quarks from GSF or hard scattering 
		bool isGluonSplitting = false;
		bool isMatrixElement = false;
		bool bFoundS3Quark = false;
		bool bFoundS2Quark = false;
		for(reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it) {
			if( abs(it->pdgId()) == abs(m_bosonPdgId) && it->status() == 3 ) {
				h1_BosonPt->Fill( it->pt(), eventWeight );
			}

			int pdgid = abs(it->pdgId());
			// skip GenParticle if not B hadron
			//if ( !((pdgid/100)%10 == 5 || (pdgid/1000)%10 == 5) ) continue;
			if( abs(it->pdgId()) == 5 && it->status()==3 ) bFoundS3Quark = true;
			if( abs(it->pdgId()) == 5 && it->status()==2 ) { bFoundS2Quark = true; bFromGSplit4v.push_back(it->p4()) ; } 
		}
		//// if no status 3 b quark but status 2
		if( (!bFoundS3Quark) && bFoundS2Quark) isGluonSplitting = true;
		//// if status 3 b quark found
		if( bFoundS3Quark ) isMatrixElement = true;

		//// skip event if it does not pass the selection
		//if( (!isGluonSplitting) || (!isMatrixElement) ) return;
	} 

	//// PV selection 
	Handle<reco::VertexCollection> primaryVertex;
	iEvent.getByLabel(m_pvTag,primaryVertex);
	const  reco::Vertex* pvtx ;
	bool pvFound = (primaryVertex->size() != 0);
	if ( pvFound ) {           
		pvtx = &(*primaryVertex->begin());
	} 
	else { 
		reco::Vertex::Error err;   
		err(0,0)=0.0015*0.0015;    
		err(1,1)=0.0015*0.0015;    
		err(2,2)=15.*15.;
		reco::Vertex::Point pt(0,0,0);
		pvtx=  new reco::Vertex(pt,err,1,1,1);
		//newvertex = true;
	}
	GlobalPoint pv_point = GlobalPoint((*pvtx).x(), (*pvtx).y(), (*pvtx).z());

	//// Track selection 
	Handle<std::vector<reco::TrackIPTagInfo> > tagInfo;
	iEvent.getByLabel(m_ipTrackTag, tagInfo);

	//// Jet selection 
	edm::Handle<PatJetCollection> jets;
	iEvent.getByLabel(m_jetsTag,jets);

	edm::Handle<PatJetCollection> groomedJets;
	iEvent.getByLabel(m_groomedJetsTag,groomedJets);

	edm::Handle<PatJetCollection> subJets;
	iEvent.getByLabel(m_subJetsTag,subJets);

	edm::LogInfo("NJets") << "NJets = " << jets->size() ;

	//// Loop over all jets 
	int njets(0) ; 
	int njetsMatched(0) ; 
	int njetsBtagged(0) ; 
	int njets_b(0) ; 
	int njets_bfromg(0) ; 
	int njets_c(0) ; 
	int njets_l(0) ; 

	int nsubjets_b(0) ; 
	int nsubjets_bfromg(0) ; 
	int nsubjets_c(0) ; 
	int nsubjets_l(0) ; 
	for(PatJetCollection::const_iterator it = jets->begin(); it != jets->end(); ++it) {

		int jetflavour(0) ; 
		bool bjetFromGSplit(false) ; 

		double jetPt = it->pt();
		double jetEta = fabs(it->eta()) ; 
		//// skip the jet if it does not pass pT and eta cuts
		if( !(jetPt > m_jetPtMin && jetEta < m_jetAbsEtaMax) ) continue;

		if (!m_isData) {
			//// Check jet flavour 
			bool isRightFlavour = false; 
			for(std::vector<int>::const_iterator pdgIdIt = m_jetFlavourPdgId.begin(); pdgIdIt != m_jetFlavourPdgId.end(); ++pdgIdIt) {
				if( abs(it->partonFlavour()) == abs(*pdgIdIt)) { 
					isRightFlavour = true; 
					break;
				}
				jetflavour = abs(it->partonFlavour()) ; 
			}
			//// skip the jet if it does not have the right flavour
			if( !isRightFlavour ) continue;
		}

		//// Find jets matched to groomed jets 
		PatJetCollection::const_iterator groomedJetMatch;
		bool groomedJetMatchFound = false;
		double dR = m_jetRadius;
		for (PatJetCollection::const_iterator gjIt = groomedJets->begin(); gjIt != groomedJets->end(); ++gjIt) {
			double dR_temp = reco::deltaR( it->p4(), gjIt->p4() );
			if( dR_temp < dR ) { 
				groomedJetMatchFound = true;
				dR = dR_temp;
				groomedJetMatch = gjIt;
			}
		}
		if( groomedJetMatchFound ) { 

			++njetsMatched ; 

			// vector of pointers to subjets
			std::vector<const pat::Jet*> subjets;
			for (unsigned ndau = 0; ndau < groomedJetMatch->numberOfDaughters(); ++ndau) {			
				pat::Jet const * subjet = dynamic_cast<pat::Jet const *>(groomedJetMatch->daughter(ndau));
				subjets.push_back(subjet) ; 
			}

			h1_All_NSubjets -> Fill (subjets.size(), eventWeight) ; 

			edm::LogInfo("NSubjets") << "Subjets found = " << subjets.size() ;
			if( subjets.size() < 2 ) edm::LogError("TooFewSubjets") << "Less than two subjets found.";
			else if( subjets.size() > 2 ) edm::LogError("TooManySubjets") << "More than two subjets found.";
			else { //// These are the jets we want to use 

				//// Do validation 

				if (!m_isData) {
					edm::LogInfo("JetFlavour") << " jet flavour is " << jetflavour ; 

					//// Check if jet is from gluon splitting 
					int nBFromGSplit(0) ; 
					for (std::vector<LorentzVector>::const_iterator itSplitBs = bFromGSplit4v.begin(); itSplitBs != bFromGSplit4v.end(); ++itSplitBs) {
						if (reco::deltaR(it->p4(), *itSplitBs) < m_dRSubjetMatch) ++nBFromGSplit ; 
					}
					if (jetflavour == 5 && nBFromGSplit >= 2) {
						bjetFromGSplit = true ; 
						edm::LogInfo("BFromGSplit") << " Identified " ; 
					}
					else {
						bjetFromGSplit = false ; 
						edm::LogInfo("BFromGSplit") << " Not Identified " ; 
					} 
					++njets ; 
					if (jetflavour == 5 && bjetFromGSplit == false) ++njets_b ; 
					if (jetflavour == 5 && bjetFromGSplit == true) ++njets_bfromg ; 
					if (jetflavour == 4 && bjetFromGSplit == false) ++njets_c ; 
					if ( (jetflavour < 4 || jetflavour == 21)  && bjetFromGSplit == false) ++njets_l ; 

				}

				histProducer.FillHisto_floatFromMap("jet_pt_all",   jetflavour, bjetFromGSplit, groomedJetMatch->correctedP4(0).pt(),  eventWeight) ; 
				histProducer.FillHisto_floatFromMap("jet_eta_all",  jetflavour, bjetFromGSplit, groomedJetMatch->correctedP4(0).eta(), eventWeight) ; 
				histProducer.FillHisto_floatFromMap("jet_phi_all",  jetflavour, bjetFromGSplit, groomedJetMatch->correctedP4(0).phi(), eventWeight) ; 
				histProducer.FillHisto_floatFromMap("jet_mass_all", jetflavour, bjetFromGSplit, groomedJetMatch->correctedP4(0).mag(), eventWeight) ; 

				h1_All_JetPt -> Fill( groomedJetMatch->correctedP4(0).pt(), eventWeight) ; 
				h1_All_JetEta -> Fill( groomedJetMatch->correctedP4(0).eta(), eventWeight) ; 
				h1_All_JetMass -> Fill (groomedJetMatch->correctedP4(0).mag(), eventWeight) ; 
				const reco::TrackRefVector tracks = groomedJetMatch->associatedTracks() ; 
				int ntracks = tracks.size() ; 
				h1_All_NTracks -> Fill(ntracks, eventWeight) ; 
				double jetBdisc = groomedJetMatch->bDiscriminator("combinedSecondaryVertexBJetTags") ;  
				h1_All_csvBDisc -> Fill(jetBdisc, eventWeight) ; 
				if (jetBdisc > 0.679) ++njetsBtagged ; 

				histProducer.FillHisto_floatFromMap("jet_csvBDisc",  jetflavour, bjetFromGSplit, jetBdisc, eventWeight) ; 

				const reco::TrackRefVector trks_jets = groomedJetMatch->associatedTracks() ; 
				histProducer.FillHisto_floatFromMap("jet_track_multi",  jetflavour, bjetFromGSplit, trks_jets.size(),  eventWeight) ;

				for (reco::track_iterator ittrk_jets = trks_jets.begin(); ittrk_jets != trks_jets.end(); ++ittrk_jets) {
					histProducer.FillHisto_floatFromMap("jet_track_chi2",  jetflavour, bjetFromGSplit, (*ittrk_jets)->normalizedChi2(),                        eventWeight) ;
					histProducer.FillHisto_floatFromMap("jet_track_charge",jetflavour, bjetFromGSplit, (*ittrk_jets)->charge(),                                eventWeight) ;
					histProducer.FillHisto_floatFromMap("jet_track_nHit",  jetflavour, bjetFromGSplit, (*ittrk_jets)->numberOfValidHits(),                     eventWeight) ;
					histProducer.FillHisto_floatFromMap("jet_track_HPix",  jetflavour, bjetFromGSplit, ((*ittrk_jets)->hitPattern()).numberOfValidPixelHits(), eventWeight) ;
				}

				////// Loop over tracks
				if ( groomedJetMatch->hasTagInfo("impactParameter") ) { 
					edm::LogInfo("impactParameterTagInfos") << "Found" ; 
					const edm::RefVector<reco::TrackCollection> &selTrks_svTagInfos_jets( (groomedJetMatch->tagInfoTrackIP("impactParameter"))->selectedTracks() ) ; 
					int nSelectedTracks = (it->hasTagInfo("impactParameter") ? groomedJetMatch->tagInfoTrackIP("impactParameter")->selectedTracks().size() : -99.);
					for (unsigned int itrk = 0; itrk < selTrks_svTagInfos_jets.size(); ++itrk) {
						double Track_dxy = selTrks_svTagInfos_jets[itrk]->dxy(pvtx->position());
						histProducer.FillHisto_floatFromMap("jet_ipTagTrack_dxy", jetflavour, bjetFromGSplit, selTrks_svTagInfos_jets[itrk]->dxy(pvtx->position()),eventWeight ) ; 
						histProducer.FillHisto_floatFromMap("jet_ipTagTrack_dz", jetflavour, bjetFromGSplit, selTrks_svTagInfos_jets[itrk]->dz(pvtx->position()),eventWeight) ; 
						histProducer.FillHisto_floatFromMap("jet_ipTagTrack_zIP", jetflavour, bjetFromGSplit, selTrks_svTagInfos_jets[itrk]->dz(pvtx->position()) - pvtx->z(),eventWeight) ; 
					}
				}
				else {
					edm::LogWarning("impactParameterTagInfos") << "Not present" ; 
				}
				if ( groomedJetMatch->hasTagInfo("secondaryVertex") ) { 
					const reco::SecondaryVertexTagInfo *svTagInfo =  groomedJetMatch->tagInfoSecondaryVertex("secondaryVertex");
					histProducer.FillHisto_floatFromMap("jet_ipTagNSV", jetflavour, bjetFromGSplit, svTagInfo->nVertices(), eventWeight) ; 
				}

				//// Loop over all subjets 
				double subjetBdisc(1) ; 
				int subjetflavour(0) ; 
				bool bsubjetFromGSplit ; 
				for (std::vector<const pat::Jet*>::const_iterator itsubjet = subjets.begin(); itsubjet != subjets.end(); ++itsubjet) {

					double subjetPtUncorr = (*itsubjet)->correctedP4(0).pt() ;  
					h1_JetPt -> Fill( subjetPtUncorr, eventWeight) ; 
					h1_JetEta -> Fill( (*itsubjet)->correctedP4(0).eta(), eventWeight) ; 
					h1_JetMass -> Fill ((*itsubjet)->correctedP4(0).mag(), eventWeight) ; 
					h1_JetMassDrop -> Fill ((*itsubjet)->correctedP4(0).mag()/groomedJetMatch->correctedP4(0).mag(), eventWeight) ; 
					//h1_JetMassDropOrg -> Fill ((*itsubjet)->correctedP4(0).mag()/it->correctedP4(0).mag(), eventWeight) ; 
					h1_JetMassDropOrg -> Fill ((*itsubjet)->correctedP4(0).mag()/it->correctedJet("Uncorrected").mass(), eventWeight) ; 
					const reco::TrackRefVector subjetTracks = (*itsubjet)->associatedTracks() ; 
					int nsubjetTracks = subjetTracks.size() ; 
					h1_NTracks -> Fill(nsubjetTracks, eventWeight) ; 
					if ( subjetBdisc > (*itsubjet)->bDiscriminator("combinedSecondaryVertexBJetTags") ) subjetBdisc = (*itsubjet)->bDiscriminator("combinedSecondaryVertexBJetTags") ; 

					subjetflavour = abs((*itsubjet)->partonFlavour()) ;
					bsubjetFromGSplit = ( bjetFromGSplit && (subjetflavour == 5) ) ? true : false ; 

					if (subjetflavour == 5 && bsubjetFromGSplit == false) ++nsubjets_b ; 
					if (subjetflavour == 5 && bsubjetFromGSplit == true)  ++nsubjets_bfromg ; 
					if (subjetflavour == 4 && bsubjetFromGSplit == false) ++nsubjets_c ; 
					if ( (subjetflavour < 4 || subjetflavour == 21)  && bsubjetFromGSplit == false) ++nsubjets_l ; 

					histProducer.FillHisto_floatFromMap("subjet_pt_all",  subjetflavour, bsubjetFromGSplit, (*itsubjet)->correctedP4(0).pt(),  eventWeight) ; 
					histProducer.FillHisto_floatFromMap("subjet_eta_all", subjetflavour, bsubjetFromGSplit, (*itsubjet)->correctedP4(0).eta(), eventWeight) ; 
					histProducer.FillHisto_floatFromMap("subjet_phi_all", subjetflavour, bsubjetFromGSplit, (*itsubjet)->correctedP4(0).phi(), eventWeight) ; 
					histProducer.FillHisto_floatFromMap("subjet_mass_all",subjetflavour, bsubjetFromGSplit, (*itsubjet)->correctedP4(0).mag(), eventWeight) ; 
					histProducer.FillHisto_floatFromMap("subjet_massdrop_all",subjetflavour, bsubjetFromGSplit, (*itsubjet)->correctedP4(0).mag()/groomedJetMatch->correctedP4(0).mag(), eventWeight) ; 
					histProducer.FillHisto_floatFromMap("subjet_massdropOrg_all",subjetflavour, bsubjetFromGSplit, (*itsubjet)->correctedP4(0).mag()/it->correctedP4(0).mag(), eventWeight) ; 

					const reco::TrackRefVector trks_subjets = (*itsubjet)->associatedTracks() ; 
					histProducer.FillHisto_floatFromMap("subjet_track_multi",  subjetflavour, bsubjetFromGSplit, trks_subjets.size(),  eventWeight) ;

					for (reco::track_iterator ittrk_subjets = trks_subjets.begin(); ittrk_subjets != trks_subjets.end(); ++ittrk_subjets) {
						histProducer.FillHisto_floatFromMap("subjet_track_chi2",  subjetflavour, bsubjetFromGSplit, (*ittrk_subjets)->normalizedChi2(),                        eventWeight) ;
						histProducer.FillHisto_floatFromMap("subjet_track_charge",subjetflavour, bsubjetFromGSplit, (*ittrk_subjets)->charge(),                                eventWeight) ;
						histProducer.FillHisto_floatFromMap("subjet_track_nHit",  subjetflavour, bsubjetFromGSplit, (*ittrk_subjets)->numberOfValidHits(),                     eventWeight) ;
						histProducer.FillHisto_floatFromMap("subjet_track_HPix",  subjetflavour, bsubjetFromGSplit, ((*ittrk_subjets)->hitPattern()).numberOfValidPixelHits(), eventWeight) ;
					}

					//const edm::RefVector<reco::TrackCollection> selTrks_svTagInfos_subjets = ( (*itsubjet)->tagInfoTrackIP("secondaryVertex") )->selectedTracks() ; 
					//int nSelTrks_subjets = selTrks_svTagInfos_subjets.size() ; 

					//histProducer.FillHisto_floatFromMap("subjet_track_multi_sel",  subjetflavour, bsubjetFromGSplit, nSelTrks_subjets,  eventWeight) ;

				}
				h1_csvBDisc->Fill(subjetBdisc,eventWeight) ; 
				histProducer.FillHisto_floatFromMap("subjet_csvBDisc",  subjetflavour, bsubjetFromGSplit, subjetBdisc, eventWeight) ; 


			} //// Groomed jet with two subjets found 

		} //// PF jet matched to groomed and packed jet 
		else  {
			edm::LogError("NoMatchingGroomedJet") << "Matching groomed jet not found. Using the original jet mass.";
		} //// PF jet not matched to any groomed and packed jet 

	} //// jet selection 

	h1_All_NJets    -> Fill(njets, eventWeight) ; 
	h1_All_NJetsMatched -> Fill(njetsMatched, eventWeight) ; 
	h1_All_NJetsBtagged -> Fill(njetsBtagged, eventWeight) ; 

	histProducer.FillHisto_floatFromMap("jet_multi",5  ,false ,njets_b      ,eventWeight) ;  
	histProducer.FillHisto_floatFromMap("jet_multi",5  ,true  ,njets_bfromg ,eventWeight) ;  
	histProducer.FillHisto_floatFromMap("jet_multi",4  ,false ,njets_c      ,eventWeight) ;  
	histProducer.FillHisto_floatFromMap("jet_multi",21 ,false ,njets_l      ,eventWeight) ;  

	histProducer.FillHisto_floatFromMap("subjet_multi",5  ,false ,nsubjets_b      ,eventWeight) ;  
	histProducer.FillHisto_floatFromMap("subjet_multi",5  ,true  ,nsubjets_bfromg ,eventWeight) ;  
	histProducer.FillHisto_floatFromMap("subjet_multi",4  ,false ,nsubjets_c      ,eventWeight) ;  
	histProducer.FillHisto_floatFromMap("subjet_multi",21 ,false ,nsubjets_l      ,eventWeight) ;  

#ifdef THIS_IS_AN_EVENT_EXAMPLE
	Handle<ExampleData> pIn;
	iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	ESHandle<SetupData> pSetup;
	iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
	void 
HiggsTagValidation::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
HiggsTagValidation::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
	void 
HiggsTagValidation::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
	void 
HiggsTagValidation::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
	void 
HiggsTagValidation::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
	void 
HiggsTagValidation::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HiggsTagValidation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HiggsTagValidation);
