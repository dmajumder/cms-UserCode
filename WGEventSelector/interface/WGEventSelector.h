#ifndef WG_EVENT_SELECTOR
#define WG_EVENT_SELECTOR 

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h" 

#include "DataFormats/BeamSpot/interface/BeamSpot.h" 

#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h" 
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include <string>
#include <fstream>

#include "TFile.h"
#include "TTree.h"

//
// constants, enums and typedefs
//
typedef std::vector<std::string> vstring;

//
// class decleration
//

class WGEventSelector : public edm::EDAnalyzer {
public:
  explicit WGEventSelector(const edm::ParameterSet&);
  ~WGEventSelector();


private:
  virtual void beginJob() ;

  void buildTree() ;
  void clearTreeVectors() ;

  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  void genInfo(const edm::Event&)     ; 
  bool trigInfo(const edm::Event&)    ;
  void beamSpotInfo(const edm::Event&);
  void muonInfo(const edm::Event&)    ;
  void photInfo(const edm::Event&)    ;
  void elecInfo(const edm::Event&)    ;
  void jetInfo(const edm::Event&)     ;
  void metInfo(const edm::Event&)     ;

  virtual void endJob() ;

  // ----------member data ---------------------------
  
  //---- configurable parameters --------
  
  bool isMC_;

  edm::InputTag genParticleTag_ ; 
  edm::InputTag l1gtrrTag_ ; 
  edm::InputTag l1ObjectMapTag_ ; 
  edm::InputTag triggerResultsTag_ ; 

  edm::InputTag genTag_      ;
  edm::InputTag beamSpotTag_ ;
  edm::InputTag vertexTag_   ; 
  edm::InputTag muoTag_      ;
  edm::InputTag phoTag_      ;
  edm::InputTag trackTag_    ;
  edm::InputTag elecTag_     ;
  edm::InputTag jetTag_      ; 
  edm::InputTag metTag_      ; 

  double pvTrkMinWt_         ; 
  double maxPtClosestTrack_  ; 
  double muPtMin_            ;
  double muEtaMax_           ; 
  double phPtMin_            ;
  double phEtaMax_           ; 

  ofstream outfile_ ; 
  std::string outfilename_ ; 
  TTree *TPat_ ;

  edm::Handle<std::vector<reco::GenParticle> > genHandle_ ; 
  edm::Handle<edm::TriggerResults> trigRes_ ;
  edm::Handle<L1GlobalTriggerReadoutRecord> l1gtrr_ ; 
  edm::Handle<L1GlobalTriggerObjectMapRecord> gtObjectMapRecord_ ; 
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle_ ;
  edm::Handle<reco::VertexCollection> recVtxs_ ;
  edm::Handle<edm::View<pat::Muon> > muonHandle_ ;
  edm::Handle<edm::View<pat::Photon> > photHandle_ ;
  edm::Handle<edm::View<reco::Track> > tracks_ ;
  edm::Handle<edm::View<pat::Electron> > electronHandle_ ;
  std::vector<edm::Handle<edm::View<pat::Jet> > > jetHandle_vect_ ;
  edm::Handle<edm::View<pat::Jet> > jetHandle_ ;
  edm::Handle<edm::View<pat::MET> > metHandle_ ;
  std::vector<math::XYZPoint> *vtx_ ; 

   //---- TREE variables -------- 

  unsigned long nEvt_;
  unsigned long runNo_ ;
  unsigned long evtNo_ ;
  unsigned long lumi_  ;
  unsigned long bunch_ ;

  vstring *l1Names_ ; 
  vstring *hltNames_ ; 
  vstring *hltWasRun_ ; 
  vstring *l1Accept_ ; 
  vstring *hltAccept_ ;
  vstring *hltError_ ; 
  bool isL1Accept_ ; 
  bool isHLTAccept_ ;
  int nHLTwasRun_ ; 
  int nHLTAccept_ ; 
  int nHLTErrors_ ; 
  int nL1Accept_ ; 
      
  reco::BeamSpot bs_ ;

  int nGenMu_ ;
  int nGenNu_ ;
  int nGenGamma_ ;
  
  double  beamSpotX0_                 ;
  double  beamSpotY0_                 ;
  double  beamSpotZ0_                 ; 

  int n1ryVertices_                   ; 
  std::vector<double> *xVtx_          ; 
  std::vector<double> *yVtx_          ; 
  std::vector<double> *vtxXError_     ; 
  std::vector<double> *vtxYError_     ; 
  std::vector<double> *xVtxBeamSpot_  ; 
  std::vector<double> *yVtxBeamSpot_  ; 
  std::vector<double> *zVtx_          ; 
  std::vector<int   > *nTrkFromVtx_   ; 
  std::vector<int   > *nTrkVtxWithWt_ ;  
  std::vector<double> *vtxchi2_       ;
  std::vector<double> *vtxndof_       ;
  std::vector<double> *vtxNormChi2_   ;  
  std::vector<double> *vtxprob_       ; 
  
  std::vector<int   > *genMuCh_    ;
  std::vector<double> *genMuPt_    ;   
  std::vector<double> *genMuEta_   ;
  std::vector<double> *genMuPhi_   ;
  std::vector<double> *genMuE_     ;
  
  std::vector<int   > *genElCh_    ;
  std::vector<double> *genElPt_    ;   
  std::vector<double> *genElEta_   ;
  std::vector<double> *genElPhi_   ;
  std::vector<double> *genElE_     ; 
  
  std::vector<double> *genNuPt_    ;   
  std::vector<double> *genNuEta_   ;
  std::vector<double> *genNuPhi_   ;
  std::vector<double> *genNuE_     ;
  
  std::vector<double> *genGammaPt_ ;
  std::vector<double> *genGammaEta_;
  std::vector<double> *genGammaPhi_;
  std::vector<double> *genGammaE_  ;

  std::vector<int   > *genWCh_ ;
  std::vector<double> *genWPt_ ;
  std::vector<double> *genWEta_;
  std::vector<double> *genWPhi_;
  std::vector<double> *genWE_  ;

  std::vector<double> *genWGammaPt_ ;
  std::vector<double> *genWGammaEta_;
  std::vector<double> *genWGammaPhi_;
  std::vector<double> *genWGammaE_  ;

  int                         nRecoMuons_            ;
  std::vector<bool  >        *isInnerTrackNonNull_   ;   
  std::vector<bool  >        *isGlobalTrackNonNull_  ;   
  std::vector<int   >        *noOfValidHits_         ;
  std::vector<int   >        *noOfLostHits_          ;
  std::vector<double>        *recoMuonPt_            ;
  std::vector<int   >        *genMuonMother_         ;
  std::vector<int   >        *genMuonGrandmother_    ;
  std::vector<double>        *genMuonPt_             ;
  std::vector<double>        *recoMuonPtError_       ;
  std::vector<double>        *recoMuonPx_            ;
  std::vector<double>        *recoMuonPy_            ;
  std::vector<double>        *recoMuonPz_            ;
  std::vector<double>        *recoMuonEta_           ;
  std::vector<double>        *genMuonEta_            ;
  std::vector<double>        *recoMuonEtaError_      ;
  std::vector<double>        *recoMuonPhi_           ;
  std::vector<double>        *genMuonPhi_            ;
  std::vector<double>        *recoMuonPhiError_      ;
  std::vector<double>        *recoMuonP_             ;
  std::vector<double>        *recoMuonCharge_        ;
  std::vector<double>        *recoMuonEnergy_        ;
  std::vector<double>        *recoMuonIsoTrack_      ;
  std::vector<double>        *recoMuonIsoEcal_       ;
  std::vector<double>        *recoMuonIsoHcal_       ;
  std::vector<double>        *recoMuondxy_           ;
  std::vector<double>        *recoMuondxyError_      ;
  std::vector<double>        *recoMuondxySigma_      ;
  std::vector<double>        *recoMuondz_            ;
  std::vector<double>        *recoMuondzError_       ;
  std::vector<double>        *recoMuonNormalizedChi2_;

  int                  nRecoPhotons_                                        ;
  std::vector<int   > *genPhotonMother_     /* 0;:e,1:mu,2:tau,3:pi0,4:W */ ;  
  std::vector<int   > *genPhotonGrandmother_/* 0;:e,1:mu,2:tau,3:pi0,4:W */ ;  
  std::vector<bool  > *recoPhotonHasPixelSeed_                              ;
  std::vector<bool  > *isEBPhoton_                                          ;
  std::vector<bool  > *isEEPhoton_                                          ;
  std::vector<bool  > *isEBGapPhoton_                                       ;
  std::vector<bool  > *isEEGapPhoton_                                       ;
  std::vector<bool  > *isEBEEGapPhoton_                                     ;

  std::vector<double> *recoPhotonHoverE_                                    ;
  std::vector<double> *recoPhotonHDepth1overE_                              ;
  std::vector<double> *recoPhotonHDepth2overE_                              ; 
  std::vector<double> *recoPhotonR9_                                        ;
  std::vector<double> *recoPhotonSigmaEtaEta_                               ; 
  std::vector<double> *recoPhotonSigmaIetaIeta_                             ;
  std::vector<double> *recoPhotonR1x5_                                      ;
  std::vector<double> *recoPhotonR2x5_                                      ;
  std::vector<double> *recoPhotonE1x5_                                      ;
  std::vector<double> *recoPhotonE2x5_                                      ;
  std::vector<double> *recoPhotonE3x3_                                      ;
  std::vector<double> *recoPhotonE5x5_                                      ;
  std::vector<double> *recoPhotonMaxEXtal_                                  ; 

  std::vector<double> *recoPhotonTrkSumPtHollowConeDR04_                    ;
  std::vector<double> *recoPhotonEcalRecHitSumEtConeDR04_                   ;
  std::vector<double> *recoPhotonHcalTowerSumEtConeDR04_                    ;
  std::vector<int   > *nTrkSolidConeDR04_                                   ; 
  std::vector<int   > *nTrkHollowConeDR04_                                  ;

  std::vector<double> *recoPhotonTrkSumPtHollowConeDR03_                    ;
  std::vector<double> *recoPhotonEcalRecHitSumEtConeDR03_                   ;
  std::vector<double> *recoPhotonHcalTowerSumEtConeDR03_                    ;
  std::vector<int   > *nTrkSolidConeDR03_                                   ;
  std::vector<int   > *nTrkHollowConeDR03_                                  ;
  std::vector<double> *recoPhotonDelRClosestTrack_                          ; 

  std::vector<double> *recoPhotonPt_                                        ;
  std::vector<double> *genPhotonPt_                                         ;
  std::vector<double> *recoPhotonEta_                                       ;
  std::vector<double> *genPhotonEta_                                        ;
  std::vector<double> *recoPhotonPhi_                                       ;
  std::vector<double> *genPhotonPhi_                                        ;
  std::vector<double> *recoPhotonPx_                                        ;
  std::vector<double> *recoPhotonPy_                                        ;
  std::vector<double> *recoPhotonPz_                                        ;
  std::vector<double> *recoPhotonCharge_                                    ;
  std::vector<double> *recoPhotonEnergy_                                    ;

  std::vector<double> recoPhotonIsoTrack_                                   ;
  std::vector<double> recoPhotonIsoEcal_                                    ;
  std::vector<double> recoPhotonIsoHcal_                                    ;

  int                  nRecoElectrons       ;
  std::vector<double> *recoElectronPt       ; 
  std::vector<double> *recoElectronPx       ;
  std::vector<double> *recoElectronPy       ;
  std::vector<double> *recoElectronPz       ;
  std::vector<double> *recoElectronEta      ;
  std::vector<double> *recoElectronPhi      ;
  std::vector<double> *recoElectronP        ;
  std::vector<double> *recoElectronCharge   ;
  std::vector<double> *recoElectronEnergy   ;
  std::vector<double> *recoElectronIsoTrack ;
  std::vector<double> *recoElectronIsoEcal  ;
  std::vector<double> *recoElectronIsoHcal  ;

  int                 nJets_         ;
  std::vector<double> *jet_px_       ;
  std::vector<double> *jet_py_       ;
  std::vector<double> *jet_pz_       ;
  std::vector<double> *jet_pt_       ;
  std::vector<double> *jet_eta_      ;
  std::vector<double> *jet_phi_      ;
  std::vector<double> *jet_E_        ;
  std::vector<double> *jet_P_        ; 
  std::vector<double> *jet_charge_   ;

  std::vector<int   > *jet_nTracks_  ; 
  std::vector<bool  > *isCaloJet_    ; 
  std::vector<bool  > *isPFJet_      ; 
  std::vector<bool  > *isBasicJet_   ; 

  std::vector<double> *maxEInEMTowers_   ; 
  std::vector<double> *maxEInHadTowers_  ;
  std::vector<double> *hadEFraction_     ;
  std::vector<double> *emEFraction_      ; 
  std::vector<double> *hadEInHB_         ; 
  std::vector<double> *hadEInHE_         ; 
  std::vector<double> *hadEInHO_         ; 
  std::vector<double> *hadEInHF_         ; 
  std::vector<double> *emEnergyInEB_     ; 
  std::vector<double> *emEnergyInEE_     ; 
  std::vector<double> *emEnergyInHF_     ; 
  std::vector<double> *towersArea_       ; 
  std::vector<double> *jet_n60_          ;
  std::vector<double> *jet_n90_          ; 
  std::vector<double> *jet_fHPD_         ; 

  int                 *nKT4Jets         ; 
  std::vector<double> *jetKT4_pt        ;
  std::vector<double> *jetKT4_px        ;
  std::vector<double> *jetKT4_py        ;
  std::vector<double> *jetKT4_pz        ;
  std::vector<double> *jetKT4_eta       ;
  std::vector<double> *jetKT4_phi       ;
  std::vector<double> *jetKT4_charge    ;

  int                 *nPFcJets         ;
  std::vector<double> *jetPFc_pt        ;
  std::vector<double> *jetPFc_px        ;
  std::vector<double> *jetPFc_py        ;
  std::vector<double> *jetPFc_pz        ;
  std::vector<double> *jetPFc_eta       ;
  std::vector<double> *jetPFc_phi       ;
  std::vector<double> *jetPFc_charge    ; 

  int                 *nPFrJets         ;  
  std::vector<double> *jetPFr_pt        ; 
  std::vector<double> *jetPFr_px        ; 
  std::vector<double> *jetPFr_py        ; 
  std::vector<double> *jetPFr_pz        ; 
  std::vector<double> *jetPFr_eta       ; 
  std::vector<double> *jetPFr_phi       ; 
  std::vector<double> *jetPFr_charge    ; 

  double              Met_x            ;  
  double              genMet_x         ; 
  double              Met_y            ; 
  double              genMet_y         ; 
  double              sumet        ;    
  double              gensumet     ; 
  double              Met              ; 
  double              genMet           ; 
  double              Met_phi          ; 
  double              genMet_phi       ; 

  double              uncorAllMet_pt   ; 
  double              uncorAllMet_phi  ; 
  double              uncorJesMet_pt   ; 
  double              uncorJesMet_phi  ; 
  double              uncorMuonMet_pt  ; 
  double              uncorMuonMet_phi ; 

};

#endif 

