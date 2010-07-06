// -*- C++ -*-
//
// Package:    WGEventSelector
// Class:      WGEventSelector
// 
/**\class WGEventSelector WGEventSelector.cc WGammaAnalyzer/WGEventSelector/src/WGEventSelector.cc

 Description: [one line class summary]

 Select events containing muons and photons and missing transverse energy 

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Devdatta Majumder,22 1-031,+41227679681,
//         Created:  Sun Jan  3 08:11:01 CET 2010
// $Id$
//
//


// user include files
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CLHEP/Units/PhysicalConstants.h" 
#include "CLHEP/Vector/ThreeVector.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/TriggerNames.h" 

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h" 

#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"

#include "WGammaAnalyzer/WGEventSelector/interface/WGEventSelector.h" 

using namespace edm ;
using namespace reco ;
using namespace std ; 

//
// constructors and destructor
//
WGEventSelector::WGEventSelector(const edm::ParameterSet& iConfig) : 

  isMC_             (iConfig.getParameter<bool>("isMC"))                      ,
  genParticleTag_   (iConfig.getParameter<edm::InputTag>("genParticleTag"))   ,
  l1gtrrTag_        (iConfig.getParameter<edm::InputTag>("l1gtrrTag"))        ,
  l1ObjectMapTag_   (iConfig.getParameter<edm::InputTag>("l1ObjectMapTag"))   ,
  triggerResultsTag_(iConfig.getParameter<edm::InputTag>("triggerResultsTag")),
  genTag_           (iConfig.getParameter<edm::InputTag>("genTag" ))          ,  
  beamSpotTag_      (iConfig.getParameter<edm::InputTag>("beamSpotTag"))      , 
  vertexTag_        (iConfig.getParameter<edm::InputTag>("vertexTag"))        , 
  muoTag_           (iConfig.getParameter<edm::InputTag>("muonTag"))          , 
  phoTag_           (iConfig.getParameter<edm::InputTag>("phoTag"))           , 
  trackTag_         (iConfig.getParameter<edm::InputTag>("trackTag"))         , 
  elecTag_          (iConfig.getParameter<edm::InputTag>("electronTag"))      ,  
  jetTag_           (iConfig.getParameter<edm::InputTag>("jetTag"))           ,  
  metTag_           (iConfig.getParameter<edm::InputTag>("metTag"))           , 
  pvTrkMinWt_       (iConfig.getParameter<double>("pvTrkMinWt"))              ,
  maxPtClosestTrack_(iConfig.getParameter<double>("maxPtClosestTrack")), 
  muPtMin_          (iConfig.getParameter<double>("muPtMin")),   
  muEtaMax_         (iConfig.getParameter<double>("muEtaMax")), 
  phPtMin_          (iConfig.getParameter<double>("phPtMin")), 
  phEtaMax_         (iConfig.getParameter<double>("phEtaMax")), 
  outfilename_      (iConfig.getUntrackedParameter<string>("outfilename")),
  nEvt_(0),
  runNo_(0) ,
  evtNo_(0) ,
  lumi_(0) ,
  bunch_(0) , 
  l1Names_(0),
  hltNames_(0),
  hltWasRun_(0),
  l1Accept_(0), 
  hltAccept_(0),
  hltError_(0),
  isL1Accept_(false), 
  isHLTAccept_(false), 
  nHLTwasRun_(0),
  nHLTAccept_(0),
  nHLTErrors_(0),
  nL1Accept_(0)

{
   //now do what ever initialization is needed

}


WGEventSelector::~WGEventSelector() {
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called to for each event  ------------  

void WGEventSelector::clearTreeVectors() {

  l1Names_->clear() ;
  hltNames_->clear() ;  
  hltWasRun_->clear() ;
  l1Accept_->clear() ; 
  hltAccept_->clear() ;
  hltError_->clear() ;
  isL1Accept_ = false ;
  isHLTAccept_ = false ; 

  vtx_          ->clear() ; 
  xVtx_         ->clear() ; 	
  yVtx_         ->clear() ; 	
  vtxXError_    ->clear() ;
  vtxYError_    ->clear() ; 
  xVtxBeamSpot_ ->clear() ; 	
  yVtxBeamSpot_ ->clear() ; 	
  zVtx_         ->clear() ; 	
  nTrkFromVtx_  ->clear() ; 
  nTrkVtxWithWt_->clear() ; 
  vtxchi2_      ->clear() ;    
  vtxNormChi2_  ->clear() ; 
  vtxndof_      ->clear() ;
  vtxprob_      ->clear() ; 

  genWCh_    ->clear() ;                           
  genWPt_    ->clear() ; 
  genWEta_   ->clear() ; 
  genWPhi_   ->clear() ; 
  genWE_     ->clear() ; 
  
  genMuCh_    ->clear() ;                           
  genMuPt_    ->clear() ; 
  genMuEta_   ->clear() ; 
  genMuPhi_   ->clear() ; 
  genMuE_     ->clear() ; 
  
  genElCh_    ->clear() ; 
  genElPt_    ->clear() ; 
  genElEta_   ->clear() ; 
  genElPhi_   ->clear() ; 
  genElE_     ->clear() ; 
  
  genNuPt_    ->clear() ; 
  genNuEta_   ->clear() ; 
  genNuPhi_   ->clear() ; 
  genNuE_     ->clear() ;  
  
  genGammaPt_ ->clear() ; 
  genGammaEta_->clear() ; 
  genGammaPhi_->clear() ; 
  genGammaE_  ->clear() ; 

  genWGammaPt_ ->clear() ; 
  genWGammaEta_->clear() ; 
  genWGammaPhi_->clear() ; 
  genWGammaE_  ->clear() ; 

  isInnerTrackNonNull_   ->clear() ; 
  isGlobalTrackNonNull_  ->clear() ; 
  noOfValidHits_         ->clear() ;          
  noOfLostHits_          ->clear() ;          
  recoMuonPt_            ->clear() ;          
  genMuonMother_         ->clear() ;          
  genMuonGrandmother_    ->clear() ; 
  genMuonPt_             ->clear() ;          
  recoMuonPtError_       ->clear() ;
  recoMuonPx_            ->clear() ;          
  recoMuonPy_            ->clear() ;          
  recoMuonPz_            ->clear() ;          
  recoMuonEta_           ->clear() ;          
  genMuonEta_            ->clear() ;          
  recoMuonEtaError_      ->clear() ; 
  recoMuonPhi_           ->clear() ;          
  genMuonPhi_            ->clear() ;          
  recoMuonPhiError_      ->clear() ; 
  recoMuonP_             ->clear() ;          
  recoMuonCharge_        ->clear() ;          
  recoMuonEnergy_        ->clear() ;          
  recoMuonIsoTrack_      ->clear() ; 
  recoMuonIsoEcal_       ->clear() ;  
  recoMuonIsoHcal_       ->clear() ;  
  recoMuondxy_           ->clear() ;  
  recoMuondxyError_      ->clear() ;  
  recoMuondxySigma_      ->clear() ;  
  recoMuondz_            ->clear() ;  
  recoMuondzError_       ->clear() ;  
  recoMuonNormalizedChi2_->clear() ;  

  genPhotonMother_       ->clear() ;  
  genPhotonGrandmother_  ->clear() ;  
  recoPhotonHasPixelSeed_->clear() ;  
  isEBPhoton_            ->clear() ;  
  isEEPhoton_            ->clear() ;  
  isEBGapPhoton_         ->clear() ;  
  isEEGapPhoton_         ->clear() ;  
  isEBEEGapPhoton_       ->clear() ;  
                                   
  recoPhotonHoverE_       ->clear() ;             
  recoPhotonHDepth1overE_ ->clear() ;             
  recoPhotonHDepth2overE_ ->clear() ;             
  recoPhotonR9_           ->clear() ;  
  recoPhotonSigmaEtaEta_  ->clear() ;  
  recoPhotonSigmaIetaIeta_->clear() ;  
  recoPhotonR1x5_         ->clear() ;  
  recoPhotonR2x5_         ->clear() ;  
  recoPhotonE1x5_         ->clear() ;  
  recoPhotonE2x5_         ->clear() ;  
  recoPhotonE3x3_         ->clear() ;  
  recoPhotonE5x5_         ->clear() ;  

  recoPhotonTrkSumPtHollowConeDR04_ ->clear() ;  
  recoPhotonEcalRecHitSumEtConeDR04_->clear() ;  
  recoPhotonHcalTowerSumEtConeDR04_ ->clear() ;  
  nTrkSolidConeDR04_                ->clear() ;  
  nTrkHollowConeDR04_               ->clear() ;  
                                   
  recoPhotonTrkSumPtHollowConeDR03_ ->clear() ;  
  recoPhotonEcalRecHitSumEtConeDR03_->clear() ;  
  recoPhotonHcalTowerSumEtConeDR03_ ->clear() ;  
  nTrkSolidConeDR03_                ->clear() ;  
  nTrkHollowConeDR03_               ->clear() ; 
  recoPhotonDelRClosestTrack_       ->clear() ;
                                    
  recoPhotonPt_    ->clear() ;                   
  genPhotonPt_     ->clear() ;                   
  recoPhotonEta_   ->clear() ;                   
  genPhotonEta_    ->clear() ;                   
  recoPhotonPhi_   ->clear() ;                   
  genPhotonPhi_    ->clear() ;                   
  recoPhotonPx_    ->clear() ;                   
  recoPhotonPy_    ->clear() ;                   
  recoPhotonPz_    ->clear() ;                   
  recoPhotonCharge_->clear() ;                   
  recoPhotonEnergy_->clear() ;                   

  recoElectronPt      ->clear() ;  
  recoElectronPx      ->clear() ;  
  recoElectronPy      ->clear() ;  
  recoElectronPz      ->clear() ;  
  recoElectronEta     ->clear() ;  
  recoElectronPhi     ->clear() ;  
  recoElectronP       ->clear() ;  
  recoElectronCharge  ->clear() ;  
  recoElectronEnergy  ->clear() ;  
  recoElectronIsoTrack->clear() ;  
  recoElectronIsoEcal ->clear() ;  
  recoElectronIsoHcal ->clear() ;  

  jet_px_    ->clear() ;  
  jet_py_    ->clear() ;  
  jet_pz_    ->clear() ;  
  jet_pt_    ->clear() ;  
  jet_eta_   ->clear() ;  
  jet_phi_   ->clear() ;  
  jet_E_     ->clear() ; 
  jet_P_     ->clear() ; 
  jet_charge_->clear() ;  

  maxEInEMTowers_  ->clear() ; 
  maxEInHadTowers_ ->clear() ; 
  hadEFraction_    ->clear() ; 
  emEFraction_     ->clear() ; 
  hadEInHB_        ->clear() ; 
  hadEInHE_        ->clear() ; 
  hadEInHO_        ->clear() ; 
  hadEInHF_        ->clear() ; 
  emEnergyInEB_    ->clear() ; 
  emEnergyInEE_    ->clear() ; 
  emEnergyInHF_    ->clear() ; 
  towersArea_      ->clear() ; 
  jet_n60_         ->clear() ; 
  jet_n90_         ->clear() ; 
  jet_fHPD_        ->clear() ; 

  jetKT4_px    ->clear() ;  
  jetKT4_py    ->clear() ;  
  jetKT4_pz    ->clear() ;  
  jetKT4_eta   ->clear() ;  
  jetKT4_phi   ->clear() ;  
  jetKT4_charge->clear() ;  

  jetPFc_pt    ->clear() ;  
  jetPFc_px    ->clear() ;  
  jetPFc_py    ->clear() ;  
  jetPFc_pz    ->clear() ;  
  jetPFc_eta   ->clear() ;  
  jetPFc_phi   ->clear() ;  
  jetPFc_charge->clear() ;  

  jetPFr_pt    ->clear() ;  
  jetPFr_px    ->clear() ;  
  jetPFr_py    ->clear() ;  
  jetPFr_pz    ->clear() ;  
  jetPFr_eta   ->clear() ;  
  jetPFr_phi   ->clear() ;  
  jetPFr_charge->clear() ;  

} // end clearTreeVectors 

//get gen Info 
void WGEventSelector::genInfo(const edm::Event &iEvent) {

  iEvent.getByLabel(genParticleTag_,genHandle_); 

  // find w and photon
  const GenParticle * photon = 0;
  const GenParticle * wBoson = 0;
 
  for( size_t ii=0; ii<genHandle_->size(); ii++ ) {
    const GenParticle * particle = &(*genHandle_)[ii];
    if( particle->status() == 3 ) {
      if( particle->pdgId() == 22 ) {
	if( photon == 0 || particle->pt() > photon->pt() ) { photon = particle; } 
      } else if( abs(particle->pdgId()) == 24 ) {
	if( wBoson == 0 || particle->pt() > wBoson->pt() ) { wBoson = particle; }  
      }
    }
  }

  if ( photon != 0 ) {

    genGammaPt_->push_back(photon->pt()) ;
    genGammaEta_->push_back(photon->eta()) ;
    genGammaPhi_->push_back(photon->phi()) ;
    genGammaE_->push_back(photon->energy()) ; 

  } 
    
  if( wBoson != 0 ) { 

    // find w decay products
    const Candidate * muon    = 0;
    const Candidate *electron = 0;  
    const Candidate * nu     = 0;
    for( GenParticle::const_iterator it = wBoson->begin(); it!=wBoson->end(); it++ ) { 
      if( abs(it->pdgId()) == 11 || abs(it->pdgId()) == 13 || abs(it->pdgId()) == 15 ) {
        if ( abs(it->pdgId()) == 11 ) muon = &(*it) ; 
        if ( abs(it->pdgId()) == 13 ) electron = &(*it) ; 
      } else if ( abs(it->pdgId()) == 12 || abs(it->pdgId()) == 14 || abs(it->pdgId()) == 16 ) {
        nu = &(*it);	    
      }
    }

    genWCh_->push_back(wBoson->charge()) ; 
    genWPt_->push_back(wBoson->pt()) ;
    genWEta_->push_back(wBoson->eta()) ;
    genWPhi_->push_back(wBoson->phi()) ;
    genWE_->push_back(wBoson->energy()) ; 
     	
    if (muon != 0 ) {
      genMuCh_->push_back(muon->charge()) ; 
      genMuPt_->push_back(muon->pt()) ;
      genMuEta_->push_back(muon->eta()) ;
      genMuPhi_->push_back(muon->phi()) ;
      genMuE_->push_back(muon->energy()) ; 
    }
     	
    if ( electron != 0 ) {
      genElCh_->push_back(electron->charge()) ;
      genElPt_->push_back(electron->pt()) ;
      genElEta_->push_back(electron->eta()) ;
      genElPhi_->push_back(electron->phi()) ;
      genElE_->push_back(electron->energy()) ; 
    }
     	
    if ( nu != 0 ) {
      genNuPt_->push_back(nu->pt()) ;
      genNuEta_->push_back(nu->eta()) ;
      genNuPhi_->push_back(nu->phi()) ;
      genNuE_->push_back(nu->energy()) ; 
    }        

    if ( photon != 0 ) { 
    
    math::XYZTLorentzVector p4_wgamma = photon->p4() + wBoson->p4();
    genWGammaPt_->push_back(p4_wgamma.pt()) ; 
    genWGammaEta_->push_back(p4_wgamma.eta()) ; 
    genWGammaPhi_->push_back(p4_wgamma.phi()) ; 
    genWGammaE_->push_back(p4_wgamma.energy()) ; 

    } 

  } 

  return ; 

} // get gen Info 

//get trig info 
bool WGEventSelector::trigInfo(const edm::Event &iEvent) {

 iEvent.getByLabel(triggerResultsTag_, trigRes_) ; 
 if ( trigRes_.failedToGet() )  return false ;
 if (trigRes_.isValid()) {
    if (trigRes_->wasrun()) ++nHLTwasRun_ ; 
    if (trigRes_->error() ) nHLTErrors_++; 
    isHLTAccept_= trigRes_->accept() ;
    if (isHLTAccept_) { 
            
      ++nHLTAccept_ ; 
      cout << " HLTAccept " << isHLTAccept_ << endl ; 

      edm::TriggerNames triggerNames;
      triggerNames.init(*trigRes_) ;
      *hltNames_=triggerNames.triggerNames(); 
      const unsigned int size(hltNames_->size()); 

      for (unsigned int ii=0; ii< size; ++ii) {
        if (trigRes_->wasrun(ii)) hltWasRun_->push_back(hltNames_->at(ii)) ; 
        if (trigRes_->accept(ii)) hltAccept_->push_back(hltNames_->at(ii)) ;
        if (trigRes_->error(ii)) hltError_->push_back(hltNames_->at(ii)) ; 
      }

      // get L1 Trigger Info
      iEvent.getByLabel(l1gtrrTag_,l1gtrr_) ;
      if ( !l1gtrr_.failedToGet() ) {

        iEvent.getByLabel(l1ObjectMapTag_, gtObjectMapRecord_);
        if ( !gtObjectMapRecord_.failedToGet() ) { 
           isL1Accept_ = l1gtrr_->decision() ; 
           if(isL1Accept_) { 
              ++nL1Accept_ ;
           }
           DecisionWord l1GtDecisionWord = l1gtrr_->decisionWord() ;
           const std::vector<L1GlobalTriggerObjectMap>& objMapVec = 
      	     gtObjectMapRecord_->gtObjectMap();
           for(std::vector<L1GlobalTriggerObjectMap>::const_iterator itMap = 
      		     objMapVec.begin(); 
      		     itMap != objMapVec.end(); 
      		     ++itMap) { 
              std::string algoName = (*itMap).algoName() ; 
              int l1index = (*itMap).algoBitNumber() ;
              if (l1GtDecisionWord[l1index]) l1Accept_->push_back(algoName) ; 
           }
            	        
        }  
      }
         
      return true ; 

    } else return false ; 


 } else {
   ++nHLTErrors_ ; 	   
   return false ; 
 }

 return false ; 

}//get trig info 

//get beam spot 
void WGEventSelector::beamSpotInfo(const edm::Event &iEvent) {

  iEvent.getByLabel(beamSpotTag_, recoBeamSpotHandle_);
  bs_ = *recoBeamSpotHandle_; 

  beamSpotX0_ = bs_.x0() ; 
  beamSpotY0_ = bs_.y0() ; 
  beamSpotZ0_ = bs_.z0() ; 

  iEvent.getByLabel(vertexTag_, recVtxs_); 

  for(reco::VertexCollection::const_iterator v=recVtxs_->begin();v!=recVtxs_->end(); ++v) {
      if (!v->isValid()) continue ; 
      math::XYZPoint vtx(v->x(),v->y(),v->z()) ;
      vtx_->push_back(vtx) ; 
      xVtx_->push_back(v->position().x()) ; 
      yVtx_->push_back(v->position().y()) ; 
      vtxXError_->push_back(v->xError()) ;
      vtxYError_->push_back(v->yError()) ; 
      xVtxBeamSpot_->push_back(v->position().x() 
	      + (bs_.dxdz() * v->position().z() - bs_.z0())) ; 
      yVtxBeamSpot_->push_back(v->position().y() 
	      + (bs_.dydz() * v->position().z() - bs_.z0())) ; 
      zVtx_->push_back(v->position().z()) ; 
      nTrkFromVtx_->push_back(v->tracksSize()) ; 
      unsigned pvntrk(0) ;
      reco::Vertex::trackRef_iterator it_trk ; 
      for(it_trk=v->tracks_begin();it_trk!=v->tracks_end();++it_trk) { 
        if (v->trackWeight(*it_trk) > pvTrkMinWt_) pvntrk++ ; 
      }
      nTrkVtxWithWt_->push_back(pvntrk) ; 
      vtxchi2_->push_back(v->chi2()) ;
      vtxndof_->push_back(v->ndof()) ; 
      vtxNormChi2_->push_back(v->normalizedChi2()) ;
      vtxprob_->push_back(ChiSquaredProbability(v->chi2() ,v->ndof())) ;

      ++n1ryVertices_ ; 
  
   }      

   return ; 

} // end beamSpotInfo 

void WGEventSelector::muonInfo(const edm::Event& iEvent) {
  
  iEvent.getByLabel(muoTag_,muonHandle_);
  edm::View<pat::Muon> muons = *muonHandle_;

  //+++++++++++++++++++++++ for PAT muons +++++++++++++++++++++++ 

  nRecoMuons_=0;

  for(edm::View<pat::Muon>::const_iterator mu_iter = muons.begin(); mu_iter!=muons.end(); ++mu_iter){

    if ( !mu_iter->isGlobalMuon()  
         || mu_iter->pt()<muPtMin_
	 || abs(mu_iter->eta())>muEtaMax_ ) continue ; 
    
    TrackRef muonGlobalTrack = mu_iter->globalTrack() ; 
    /* Muon quality variables */ 
    if ( muonGlobalTrack.isNonnull() ) { 
      isGlobalTrackNonNull_->push_back(1) ; 
      noOfValidHits_->push_back(muonGlobalTrack->numberOfValidHits()) ; 
      noOfLostHits_->push_back(muonGlobalTrack->numberOfLostHits()) ; 
      double dxy_min(0.);
      double dz_min(0.) ; 
      std::vector<math::XYZPoint>::const_iterator vtx_iter ; 
      for (vtx_iter=vtx_->begin();vtx_iter!=vtx_->end();++vtx_iter) {
        double dxy = muonGlobalTrack->dxy(*vtx_iter) ; 
	if (dxy < dxy_min) dxy_min = dxy ;
	double dz = muonGlobalTrack->dz(*vtx_iter) ; 
	if (dz < dz_min) dz_min = dz ; 
      }
      recoMuondxy_->push_back(dxy_min) ; 
      recoMuondz_->push_back(dz_min) ; 
      recoMuondxyError_->push_back(muonGlobalTrack->dxyError()) ; 
      recoMuondxySigma_->push_back(
		      sqrt( muonGlobalTrack->dxyError()*muonGlobalTrack->dxyError() + 
			      bs_.BeamWidthX()*bs_.BeamWidthX() ) ); 
      recoMuondzError_->push_back(muonGlobalTrack->dzError() );
      recoMuonNormalizedChi2_->push_back(muonGlobalTrack->normalizedChi2() ); 
      recoMuonPtError_->push_back(muonGlobalTrack->ptError());
      recoMuonEtaError_->push_back(muonGlobalTrack->etaError());
      recoMuonPhiError_->push_back(muonGlobalTrack->phiError()); 
    } else {
      isGlobalTrackNonNull_->push_back(0) ;
    }
    
    /* Muon kinematic variables */	  
    recoMuonPt_->push_back(mu_iter->pt());
    recoMuonPx_->push_back(mu_iter->px());
    recoMuonPy_->push_back(mu_iter->py());
    recoMuonPz_->push_back(mu_iter->pz());
    recoMuonEta_->push_back(mu_iter->eta());
    recoMuonPhi_->push_back(mu_iter->phi());
    recoMuonCharge_->push_back(mu_iter->charge());
    recoMuonEnergy_->push_back(mu_iter->energy());
    
    /* Muon Isolation variables */ 
    recoMuonIsoTrack_->push_back(mu_iter->trackIso()) ; 
    recoMuonIsoEcal_->push_back(mu_iter->ecalIso()) ; 
    recoMuonIsoHcal_->push_back(mu_iter->hcalIso()) ; 

    if(isMC_) {
      if(mu_iter->genParticleRef().isNonnull() && mu_iter->genParticlesSize()!=0){

        genMuonPt_->push_back(mu_iter->genLepton()->pt()) ; 
        genMuonEta_->push_back(mu_iter->genLepton()->eta());
        genMuonPhi_->push_back(mu_iter->genLepton()->phi()) ; 
        for(size_t i = 0; i < mu_iter->genParticlesSize(); ++ i) { 
          mu_iter->genParticleRef()->status() ; 	       
        }
              
        genMuonGrandmother_->push_back(mu_iter->genLepton()->mother()->mother()->pdgId()) ; 
        genMuonMother_->push_back(mu_iter->genLepton()->mother()->pdgId()) ; 
      }
      else {
        genMuonPt_->push_back(-99.) ; 
        genMuonEta_->push_back(-99.) ;
        genMuonPhi_->push_back(-99.) ; 
        genMuonGrandmother_->push_back(0) ; 
        genMuonMother_->push_back(0) ; 
      }
    }

    nRecoMuons_++;
  }

  return ; 

} //end muonInfo

void WGEventSelector::photInfo(const edm::Event& iEvent) {


  iEvent.getByLabel(phoTag_,photHandle_);
  edm::View<pat::Photon> photons = *photHandle_;

  nRecoPhotons_=0;

  for(edm::View<pat::Photon>::const_iterator ph_iter = photons.begin(); ph_iter!= photons.end() ; ++ph_iter){

      if (ph_iter->pt()<phPtMin_ 
          || abs(ph_iter->eta()) > phEtaMax_) continue ; 
    
      /* Photon flags */ 
      isEBPhoton_->push_back(ph_iter->isEB()) ; 
      isEEPhoton_->push_back(ph_iter->isEE()) ; 
      isEBGapPhoton_->push_back(ph_iter->isEBGap()) ;
      isEEGapPhoton_->push_back(ph_iter->isEEGap()) ;
      isEBEEGapPhoton_->push_back(ph_iter->isEBEEGap()) ;
      if (isMC_) recoPhotonHasPixelSeed_->push_back(ph_iter->hasPixelSeed()) ;
      
      /* Photon kinematic variables */ 
      recoPhotonPt_->push_back(ph_iter->pt());
      recoPhotonEta_->push_back(ph_iter->eta());
      recoPhotonPhi_->push_back(ph_iter->phi());
      recoPhotonPx_->push_back(ph_iter->px());
      recoPhotonPy_->push_back(ph_iter->py());
      recoPhotonPz_->push_back(ph_iter->pz());
      recoPhotonCharge_->push_back(ph_iter->charge());
      recoPhotonEnergy_->push_back(ph_iter->energy());

      /* Photon shower shape variable */ 
      recoPhotonHoverE_->push_back(ph_iter->hadronicOverEm()) ; 
      recoPhotonHDepth1overE_->push_back(ph_iter->hadronicDepth1OverEm()) ; 
      recoPhotonHDepth2overE_->push_back(ph_iter->hadronicDepth2OverEm()) ; 
      recoPhotonR9_->push_back(ph_iter->r9()) ; 
      recoPhotonSigmaEtaEta_->push_back(ph_iter->sigmaEtaEta()) ; 
      recoPhotonSigmaIetaIeta_->push_back(ph_iter->sigmaIetaIeta()) ; 
      recoPhotonR1x5_->push_back(ph_iter->r1x5());
      recoPhotonR2x5_->push_back(ph_iter->r2x5());
      recoPhotonE1x5_->push_back(ph_iter->e1x5());
      recoPhotonE2x5_->push_back(ph_iter->e2x5());
      recoPhotonE3x3_->push_back(ph_iter->e3x3());
      recoPhotonE3x3_->push_back(ph_iter->e3x3());
      recoPhotonE5x5_->push_back(ph_iter->e5x5());      
      recoPhotonMaxEXtal_->push_back(ph_iter->maxEnergyXtal()); 
      
      /* PhotonID variables in 0.4 cone */      
      recoPhotonTrkSumPtHollowConeDR04_->push_back(ph_iter->trkSumPtHollowConeDR04()) ; 
      recoPhotonEcalRecHitSumEtConeDR04_->push_back(ph_iter->ecalRecHitSumEtConeDR04()) ; 
      recoPhotonHcalTowerSumEtConeDR04_->push_back(ph_iter->hcalTowerSumEtConeDR04()) ; 
      nTrkSolidConeDR04_->push_back(ph_iter->nTrkSolidConeDR04()) ; 
      nTrkHollowConeDR04_->push_back(ph_iter->nTrkHollowConeDR04()) ; 

      /* PhotonID variables in 0.3 cone */ 
      recoPhotonTrkSumPtHollowConeDR03_->push_back(ph_iter->trkSumPtHollowConeDR03()) ; 
      recoPhotonEcalRecHitSumEtConeDR03_->push_back(ph_iter->ecalRecHitSumEtConeDR03()) ; 
      recoPhotonHcalTowerSumEtConeDR03_->push_back(ph_iter->hcalTowerSumEtConeDR03()) ; 
      nTrkSolidConeDR03_->push_back(ph_iter->nTrkSolidConeDR03()) ; 
      nTrkHollowConeDR03_->push_back(ph_iter->nTrkHollowConeDR03()) ; 

      CLHEP::Hep3Vector ph3v(ph_iter->px(),ph_iter->py(),ph_iter->pz()) ;
      math ::XYZPoint photonCaloPosition = ph_iter->caloPosition() ; 

      /* To find DeltaR of photon from nearest track */
      iEvent.getByLabel(trackTag_,tracks_) ; 
      double delR_max(0.) ;
      for (edm::View<reco::Track>::const_iterator track_iter = tracks_->begin();
		                                  track_iter != tracks_->end(); 
						  ++track_iter) {
        CLHEP::Hep3Vector trk3v(track_iter->px(),track_iter->py(),track_iter->pz()); 
        double delR = deltaR(ph3v,trk3v) ; 
	if (trk3v.mag() > maxPtClosestTrack_ && delR > delR_max) {
	  delR_max = delR ; 
	}        
      }
      recoPhotonDelRClosestTrack_->push_back(delR_max) ; 
      
      /* Matching to MC truth */ 
      if (isMC_) {
        if(ph_iter->genParticleRef().isNonnull() && ph_iter->genParticlesSize()!=0){

          genPhotonPt_->push_back(ph_iter->genPhoton()->pt()) ; 
          genPhotonEta_->push_back(ph_iter->genPhoton()->eta()) ; 
          genPhotonPhi_->push_back(ph_iter->genPhoton()->phi()) ; 
          //if (ph_iter->genParticleRef()->motherRef()->status() == 3) {	      
            genPhotonMother_->push_back(ph_iter->genParticleRef()->motherRef()->pdgId()) ; 
          //}
          //if (ph_iter->genParticleRef()->motherRef()->motherRef()->pdgId() == 3) { 
            genPhotonGrandmother_->push_back(ph_iter->genParticleRef()->motherRef()->motherRef()->pdgId()) ;  
          //}
        }
        else {
          genPhotonPt_->push_back(-99.) ; 
          genPhotonEta_->push_back(-99.) ; 
          genPhotonPhi_->push_back(-99.) ;
          genPhotonMother_->push_back(0) ; 
          genPhotonGrandmother_->push_back(0) ; 
        }
      }
      	  
      nRecoPhotons_++;
    
  }

  return ; 

} // end photInfo
  
void WGEventSelector::elecInfo(const edm::Event& iEvent) {

  iEvent.getByLabel(elecTag_,electronHandle_);
  edm::View<pat::Electron> electrons = *electronHandle_;

  nRecoElectrons=0;
  
  for(edm::View<pat::Electron>::const_iterator el_iter = electrons.begin(); el_iter!=electrons.end() ; ++el_iter){
    
    recoElectronPt->push_back(el_iter->pt());
    recoElectronPx->push_back(el_iter->px());
    recoElectronPy->push_back(el_iter->py());
    recoElectronPz->push_back(el_iter->pz());
    recoElectronEta->push_back(el_iter->eta());
    recoElectronPhi->push_back(el_iter->phi());
    recoElectronCharge->push_back(el_iter->charge());
    recoElectronEnergy->push_back(el_iter->energy());
    recoElectronIsoTrack->push_back(el_iter->trackIso());
    recoElectronIsoEcal->push_back(el_iter->ecalIso());
    recoElectronIsoHcal->push_back(el_iter->hcalIso());
    nRecoElectrons++;
    
  }

  return ; 

} // end elecInfo

void WGEventSelector::jetInfo(const edm::Event& iEvent) {

  jetHandle_vect_.clear() ;

  iEvent.getByLabel(jetTag_,jetHandle_);
  jetHandle_vect_.push_back(jetHandle_) ;
  edm::View<pat::Jet> jets = *jetHandle_ ;
  
   nJets_=0;
  for(edm::View<pat::Jet>::const_iterator jet_iter = jets.begin(); jet_iter!=jets.end(); ++jet_iter){
    
    jet_px_->push_back(jet_iter->px());
    jet_py_->push_back(jet_iter->py());
    jet_pz_->push_back(jet_iter->pz());
    jet_pt_->push_back(jet_iter->pt());
    jet_eta_->push_back(jet_iter->eta());
    jet_phi_->push_back(jet_iter->phi());
    jet_E_->push_back(jet_iter->energy()) ; 
    jet_P_->push_back(jet_iter->p()) ; 
    jet_charge_->push_back(jet_iter->charge()); 

    maxEInEMTowers_  ->push_back(jet_iter->maxEInEmTowers()) ; 
    maxEInHadTowers_ ->push_back(jet_iter->maxEInHadTowers()) ;
    hadEFraction_    ->push_back(jet_iter->energyFractionHadronic()) ; 
    emEFraction_     ->push_back(jet_iter->emEnergyFraction()) ;
    hadEInHB_        ->push_back(jet_iter->hadEnergyInHB()) ; 
    hadEInHE_        ->push_back(jet_iter->hadEnergyInHE()) ; 
    hadEInHO_        ->push_back(jet_iter->hadEnergyInHO()) ; 
    hadEInHF_        ->push_back(jet_iter->hadEnergyInHF()) ; 
    emEnergyInEB_    ->push_back(jet_iter->emEnergyInEB()) ; 
    emEnergyInEE_    ->push_back(jet_iter->emEnergyInEE()) ; 
    emEnergyInHF_    ->push_back(jet_iter->emEnergyInHF()) ; 
    towersArea_      ->push_back(jet_iter->towersArea()) ; 
    jet_n60_         ->push_back(jet_iter->n90()) ;
    jet_n90_         ->push_back(jet_iter->n60()) ; 
    jet_fHPD_        ->push_back(jet_iter->jetID().fHPD) ; 

    ++nJets_; 
  }

  return ; 
    
} // end jetInfo

void WGEventSelector::metInfo(const edm::Event& iEvent) {

  iEvent.getByLabel(metTag_,metHandle_);
  edm::View<pat::MET> mets = *metHandle_;

  Met_x     = mets[0].corEx();
  Met_y     = mets[0].corEy();
  Met       = mets[0].pt();
  Met_phi   = mets[0].phi();
  sumet     = mets[0].corSumEt(); 

  // For getting gen level info
  const GenMET* genmet = mets[0].genMET() ; 
  genMet_x = genmet->px(); 
  genMet_y = genmet->py() ;
  genMet = genmet->pt() ; 
  genMet_phi = genmet->phi() ; 
  gensumet = genmet->sumEt() ; 
  
  // For undoing all corrections
  uncorAllMet_pt  = mets[0].uncorrectedPt();
  uncorAllMet_phi = mets[0].uncorrectedPhi();
  
  // For undoing JES correction only
  uncorJesMet_pt  = mets[0].uncorrectedPt(pat::MET::uncorrJES);             
  uncorJesMet_phi = mets[0].uncorrectedPhi(pat::MET::uncorrJES);    //------do-------
  
  // For undoing muon correction
  uncorMuonMet_pt  = mets[0].uncorrectedPt(pat::MET::uncorrMUON);             
  uncorMuonMet_phi = mets[0].uncorrectedPhi(pat::MET::uncorrMUON);

  return ; 
  
}// end metInfo

void WGEventSelector::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  nEvt_++ ;

  clearTreeVectors() ; 

  runNo_ = iEvent.id().run()        ;
  evtNo_ = iEvent.id().event()      ;
  lumi_  = iEvent.luminosityBlock() ;
  bunch_ = iEvent.bunchCrossing()   ; 

  //if ( trigInfo(iEvent) ) { 
  //  if (isMC_) genInfo(iEvent) ; 
  //  beamSpotInfo(iEvent) ;  
  //  if (n1ryVertices_ > 0) {
  //    muonInfo(iEvent) ;
  //    if (nRecoMuons_ > 0) {
  //      photInfo(iEvent) ; 
  //      if (nRecoPhotons_ > 0) {
  //        elecInfo(iEvent) ; 
  //        jetInfo(iEvent) ; 
  //        metInfo(iEvent) ; 

  //        TPat_->Fill() ;
  //      }
  //    }
  //  }
  //} 

    beamSpotInfo(iEvent) ;  
    muonInfo(iEvent) ;
    photInfo(iEvent) ; 
    elecInfo(iEvent) ; 
    jetInfo(iEvent) ; 
    //metInfo(iEvent) ; 

    TPat_->Fill() ;
    
  return ; 
       
}

// ------------ method called once each job just before starting event loop  ------------
void WGEventSelector::beginJob() { 

   edm::Service<TFileService> fs;
   TPat_ = fs->make<TTree>("patTree", "patTree");
   buildTree() ;

   outfile_.open(outfilename_.c_str()) ; 

}

void WGEventSelector::buildTree() {

  l1Names_ = new std::vector<std::string>(); 
  hltNames_ = new std::vector<std::string>(); 
  hltWasRun_ = new std::vector<std::string>();
  l1Accept_ = new std::vector<std::string>(); 
  hltAccept_ = new std::vector<std::string>();
  hltError_ = new std::vector<std::string>();

  vtx_          = new std::vector<math::XYZPoint>() ;

  xVtx_         = new std::vector<double>() ; 
  yVtx_         = new std::vector<double>() ;
  vtxXError_    = new std::vector<double>() ;
  vtxYError_    = new std::vector<double>() ;  
  xVtxBeamSpot_ = new std::vector<double>() ;
  yVtxBeamSpot_ = new std::vector<double>() ;
  zVtx_         = new std::vector<double>() ;
  nTrkFromVtx_  = new std::vector<int   >() ;
  nTrkVtxWithWt_= new std::vector<int   >() ;
  vtxchi2_      = new std::vector<double>() ;
  vtxndof_      = new std::vector<double>() ;
  vtxNormChi2_  = new std::vector<double>() ; 
  vtxprob_      = new std::vector<double>() ; 
  
  genWCh_ = new std::vector<int   >() ; 
  genWPt_ = new std::vector<double>() ; 
  genWEta_= new std::vector<double>() ; 
  genWPhi_= new std::vector<double>() ; 
  genWE_  = new std::vector<double>() ; 

  genMuCh_    = new std::vector<int   >() ; 
  genMuPt_    = new std::vector<double>() ; 
  genMuEta_   = new std::vector<double>() ; 
  genMuPhi_   = new std::vector<double>() ; 
  genMuE_     = new std::vector<double>() ; 
  
  genElCh_    = new std::vector<int   >() ;   
  genElPt_    = new std::vector<double>() ; 
  genElEta_   = new std::vector<double>() ; 
  genElPhi_   = new std::vector<double>() ; 
  genElE_     = new std::vector<double>() ; 
  
  genNuPt_    = new std::vector<double>() ; 
  genNuEta_   = new std::vector<double>() ; 
  genNuPhi_   = new std::vector<double>() ; 
  genNuE_     = new std::vector<double>() ; 
  
  genGammaPt_ = new std::vector<double>() ; 
  genGammaEta_= new std::vector<double>() ; 
  genGammaPhi_= new std::vector<double>() ; 
  genGammaE_  = new std::vector<double>() ; 

  genWGammaPt_ = new std::vector<double>() ; 
  genWGammaEta_= new std::vector<double>() ; 
  genWGammaPhi_= new std::vector<double>() ; 
  genWGammaE_  = new std::vector<double>() ; 

  isInnerTrackNonNull_    = new std::vector<bool  >() ; 
  isGlobalTrackNonNull_   = new std::vector<bool  >() ; 
  noOfValidHits_          = new std::vector<int   >() ;
  noOfLostHits_           = new std::vector<int   >() ;
  recoMuonPt_             = new std::vector<double>() ;
  genMuonMother_          = new std::vector<int   >() ;
  genMuonGrandmother_     = new std::vector<int   >() ;
  genMuonPt_              = new std::vector<double>() ;
  recoMuonPtError_        = new std::vector<double>() ;
  recoMuonPx_             = new std::vector<double>() ;
  recoMuonPy_             = new std::vector<double>() ;
  recoMuonPz_             = new std::vector<double>() ;
  recoMuonEta_            = new std::vector<double>() ;
  genMuonEta_             = new std::vector<double>() ;
  recoMuonEtaError_       = new std::vector<double>() ;
  recoMuonPhi_            = new std::vector<double>() ;
  genMuonPhi_             = new std::vector<double>() ;
  recoMuonPhiError_       = new std::vector<double>() ;
  recoMuonP_              = new std::vector<double>() ;
  recoMuonCharge_         = new std::vector<double>() ;
  recoMuonEnergy_         = new std::vector<double>() ;
  recoMuonIsoTrack_       = new std::vector<double>() ;  
  recoMuonIsoEcal_        = new std::vector<double>() ;
  recoMuonIsoHcal_        = new std::vector<double>() ;
  recoMuondxy_            = new std::vector<double>() ;
  recoMuondxyError_       = new std::vector<double>() ;
  recoMuondxySigma_       = new std::vector<double>() ;
  recoMuondz_             = new std::vector<double>() ;
  recoMuondzError_        = new std::vector<double>() ;
  recoMuonNormalizedChi2_ = new std::vector<double>() ;

  genPhotonMother_        = new std::vector<int >() ;   
  genPhotonGrandmother_   = new std::vector<int >() ;   
  recoPhotonHasPixelSeed_ = new std::vector<bool>() ;   
  isEBPhoton_             = new std::vector<bool>() ;   
  isEEPhoton_             = new std::vector<bool>() ;   
  isEBGapPhoton_          = new std::vector<bool>() ;   
  isEEGapPhoton_          = new std::vector<bool>() ;   
  isEBEEGapPhoton_        = new std::vector<bool>() ;   
                                     
  recoPhotonHoverE_        = new std::vector<double>() ;   
  recoPhotonHDepth1overE_  = new std::vector<double>() ;   
  recoPhotonHDepth2overE_  = new std::vector<double>() ;   
  recoPhotonR9_            = new std::vector<double>() ;   
  recoPhotonSigmaEtaEta_   = new std::vector<double>() ;   
  recoPhotonSigmaIetaIeta_ = new std::vector<double>() ;   
  recoPhotonR1x5_          = new std::vector<double>() ;   
  recoPhotonR2x5_          = new std::vector<double>() ;   
  recoPhotonE1x5_          = new std::vector<double>() ;   
  recoPhotonE2x5_          = new std::vector<double>() ;   
  recoPhotonE3x3_          = new std::vector<double>() ;   
  recoPhotonE5x5_          = new std::vector<double>() ;   
  recoPhotonMaxEXtal_      = new std::vector<double>() ;

  recoPhotonTrkSumPtHollowConeDR04_  = new std::vector<double>() ;   
  recoPhotonEcalRecHitSumEtConeDR04_ = new std::vector<double>() ;   
  recoPhotonHcalTowerSumEtConeDR04_  = new std::vector<double>() ;   
  nTrkSolidConeDR04_                 = new std::vector<int   >() ;   
  nTrkHollowConeDR04_                = new std::vector<int   >() ;   
                                     
  recoPhotonTrkSumPtHollowConeDR03_  = new std::vector<double>() ;   
  recoPhotonEcalRecHitSumEtConeDR03_ = new std::vector<double>() ;   
  recoPhotonHcalTowerSumEtConeDR03_  = new std::vector<double>() ;   
  nTrkSolidConeDR03_                 = new std::vector<int   >() ;   
  nTrkHollowConeDR03_                = new std::vector<int   >() ;   
  recoPhotonDelRClosestTrack_        = new std::vector<double>() ;
                                    
  recoPhotonPt_     = new std::vector<double>() ;   
  genPhotonPt_      = new std::vector<double>() ;   
  recoPhotonEta_    = new std::vector<double>() ;   
  genPhotonEta_     = new std::vector<double>() ;   
  recoPhotonPhi_    = new std::vector<double>() ;   
  genPhotonPhi_     = new std::vector<double>() ;   
  recoPhotonPx_     = new std::vector<double>() ;   
  recoPhotonPy_     = new std::vector<double>() ;   
  recoPhotonPz_     = new std::vector<double>() ;   
  recoPhotonCharge_ = new std::vector<double>() ;   
  recoPhotonEnergy_ = new std::vector<double>() ;   

  recoElectronPt       = new std::vector<double>() ;  
  recoElectronPx       = new std::vector<double>() ;  
  recoElectronPy       = new std::vector<double>() ;  
  recoElectronPz       = new std::vector<double>() ;  
  recoElectronEta      = new std::vector<double>() ;  
  recoElectronPhi      = new std::vector<double>() ;  
  recoElectronP        = new std::vector<double>() ;  
  recoElectronCharge   = new std::vector<double>() ;  
  recoElectronEnergy   = new std::vector<double>() ;  
  recoElectronIsoTrack = new std::vector<double>() ;  
  recoElectronIsoEcal  = new std::vector<double>() ;  
  recoElectronIsoHcal  = new std::vector<double>() ;  

  jet_pt_     = new std::vector<double>() ;   
  jet_px_     = new std::vector<double>() ;   
  jet_py_     = new std::vector<double>() ;   
  jet_pz_     = new std::vector<double>() ;   
  jet_eta_    = new std::vector<double>() ;   
  jet_phi_    = new std::vector<double>() ;   
  jet_E_      = new std::vector<double>() ;
  jet_P_      = new std::vector<double>() ;
  jet_charge_ = new std::vector<double>() ;   

  maxEInEMTowers_  = new std::vector<double>() ;
  maxEInHadTowers_ = new std::vector<double>() ;  
  hadEFraction_    = new std::vector<double>() ;
  emEFraction_     = new std::vector<double>() ; 
  hadEInHB_        = new std::vector<double>() ;
  hadEInHE_        = new std::vector<double>() ;
  hadEInHO_        = new std::vector<double>() ;
  hadEInHF_        = new std::vector<double>() ;
  emEnergyInEB_    = new std::vector<double>() ;
  emEnergyInEE_    = new std::vector<double>() ;
  emEnergyInHF_    = new std::vector<double>() ;
  towersArea_      = new std::vector<double>() ;
  jet_n60_         = new std::vector<double>() ;
  jet_n90_         = new std::vector<double>() ;
  jet_fHPD_        = new std::vector<double>() ; 

  jetKT4_pt     = new std::vector<double>() ;   
  jetKT4_px     = new std::vector<double>() ;   
  jetKT4_py     = new std::vector<double>() ;   
  jetKT4_pz     = new std::vector<double>() ;   
  jetKT4_eta    = new std::vector<double>() ;   
  jetKT4_phi    = new std::vector<double>() ;   
  jetKT4_charge = new std::vector<double>() ;   

  jetPFc_pt     = new std::vector<double>() ;   
  jetPFc_px     = new std::vector<double>() ;   
  jetPFc_py     = new std::vector<double>() ;   
  jetPFc_pz     = new std::vector<double>() ;   
  jetPFc_eta    = new std::vector<double>() ;   
  jetPFc_phi    = new std::vector<double>() ;   
  jetPFc_charge = new std::vector<double>() ;   

  jetPFr_pt     = new std::vector<double>() ;   
  jetPFr_px     = new std::vector<double>() ;   
  jetPFr_py     = new std::vector<double>() ;   
  jetPFr_pz     = new std::vector<double>() ;   
  jetPFr_eta    = new std::vector<double>() ;   
  jetPFr_phi    = new std::vector<double>() ;   
  jetPFr_charge = new std::vector<double>() ;   
  
  TPat_->Branch("runNo_", &runNo_, "runNo_/I") ;   
  TPat_->Branch("evtNo_", &evtNo_, "evtNo_/I") ;   
  TPat_->Branch("lumi_",  &lumi_,  "lumi_/I")  ;   
  TPat_->Branch("bunch_", &bunch_, "bunch_/I") ; 
  
  TPat_->Branch("nEvt_",        &nEvt_,           "nEvt_/I:runNo_/I:evtNo_/I:lumi_/I:bunch_/I") ; 
  TPat_->Branch("l1Names_",     "vector<string>", &l1Names_);
  TPat_->Branch("hltNames_",    "vector<string>", &hltNames_) ; 
  TPat_->Branch("hltWasRun_",   "vector<string>", &hltWasRun_) ;
  TPat_->Branch("l1Accept_",    "vector<string>", &l1Accept_) ;
  TPat_->Branch("hltAccept_",   "vector<string>", &hltAccept_) ;
  TPat_->Branch("hltError_",    "vector<string>", &hltError_) ; 
  TPat_->Branch("isHLTAccept_", &isHLTAccept_,    "isHLTAccept_/O") ; 
  TPat_->Branch("isL1Accept_",  &isL1Accept_,     "isL1Accept_/O") ; 
  TPat_->Branch("nHLTwasRun_",  &nHLTwasRun_,     "nHLTwasRun_/I:nHLTAccept_/I:nHLTErrors_/I") ; 

  TPat_->Branch("beamSpotX0_", &beamSpotX0_ , "beamSpotX0_/D" ) ;
  TPat_->Branch("beamSpotY0_", &beamSpotY0_ , "beamSpotY0_/D" ) ; 
  TPat_->Branch("beamSpotZ0_", &beamSpotZ0_ , "beamSpotZ0_/D" ) ; 

  TPat_->Branch("n1ryVertices_", &n1ryVertices_,       "n1ryVertices_/I") ; 
  
  TPat_->Branch("xVtx_",         "std::vector<double>", &xVtx_         ) ;
  TPat_->Branch("yVtx_",         "std::vector<double>", &vtxXError_    ) ;   
  TPat_->Branch("vtxXError_",    "std::vector<double>", &vtxYError_    ) ; 
  TPat_->Branch("vtxYError_",    "std::vector<double>", &yVtx_         ) ; 
  TPat_->Branch("xVtxBeamSpot_", "std::vector<double>", &xVtxBeamSpot_ ) ;
  TPat_->Branch("yVtxBeamSpot_", "std::vector<double>", &yVtxBeamSpot_ ) ; 
  TPat_->Branch("zVtx_",         "std::vector<double>", &zVtx_         ) ; 
  TPat_->Branch("nTrkFromVtx_",  "std::vector<int>",    &nTrkFromVtx_  ) ; 
  TPat_->Branch("nTrkVtxWithWt_","std::vector<int>",    &nTrkVtxWithWt_) ; 
  TPat_->Branch("vtxchi2_",      "std::vector<double>", &vtxchi2_      ) ;
  TPat_->Branch("vtxNormChi2_",  "std::vector<double>", &vtxNormChi2_  ) ;  
  TPat_->Branch("vtxndof_",      "std::vector<double>", &vtxndof_      ) ; 
  TPat_->Branch("vtxprob_",      "std::vector<double>", &vtxprob_      ) ;  

  TPat_->Branch("genWCh_",     "std::vector<int   >", &genWCh_    ) ;     
  TPat_->Branch("genWPt_",     "std::vector<double>", &genWPt_    ) ;
  TPat_->Branch("genWEta_",    "std::vector<double>", &genWEta_   ) ; 
  TPat_->Branch("genWPhi_",    "std::vector<double>", &genWPhi_   ) ; 
  TPat_->Branch("genWE_",      "std::vector<double>", &genWE_     ) ; 

  TPat_->Branch("genMuCh_",     "std::vector<int   >", &genMuCh_    ) ;     
  TPat_->Branch("genMuPt_",     "std::vector<double>", &genMuPt_    ) ;
  TPat_->Branch("genMuEta_",    "std::vector<double>", &genMuEta_   ) ; 
  TPat_->Branch("genMuPhi_",    "std::vector<double>", &genMuPhi_   ) ; 
  TPat_->Branch("genMuE_",      "std::vector<double>", &genMuE_     ) ; 

  TPat_->Branch("genElCh_",     "std::vector<int   >", &genElCh_    ) ;     
  TPat_->Branch("genElPt_",     "std::vector<double>", &genElPt_    ) ;
  TPat_->Branch("genElEta_",    "std::vector<double>", &genElEta_   ) ; 
  TPat_->Branch("genElPhi_",    "std::vector<double>", &genElPhi_   ) ; 
  TPat_->Branch("genElE_",      "std::vector<double>", &genElE_     ) ; 

  TPat_->Branch("genNuPt_",     "std::vector<double>", &genNuPt_    ) ;
  TPat_->Branch("genNuEta_",    "std::vector<double>", &genNuEta_   ) ; 
  TPat_->Branch("genNuPhi_",    "std::vector<double>", &genNuPhi_   ) ; 
  TPat_->Branch("genNuE_",      "std::vector<double>", &genNuE_     ) ; 

  TPat_->Branch("genGammaPt_",  "std::vector<double>", &genGammaPt_ ) ;
  TPat_->Branch("genGammaEta_", "std::vector<double>", &genGammaEta_) ; 
  TPat_->Branch("genGammaPhi_", "std::vector<double>", &genGammaPhi_) ; 
  TPat_->Branch("genGammaE_",   "std::vector<double>", &genGammaE_  ) ; 

  TPat_->Branch("genWGammaPt_",  "std::vector<double>", &genWGammaPt_ ) ;
  TPat_->Branch("genWGammaEta_", "std::vector<double>", &genWGammaEta_) ; 
  TPat_->Branch("genWGammaPhi_", "std::vector<double>", &genWGammaPhi_) ; 
  TPat_->Branch("genWGammaE_",   "std::vector<double>", &genWGammaE_  ) ; 

  TPat_->Branch("nRecoMuons_",             &nRecoMuons_,          "nRecoMuons_/I"         ); 
  TPat_->Branch("isInnerTrackNonNull_",    "std::vector<bool  >", &isInnerTrackNonNull_   ); 
  TPat_->Branch("isGlobalTrackNonNull_",   "std::vector<bool  >", &isGlobalTrackNonNull_  ); 
  TPat_->Branch("noOfValidHits_",          "std::vector<int   >", &noOfValidHits_         ); 
  TPat_->Branch("noOfLostHits_",           "std::vector<int   >", &noOfLostHits_          ); 
  TPat_->Branch("recoMuonPt_",             "std::vector<double>", &recoMuonPt_            );
  TPat_->Branch("genMuonMother_",          "std::vector<int   >", &genMuonMother_         );
  TPat_->Branch("genMuonPt_",              "std::vector<double>", &genMuonPt_             );
  TPat_->Branch("recoMuonPtError_",        "std::vector<double>", &recoMuonPtError_       );
  TPat_->Branch("recoMuonPx_",             "std::vector<double>", &recoMuonPx_            );
  TPat_->Branch("recoMuonPy_",             "std::vector<double>", &recoMuonPy_            );
  TPat_->Branch("recoMuonPz_",             "std::vector<double>", &recoMuonPz_            );
  TPat_->Branch("recoMuonEta_",            "std::vector<double>", &recoMuonEta_           );
  TPat_->Branch("genMuonEta_",             "std::vector<double>", &genMuonEta_            );
  TPat_->Branch("recoMuonEtaError_",       "std::vector<double>", &recoMuonEtaError_      );
  TPat_->Branch("recoMuonPhi_",            "std::vector<double>", &recoMuonPhi_           );
  TPat_->Branch("genMuonPhi_",             "std::vector<double>", &genMuonPhi_            );
  TPat_->Branch("recoMuonPhiError_",       "std::vector<double>", &recoMuonPhiError_      );
  TPat_->Branch("recoMuonCharge_",         "std::vector<double>", &recoMuonCharge_        );
  TPat_->Branch("recoMuonEnergy_",         "std::vector<double>", &recoMuonEnergy_        );
  TPat_->Branch("recoMuonIsoTrack_",       "std::vector<double>", &recoMuonIsoTrack_      );
  TPat_->Branch("recoMuonIsoEcal_",        "std::vector<double>", &recoMuonIsoEcal_       );
  TPat_->Branch("recoMuonIsoHcal_",        "std::vector<double>", &recoMuonIsoHcal_       );  
  TPat_->Branch("recoMuondxy_",            "std::vector<double>", &recoMuondxy_           );  
  TPat_->Branch("recoMuondxyError_",       "std::vector<double>", &recoMuondxyError_      );  
  TPat_->Branch("recoMuondxySigma_",       "std::vector<double>", &recoMuondxySigma_      );  
  TPat_->Branch("recoMuondz_",             "std::vector<double>", &recoMuondz_            );  
  TPat_->Branch("recoMuondzError_",        "std::vector<double>", &recoMuondzError_       );  
  TPat_->Branch("recoMuonNormalizedChi2_", "std::vector<double>", &recoMuonNormalizedChi2_);  

  nRecoMuons_ = 0 ;

  TPat_->Branch("nRecoPhotons_",           &nRecoPhotons_,         "nRecoPhotons_/I"      ); 
  TPat_->Branch("genPhotonMother_",        "std::vector<int >",   &genPhotonMother_       ); 
  TPat_->Branch("genPhotonGrandmother_",   "std::vector<int >",   &genPhotonGrandmother_  ); 
  TPat_->Branch("recoPhotonHasPixelSeed_", "std::vector<bool>",   &recoPhotonHasPixelSeed_);
  TPat_->Branch("isEBPhoton_",             "std::vector<bool>",   &isEBPhoton_            );
  TPat_->Branch("isEEPhoton_",             "std::vector<bool>",   &isEEPhoton_            );
  TPat_->Branch("isEBGapPhoton_",          "std::vector<bool>",   &isEBGapPhoton_         );
  TPat_->Branch("isEEGapPhoton_",          "std::vector<bool>",   &isEEGapPhoton_         );
  TPat_->Branch("isEBEEGapPhoton_",        "std::vector<bool>",   &isEBEEGapPhoton_       );

  TPat_->Branch("recoPhotonHoverE_",       "std::vector<double>", &recoPhotonHoverE_        ); 
  TPat_->Branch("recoPhotonHDepth1overE_", "std::vector<double>", &recoPhotonHDepth1overE_  ); 
  TPat_->Branch("recoPhotonHDepth2overE_", "std::vector<double>", &recoPhotonHDepth2overE_  );  
  TPat_->Branch("recoPhotonR9_",           "std::vector<double>", &recoPhotonR9_            );
  TPat_->Branch("recoPhotonSigmaEtaEta_",  "std::vector<double>", &recoPhotonSigmaEtaEta_   );
  TPat_->Branch("recoPhotonSigmaIetaIeta_","std::vector<double>", &recoPhotonSigmaIetaIeta_ ); 
  TPat_->Branch("recoPhotonR1x5_",         "std::vector<double>", &recoPhotonR1x5_          ); 
  TPat_->Branch("recoPhotonR2x5_",         "std::vector<double>", &recoPhotonR2x5_          ); 
  TPat_->Branch("recoPhotonE1x5_",         "std::vector<double>", &recoPhotonE1x5_          ); 
  TPat_->Branch("recoPhotonE2x5_",         "std::vector<double>", &recoPhotonE2x5_          ); 
  TPat_->Branch("recoPhotonE3x3_",         "std::vector<double>", &recoPhotonE3x3_          ); 
  TPat_->Branch("recoPhotonE5x5_",         "std::vector<double>", &recoPhotonE5x5_          ); 
  TPat_->Branch("recoPhotonMaxEXtal_",     "std::vector<double>", &recoPhotonMaxEXtal_      ); 


  TPat_->Branch("recoPhotonTrkSumPtHollowConeDR04_",  "std::vector<double>", &recoPhotonTrkSumPtHollowConeDR04_ );
  TPat_->Branch("recoPhotonEcalRecHitSumEtConeDR04_", "std::vector<double>", &recoPhotonEcalRecHitSumEtConeDR04_);
  TPat_->Branch("recoPhotonHcalTowerSumEtConeDR04_",  "std::vector<double>", &recoPhotonHcalTowerSumEtConeDR04_ );
  TPat_->Branch("nTrkSolidConeDR04_",                 "std::vector<int >",   &nTrkSolidConeDR04_                );
  TPat_->Branch("nTrkHollowConeDR04_",                "std::vector<int >",   &nTrkHollowConeDR04_               );

  TPat_->Branch("recoPhotonTrkSumPtHollowConeDR03_",  "std::vector<double>", &recoPhotonTrkSumPtHollowConeDR03_ );
  TPat_->Branch("recoPhotonEcalRecHitSumEtConeDR03_", "std::vector<double>", &recoPhotonEcalRecHitSumEtConeDR03_);
  TPat_->Branch("recoPhotonHcalTowerSumEtConeDR03_",  "std::vector<double>", &recoPhotonHcalTowerSumEtConeDR03_ );
  TPat_->Branch("nTrkSolidConeDR03_",                 "std::vector<int >",   &nTrkSolidConeDR03_                );
  TPat_->Branch("nTrkHollowConeDR03_",                "std::vector<int >",   &nTrkHollowConeDR03_               );
  TPat_->Branch("recoPhotonDelRClosestTrack_",        "std::vector<double>", &recoPhotonDelRClosestTrack_       );

  TPat_->Branch("recoPhotonPt_",           "std::vector<double>", &recoPhotonPt_           );
  TPat_->Branch("genPhotonPt_",            "std::vector<double>", &genPhotonPt_            );
  TPat_->Branch("recoPhotonEta_",          "std::vector<double>", &recoPhotonEta_          );
  TPat_->Branch("genPhotonEta_",           "std::vector<double>", &genPhotonEta_           );
  TPat_->Branch("recoPhotonPhi_",          "std::vector<double>", &recoPhotonPhi_          );
  TPat_->Branch("genPhotonPhi_",           "std::vector<double>", &genPhotonPhi_           );
  TPat_->Branch("recoPhotonPx_",           "std::vector<double>", &recoPhotonPx_           );
  TPat_->Branch("recoPhotonPy_",           "std::vector<double>", &recoPhotonPy_           );
  TPat_->Branch("recoPhotonPz_",           "std::vector<double>", &recoPhotonPz_           );
  TPat_->Branch("recoPhotonEnergy_",       "std::vector<double>", &recoPhotonEnergy_       );
  TPat_->Branch("recoPhotonCharge_",       "std::vector<double>", &recoPhotonCharge_       );

  nRecoPhotons_ = 0 ;
        
  TPat_->Branch("nRecoElectrons",      &nRecoElectrons,       "nRecoElectrons/I"   ); 
  TPat_->Branch("recoElectronPt",      "std::vector<double>", &recoElectronPt      );
  TPat_->Branch("recoElectronPx",      "std::vector<double>", &recoElectronPx      );
  TPat_->Branch("recoElectronPy",      "std::vector<double>", &recoElectronPy      );
  TPat_->Branch("recoElectronPz",      "std::vector<double>", &recoElectronPz      );
  TPat_->Branch("recoElectronEta",     "std::vector<double>", &recoElectronEta     );
  TPat_->Branch("recoElectronPhi",     "std::vector<double>", &recoElectronPhi     );
  TPat_->Branch("recoElectronCharge",  "std::vector<double>", &recoElectronCharge  );
  TPat_->Branch("recoElectronEnergy",  "std::vector<double>", &recoElectronEnergy  );
  TPat_->Branch("recoElectronIsoTrack","std::vector<double>", &recoElectronIsoTrack);
  TPat_->Branch("recoElectronIsoEcal", "std::vector<double>", &recoElectronIsoEcal );
  TPat_->Branch("recoElectronIsoHcal", "std::vector<double>", &recoElectronIsoHcal );
  
  nRecoElectrons = 0 ;
  
  TPat_->Branch("nJets_",      &nJets_,               "nJets_/I"  ) ;
  TPat_->Branch("jet_px_",     "std::vector<double>", &jet_px_   ) ;
  TPat_->Branch("jet_py_",     "std::vector<double>", &jet_py_   ) ;
  TPat_->Branch("jet_pz_",     "std::vector<double>", &jet_pz_   ) ;
  TPat_->Branch("jet_pt_",     "std::vector<double>", &jet_pt_   ) ;
  TPat_->Branch("jet_eta_",    "std::vector<double>", &jet_eta_  ) ;
  TPat_->Branch("jet_phi_",    "std::vector<double>", &jet_phi_  ) ;
  TPat_->Branch("jet_E_",      "std::vector<double>", &jet_E_    ) ; 
  TPat_->Branch("jet_P_",      "std::vector<double>", &jet_P_    ) ; 
  TPat_->Branch("jet_charge_", "std::vector<double>", &jet_charge_) ;

  TPat_->Branch("maxEInEMTowers_ ", "std::vector<double>", &maxEInEMTowers_ ) ; 
  TPat_->Branch("maxEInHadTowers_", "std::vector<double>", &maxEInHadTowers_) ; 
  TPat_->Branch("hadEFraction_   ", "std::vector<double>", &hadEFraction_   ) ; 
  TPat_->Branch("emEFraction_    ", "std::vector<double>", &emEFraction_    ) ; 
  TPat_->Branch("hadEInHB_       ", "std::vector<double>", &hadEInHB_       ) ; 
  TPat_->Branch("hadEInHE_       ", "std::vector<double>", &hadEInHE_       ) ; 
  TPat_->Branch("hadEInHO_       ", "std::vector<double>", &hadEInHO_       ) ; 
  TPat_->Branch("hadEInHF_       ", "std::vector<double>", &hadEInHF_       ) ; 
  TPat_->Branch("emEnergyInEB_   ", "std::vector<double>", &emEnergyInEB_   ) ; 
  TPat_->Branch("emEnergyInEE_   ", "std::vector<double>", &emEnergyInEE_   ) ; 
  TPat_->Branch("emEnergyInHF_   ", "std::vector<double>", &emEnergyInHF_   ) ; 
  TPat_->Branch("towersArea_     ", "std::vector<double>", &towersArea_     ) ; 
  TPat_->Branch("jet_n60_        ", "std::vector<double>", &jet_n60_        ) ; 
  TPat_->Branch("jet_n90_        ", "std::vector<double>", &jet_n90_        ) ; 
  TPat_->Branch("jet_fHPD_       ", "std::vector<double>", &jet_fHPD_       ) ;

  nJets_ = 0 ;
   
  TPat_->Branch("Met_x",        &Met_x        ,"Met_x       /D");
  TPat_->Branch("Met_y",        &Met_y        ,"Met_y       /D");
  TPat_->Branch("Met",          &Met          ,"Met         /D");
  TPat_->Branch("Met_phi",      &Met_phi      ,"Met_phi     /D");
  TPat_->Branch("sumet",        &sumet        ,"sumet   /D");
  TPat_->Branch("genMet_x",     &genMet_x     ,"genMet_x    /D");
  TPat_->Branch("genMet_y",     &genMet_y     ,"genMet_y    /D");
  TPat_->Branch("genMet",       &genMet       ,"genMet      /D");
  TPat_->Branch("genMet_phi",   &genMet_phi   ,"genMet_phi  /D");
  TPat_->Branch("gensumet",     &gensumet     ,"gensumet/D"); 
  
  TPat_->Branch("uncorAllMet_pt",   &uncorAllMet_pt  ,  "uncorAllMet_pt  /D");
  TPat_->Branch("uncorAllMet_phi",  &uncorAllMet_phi ,  "uncorAllMet_phi /D");
  TPat_->Branch("uncorJesMet_pt",   &uncorJesMet_pt  ,  "uncorJesMet_pt  /D");
  TPat_->Branch("uncorJesMet_phi",  &uncorJesMet_phi ,  "uncorJesMet_phi /D");
  TPat_->Branch("uncorMuonMet_pt",  &uncorMuonMet_pt ,  "uncorMuonMet_pt /D");
  TPat_->Branch("uncorMuonMet_phi", &uncorMuonMet_phi,  "uncorMuonMet_phi/D");
  
  Met_x = Met_y = Met = Met_phi = sumet = -100. ;
  genMet_x = genMet_y = genMet = genMet_phi = gensumet = -100. ;
  uncorAllMet_pt = uncorAllMet_phi = -100. ;
  uncorJesMet_pt = uncorJesMet_phi = -100. ;
  uncorMuonMet_pt = uncorMuonMet_phi = -100. ; 

} // end buildTree 

// ------------ method called once each job just after ending the event loop  ------------
void WGEventSelector::endJob() {
  outfile_ << "Total events analysed: " << nEvt_        << "\n" ; 
  outfile_ << "nL1Accept_ = "           << nL1Accept_   << "\n" ; 
  outfile_ << "nHLTAccept_ = "          << nHLTAccept_  << "\n" ; 
  outfile_ << "nHLTwasRun_ = "          << nHLTwasRun_  << "\n" ; 
  outfile_ << "nHLTErrors_ = "          << nHLTErrors_  << "\n" ; 
  outfile_ << "nRecoMuons_ = "          << nRecoMuons_  << "\n" ; 
  outfile_ << "nRecoPhotons_ = "          << nRecoPhotons_  << "\n" ; 
  outfile_ << "nRecoElectrons = "          << nRecoElectrons  << "\n" ; 
  outfile_.close() ; 
}

//define this as a plug-in
DEFINE_FWK_MODULE(WGEventSelector);
