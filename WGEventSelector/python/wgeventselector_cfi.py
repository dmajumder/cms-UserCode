import FWCore.ParameterSet.Config as cms

wgeventselector = cms.EDAnalyzer('WGEventSelector',
	isMC              = cms.bool(True), 
	genParticleTag    = cms.InputTag("genParticles"),
	l1gtrrTag         = cms.InputTag('gtDigis'),
	l1ObjectMapTag    = cms.InputTag('hltL1GtObjectMap','','HLT'),
	triggerResultsTag = cms.InputTag('TriggerResults','','HLT'),
	genTag            = cms.InputTag("genParticles"), 
	beamSpotTag       = cms.InputTag("offlineBeamSpot"), 
	vertexTag         = cms.InputTag("offlinePrimaryVertices"),
	muonTag           = cms.InputTag("cleanPatMuons"),
	phoTag            = cms.InputTag("cleanPatPhotons"),
	trackTag          = cms.InputTag("generalTracks"),
	electronTag       = cms.InputTag("cleanPatElectrons"),
	jetTag            = cms.InputTag("cleanPatJets"),
	metTag            = cms.InputTag("patMETs"), 
        pvTrkMinWt        = cms.double("0.5"), 
	maxPtClosestTrack = cms.double("1.5"), 
	phPtMin           = cms.double("20."),      
	phEtaMax          = cms.double("1.5"),    
	muPtMin           = cms.double("15."), 
	muEtaMax          = cms.double("2."), 
	outfilename       = cms.untracked.string("pat342_wgeventselector.out") 
	)
