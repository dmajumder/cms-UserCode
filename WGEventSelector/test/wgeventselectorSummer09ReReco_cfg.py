import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.options   = cms.untracked.PSet( 
		wantSummary = cms.untracked.bool(False),
		SkipEvent = cms.untracked.vstring('ProductNotFound') 
		)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) 

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    #/RelValGammaJets_Pt_80_120/CMSSW_3_3_2-STARTUP31X_V8-v1/GEN-SIM-RECO
    #'/store/relval/CMSSW_3_3_2/RelValGammaJets_Pt_80_120/GEN-SIM-RECO/STARTUP31X_V8-v1/0005/48280311-BFC7-DE11-A35F-001D09F24691.root' 
    #'file:/tmp/devdatta/DCBDDC34-A9CA-DE11-AE97-001E4F3F28E4.root'
    'file:Summer09PhotonJet_Pt30to50GEN-SIM-RECOMC_31X_V9_7TeV-v1.root' 
    )
)

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP3X_V11::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.tools.jetTools import *

from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
run33xOnReRecoMC( process, "ak5GenJets" ) 

#from PhysicsTools.PatAlgos.tools.coreTools import *
#print "*********************************************************************"
#print "Switching all processes to use the anti-kT algorithm by default."
#print "Switch the jet collection to your desired algorithm if this is not"
#print "what you want to use. Note that L7Parton correction are taken from"
#print "SC5 instead of AK5. This is an intermediate solution for the time "
#print "being."
#print "*********************************************************************"
#switchJetCollection(process, 
#    	cms.InputTag('ak5CaloJets'),   
#    	doJTA            = True,            
#    	doBTagging       = True,            
#    	jetCorrLabel     = ('AK5','Calo'),  
#    	doType1MET       = True, 
#    	genJetCollection = cms.InputTag("ak5GenJets"), 
#	doJetID          = True,
#	jetIdLabel       = "ak5"
#    	) 

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent 

process.load('WGammaAnalyzer.TriggerSelector.triggerselector_cfi') 
process.load('WGammaAnalyzer.WGEventSelector.wgeventselector_cfi') 

process.TFileService = cms.Service("TFileService", fileName = cms.string("/tmp/devdatta/pat342_wgeventselector_Summer09PhotonJetPt30to50.root") ) 
#process.TFileService = cms.Service("TFileService", fileName = cms.string("pat342_wgeventselector.root") )  

#process.p = cms.Path(process.triggerselector+process.patDefaultSequence+process.wgeventselector)  
#process.p = cms.Path(process.patDefaultSequence)  
process.p = cms.Path(process.patDefaultSequence+process.wgeventselector)  
