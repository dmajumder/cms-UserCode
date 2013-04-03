import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import copy

process = cms.Process("Demo")

options = VarParsing ('python')

options.register('outFilename',
		'outfile_Zbb.root',
		VarParsing.multiplicity.singleton,
		VarParsing.varType.string,
		"Output file name"
		)

options.register('reportEvery',
		1,
		VarParsing.multiplicity.singleton,
		VarParsing.varType.int,
		"Report every N events (default is N=10)"
		)

options.register('wantSummary',
		True,
		VarParsing.multiplicity.singleton,
		VarParsing.varType.bool,
		"Print out trigger and timing summary"
		)

options.register ('useData',
		False,
		VarParsing.multiplicity.singleton,
		VarParsing.varType.int,
		'Run this on real data')

options.setDefault('maxEvents', -1)

options.register('jetRadius',
		0.8,
		VarParsing.multiplicity.singleton,
		VarParsing.varType.float,
		"Distance parameter R for jet clustering (default is 0.8)"
		)

options.parseArguments() 

print options

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger = cms.Service("MessageLogger",
		destinations = cms.untracked.vstring(
			'detailedInfo',
			),
		detailedInfo = cms.untracked.PSet(
			threshold = cms.untracked.string('INFO'), 
			), 
		categories = cms.untracked.vstring('NSubjets',
			'NJets',
			'JetFlavour', 
			'BFromGSplit'
			),
		JetFlavour  = cms.untracked.PSet(limit=cms.untracked.int32(0)), 
		BFromGSplit = cms.untracked.PSet(limit=cms.untracked.int32(0)),    
		suppressInfo = cms.untracked.vstring("HiggsTagValidation") 
		)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(options.wantSummary) )

listOfFiles = cms.untracked.vstring() 

listOfFiles.extend( [
	'file:/afs/cern.ch/work/d/devdatta/CMSREL/CMSSW_5_3_8_patch1_B2GTopLikeBSM53X_v3/src/TopQuarkAnalysis/TopPairBSM/test/tlbsm_53x_v3_Slim_mc_fat.root'
		] ) 

process.source = cms.Source("PoolSource",
		# replace 'myfile.root' with the source file you want to use
		fileNames = listOfFiles 
		)

## Output file
process.TFileService = cms.Service("TFileService",
		fileName = cms.string(options.outFilename)
		)

process.HiggsTagValidation = cms.EDAnalyzer('HiggsTagValidation', 
		IsData                    = cms.bool(options.useData), 
		UseEventWeight            = cms.bool(False),
		BosonPdgId                = cms.int32(25),
		IpTrackTag                = cms.InputTag('impactParameterTagInfos'),
		GenParticleTag            = cms.InputTag('prunedGenParticles'), 
		JetsTag                   = cms.InputTag('goodPatJetsCA8PF'), 
		GroomedJetsTag            = cms.InputTag('goodPatJetsCA8PrunedPFPacked'),
		#GroomedBasicJetsTag       = cms.InputTag('caPFJetsPrunedCompound'),
		SubJetsTag                = cms.InputTag('selectedPatJetsCA8PrunedSubjetsPF'),
		PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
		JetRadius                 = cms.double(options.jetRadius), 
		JetFlavourPdgId           = cms.vint32(1,2,3,4,5,21), 
		JetPtMin                  = cms.double(50.), 
		JetAbsEtaMax              = cms.double(2.4), 
		JetMassMin                = cms.double(75.), 
		JetMassMax                = cms.double(155.), 
		DRSubjetMatch             = cms.double(0.8) 
		)


process.p = cms.Path(process.HiggsTagValidation)
