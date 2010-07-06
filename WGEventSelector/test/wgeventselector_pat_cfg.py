from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *

process.source.fileNames = [
		#'file:/tmp/devdatta/BaurWGAMMALO_WpMunuG_ATGC_K_01_l_02_GEN_FASTSIM.root'
		'rfio:/castor/cern.ch/user/m/mazumdar/wgamma/LO_FASTSIM.root'
		]

process.maxEvents.input = -1 

process.options.wantSummary = False 

process.load('WGammaAnalyzer.WGEventSelector.wgeventselector_cfi')

process.p = cms.Path(
		process.patDefaultSequence+process.wgeventselector
		)

process.TFileService = cms.Service("TFileService", fileName = cms.string("pat356_baurLO_KM.root") ) 

