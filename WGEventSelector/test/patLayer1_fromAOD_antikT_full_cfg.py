from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.source.fileNames = [
		#"/store/relval/CMSSW_3_4_2/RelValTTbar/GEN-SIM-RECO/STARTUP3X_V15-v1/0011/CE7222DF-8213-DF11-A1EA-001A92810AD6.root"
		"file:/tmp/devdatta/DCBDDC34-A9CA-DE11-AE97-001E4F3F28E4.root"
		] 

from PhysicsTools.PatAlgos.tools.jetTools import *

# uncomment the following lines to switch the standard jet collection
# in PAT (iterativeCone5) to ak5 (anti-kt cone = 0.5)
switchJetCollection(process, 
                    cms.InputTag('ak5CaloJets'),   
                    doJTA            = True,            
                    doBTagging       = True,            
                    jetCorrLabel     = ('AK5','Calo'),  
                    doType1MET       = True,
                    genJetCollection = cms.InputTag("kt4GenJets"),
                    doJetID          = True,
                    jetIdLabel       = "ak5"
                    ) 

process.p = cms.Path(
                process.patDefaultSequence  
            )

# In addition you usually want to change the following parameters:
#
#   process.GlobalTag.globaltag =  ...     ## (according to https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions)
#   process.source.fileNames = [ ... ]     ## (e.g. 'file:AOD.root')
process.maxEvents.input = 20               ## (e.g. -1 to run on all events)
#   process.out.outputCommands = [ ... ]   ## (e.g. taken from PhysicsTools/PatAlgos/python/patEventContent_cff.py)
#   process.out.fileName = ...             ## (e.g. 'myTuple.root')
process.options.wantSummary = True        ## (to suppress the long output at the end of the job)

