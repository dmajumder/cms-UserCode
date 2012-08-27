# Auto generated configuration file
# using: 
# Revision: 1.381.2.2 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/EightTeV/Bprime/BprimeBprimeToBHBH_HToBB_M_800_TuneZ2star_8TeV_madgraph_cff.py --filein lhe:6159 --step GEN,FASTSIM,HLT:7E33v2 --conditions START53_V7A::All --pileup 2012_Summer_inTimeOnly --datamix NODATAMIXER --beamspot Realistic8TeVCollision --eventcontent AODSIM --datatier AODSIM --no_exec -n -1
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('FastSimulation.Configuration.EventContent_cff')
process.load('FastSimulation.PileUpProducer.PileUpSimulator_2012_Summer_inTimeOnly_cff')
process.load('FastSimulation.Configuration.Geometries_START_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('FastSimulation.Configuration.FamosSequences_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedParameters_cfi')
process.load('HLTrigger.Configuration.HLT_7E33v2_Famos_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("LHESource",
    fileNames = cms.untracked.vstring('/store/lhe/6159/8TeV_bpbpbar_800_8._run10466_unweighted_events_qcut70_mgPostv2.lhe', 
        '/store/lhe/6159/8TeV_bpbpbar_800_8._run13986_unweighted_events_qcut70_mgPostv2.lhe', 
        '/store/lhe/6159/8TeV_bpbpbar_800_8._run14150_unweighted_events_qcut70_mgPostv2.lhe', 
        '/store/lhe/6159/8TeV_bpbpbar_800_8._run3977_unweighted_events_qcut70_mgPostv2.lhe', 
        '/store/lhe/6159/8TeV_bpbpbar_800_8._run913_unweighted_events_qcut70_mgPostv2.lhe')
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.381.2.2 $'),
    annotation = cms.untracked.string('Configuration/GenProduction/python/EightTeV/Bprime/BprimeBprimeToBHBH_HToBB_M_800_TuneZ2star_8TeV_madgraph_cff.py nevts:-1'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.AODSIMoutput = cms.OutputModule("PoolOutputModule",
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    outputCommands = process.AODSIMEventContent.outputCommands,
    fileName = cms.untracked.string('BprimeBprimeToBHBH_HToBB_M_800_TuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('AODSIM')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.famosSimHits.SimulateCalorimetry = True
process.famosSimHits.SimulateTracking = True
process.simulation = cms.Sequence(process.simulationWithFamos)
process.HLTEndSequence = cms.Sequence(process.reconstructionWithFamos)
process.Realistic8TeVCollisionVtxSmearingParameters.type = cms.string("BetaFunc")
process.famosSimHits.VertexGenerator = process.Realistic8TeVCollisionVtxSmearingParameters
process.famosPileUp.VertexGenerator = process.Realistic8TeVCollisionVtxSmearingParameters
# customise the HLT menu for running on MC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC
process = customizeHLTforMC(process)

process.GlobalTag.globaltag = 'START53_V7A::All'

process.generator = cms.EDFilter("Pythia6HadronizerFilter",
    jetMatching = cms.untracked.PSet(
        MEMAIN_showerkt = cms.double(0),
        MEMAIN_nqmatch = cms.int32(-1),
        MEMAIN_minjets = cms.int32(-1),
        MEMAIN_qcut = cms.double(-1),
        MEMAIN_excres = cms.string(''),
        MEMAIN_etaclmax = cms.double(5.0),
        outTree_flag = cms.int32(0),
        scheme = cms.string('Madgraph'),
        MEMAIN_maxjets = cms.int32(-1),
        mode = cms.string('auto')
    ),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(8000.0),
    maxEventsToPrint = cms.untracked.int32(0),
    PythiaParameters = cms.PSet(
        pythiaUESettings = cms.vstring('MSTU(21)=1     ! Check on possible errors during program execution', 
            'MSTJ(22)=2     ! Decay those unstable particles', 
            'PARJ(71)=10 .  ! for which ctau  10 mm', 
            'MSTP(33)=0     ! no K factors in hard cross sections', 
            'MSTP(2)=1      ! which order running alphaS', 
            'MSTP(51)=10042 ! structure function chosen (external PDF CTEQ6L1)', 
            'MSTP(52)=2     ! work with LHAPDF', 
            'PARP(82)=1.921 ! pt cutoff for multiparton interactions', 
            'PARP(89)=1800. ! sqrts for which PARP82 is set', 
            'PARP(90)=0.227 ! Multiple interactions: rescaling power', 
            'MSTP(95)=6     ! CR (color reconnection parameters)', 
            'PARP(77)=1.016 ! CR', 
            'PARP(78)=0.538 ! CR', 
            'PARP(80)=0.1   ! Prob. colored parton from BBR', 
            'PARP(83)=0.356 ! Multiple interactions: matter distribution parameter', 
            'PARP(84)=0.651 ! Multiple interactions: matter distribution parameter', 
            'PARP(62)=1.025 ! ISR cutoff', 
            'MSTP(91)=1     ! Gaussian primordial kT', 
            'PARP(93)=10.0  ! primordial kT-max', 
            'MSTP(81)=21    ! multiple parton interactions 1 is Pythia default', 
            'MSTP(82)=4     ! Defines the multi-parton model'),
        processParameters = cms.vstring('PMAS(25,1)=125.00D0    !mass of Higgs', 
            'MSTP(1) = 4', 
            'MSEL=7         ! User defined processes', 
            'MWID(7)=2', 
            'MSTJ(1)=1       ! Fragmentation/hadronization on or off', 
            'MSTP(61)=1      ! Parton showering on or off', 
            'PMAS(5,1)=4.8   ! b quark mass', 
            'PMAS(6,1)=172.5 ! t quark mass', 
            'PMAS(7,1) = 800.0D0     ! bprime quarks mass', 
            'PMAS(7,2) = 8.000D0    ! bprime quark width', 
            'PMAS(7,3) = 80.00D0 ! Max value above which the BW shape is truncated', 
            'VCKM(1,1) = 0.97414000D0', 
            'VCKM(1,2) = 0.22450000D0', 
            'VCKM(1,3) = 0.00420000D0', 
            'VCKM(1,4) = 0.02500000D0', 
            'VCKM(2,1) = 0.22560000D0', 
            'VCKM(2,2) = 0.97170000D0', 
            'VCKM(2,3) = 0.04109000D0', 
            'VCKM(2,4) = 0.05700000D0', 
            'VCKM(3,1) = 0.00100000D0', 
            'VCKM(3,2) = 0.06200000D0', 
            'VCKM(3,3) = 0.91000000D0', 
            'VCKM(3,4) = 0.41000000D0', 
            'VCKM(4,1) = 0.01300000D0', 
            'VCKM(4,2) = 0.04000000D0', 
            'VCKM(4,3) = 0.41000000D0', 
            'VCKM(4,4) = 0.91000000D0', 
            'MDME(56,1)=0     ! g b4', 
            'MDME(57,1)=0     ! gamma b4', 
            'MDME(58,1)=0     ! Z0 b', 
            'MDME(59,1)=0     ! W u', 
            'MDME(60,1)=0     ! W c', 
            'MDME(61,1)=0     ! W t', 
            'MDME(62,1)=0     ! W t4', 
            'KFDP(63,2)=5     ! defines H0 b', 
            'MDME(63,1)=1     ! h0 b4', 
            'MDME(64,1)=-1    ! H- c', 
            'MDME(65,1)=-1    ! H- t', 
            'BRAT(56)  = 0.0D0', 
            'BRAT(57)  = 0.0D0', 
            'BRAT(58)  = 0.0D0', 
            'BRAT(59)  = 0.0D0', 
            'BRAT(60)  = 0.0D0', 
            'BRAT(61)  = 0.0D0', 
            'BRAT(62)  = 0.0D0', 
            'BRAT(63)  = 1.0D0', 
            'BRAT(64)  = 0.0D0', 
            'BRAT(65)  = 0.0D0', 
            'MDME(210,1)=0     !Higgs decay into dd', 
            'MDME(211,1)=0     !Higgs decay into uu', 
            'MDME(212,1)=0     !Higgs decay into ss', 
            'MDME(213,1)=0     !Higgs decay into cc', 
            'MDME(214,1)=1     !Higgs decay into bb', 
            'MDME(215,1)=0     !Higgs decay into tt', 
            'MDME(216,1)=0     !Higgs decay into', 
            'MDME(217,1)=0     !Higgs decay into Higgs decay', 
            'MDME(218,1)=0     !Higgs decay into e nu e', 
            'MDME(219,1)=0     !Higgs decay into mu nu mu', 
            'MDME(220,1)=0     !Higgs decay into tau nu tau', 
            'MDME(221,1)=0     !Higgs decay into Higgs decay', 
            'MDME(222,1)=0     !Higgs decay into g g', 
            'MDME(223,1)=0     !Higgs decay into gam gam', 
            'MDME(224,1)=0     !Higgs decay into gam Z', 
            'MDME(225,1)=0     !Higgs decay into Z Z', 
            'MDME(226,1)=0     !Higgs decay into W W'),
        parameterSets = cms.vstring('pythiaUESettings', 
            'processParameters')
    )
)


process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen_genonly)
process.reconstruction = cms.Path(process.reconstructionWithFamos)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.AODSIMoutput_step = cms.EndPath(process.AODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.reconstruction,process.AODSIMoutput_step])
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 

