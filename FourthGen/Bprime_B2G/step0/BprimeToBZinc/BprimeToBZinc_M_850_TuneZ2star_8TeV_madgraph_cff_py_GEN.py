# Auto generated configuration file
# using: 
# Revision: 1.372.2.3 
# Source: /local/reps/CMSSW.admin/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/EightTeV/Bprime_B2G/BprimeToBZinc_M_850_TuneZ2star_8TeV_madgraph_cff.py --filein lhe:6476 -s GEN --conditions START52_V9::All --beamspot Realistic8TeVCollision --datatier GEN-SIM --eventcontent RAWSIM -n 1000 --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('GEN')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# Input source
process.source = cms.Source("LHESource",
    fileNames = cms.untracked.vstring('/store/lhe/6476/8TeV_bp_850_8.5_run17155_unweighted_events_qcut0_mgPostv2.lhe')
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.372.2.3 $'),
    annotation = cms.untracked.string('Configuration/GenProduction/python/EightTeV/Bprime_B2G/BprimeToBZinc_M_850_TuneZ2star_8TeV_madgraph_cff.py nevts:1000'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('BprimeToBZinc_M_850_TuneZ2star_8TeV_madgraph_cff_py_GEN.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'START52_V9::All'

process.generator = cms.EDFilter("Pythia6HadronizerFilter",
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    comEnergy = cms.double(8000.0),
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
        processParameters = cms.vstring('MSTP(1) = 4', 
            'MSEL=7               ! User defined processes', 
            'MWID(7)=2            !use width of bprime as defined by PMAS', 
            'MSTJ(1)=1            ! Fragmentation/hadronization on or off', 
            'MSTP(61)=1           ! Parton showering on or off', 
            'PMAS(5,1)=4.8        ! b quark mass', 
            'PMAS(6,1)=172.5      ! t quark mass', 
            'PMAS(7,1) = 850.0D0     ! bprime quark mass', 
            'PMAS(7,2) = 8.500D0                       ', 
            'PMAS(7,3) = 85.00D0                ', 
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
            'MDME(56,1)=0     ! g b4 on/off switches for individual decay modes', 
            'MDME(57,1)=0     ! gamma b4', 
            'KFDP(58,2)=5     ! defines Z0 b', 
            'MDME(58,1)=1     ! Z0 b', 
            'MDME(59,1)=0     ! W u', 
            'MDME(60,1)=0     ! W c', 
            'MDME(61,1)=0     ! W t', 
            'MDME(62,1)=0     ! W t4', 
            'MDME(63,1)=0     ! h0 b4', 
            'MDME(64,1)=-1    ! H- c', 
            'MDME(65,1)=-1    ! H- t', 
            'BRAT(56)  = 0.0D0', 
            'BRAT(57)  = 0.0D0', 
            'BRAT(58)  = 1.0D0', 
            'BRAT(59)  = 0.0D0', 
            'BRAT(60)  = 0.0D0', 
            'BRAT(61)  = 0.0D0', 
            'BRAT(62)  = 0.0D0', 
            'BRAT(63)  = 0.0D0', 
            'BRAT(64)  = 0.0D0', 
            'BRAT(65)  = 0.0D0', 
            'MDME(174,1)=1     !Z decay into d dbar', 
            'MDME(175,1)=1     !Z decay into u ubar', 
            'MDME(176,1)=1     !Z decay into s sbar', 
            'MDME(177,1)=1     !Z decay into c cbar', 
            'MDME(178,1)=1     !Z decay into b bbar', 
            'MDME(179,1)=1     !Z decay into t tbar', 
            'MDME(180,1)=-1    !Z decay into b4 b4bar', 
            'MDME(181,1)=-1    !Z decay into t4 t4bar', 
            'MDME(182,1)=1     !Z decay into e- e+', 
            'MDME(183,1)=1     !Z decay into nu_e nu_ebar', 
            'MDME(184,1)=1     !Z decay into mu- mu+', 
            'MDME(185,1)=1     !Z decay into nu_mu nu_mubar', 
            'MDME(186,1)=1     !Z decay into tau- tau+', 
            'MDME(187,1)=1     !Z decay into nu_tau nu_taubar', 
            'MDME(188,1)=-1    !Z decay into tau4 tau4bar', 
            'MDME(189,1)=-1    !Z decay into nu_tau4 nu_tau4bar'),
        parameterSets = cms.vstring('pythiaUESettings', 
            'processParameters')
    )
)


process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 

