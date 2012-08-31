# Auto generated configuration file
# using: 
# Revision: 1.372.2.3 
# Source: /local/reps/CMSSW.admin/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/EightTeV/Bprime_B2G/BprimeToBHinc_M_575_TuneZ2star_8TeV_madgraph_cff.py -s GEN,SIM --filetype=LHE --filein=6465 --beamspot Realistic8TeVCollision --conditions START52_V9::All --pileup NoPileUp --datatier GEN-SIM --eventcontent RAWSIM -n 30 --no_exe
import FWCore.ParameterSet.Config as cms

process = cms.Process('SIM')

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
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(30)
)

# Input source
process.source = cms.Source("LHESource",
    fileNames = cms.untracked.vstring('6465')
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.372.2.3 $'),
    annotation = cms.untracked.string('Configuration/GenProduction/python/EightTeV/Bprime_B2G/BprimeToBHinc_M_575_TuneZ2star_8TeV_madgraph_cff.py nevts:30'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('BprimeToBHinc_M_575_TuneZ2star_8TeV_madgraph_cff_py_GEN_SIM.root'),
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
        processParameters = cms.vstring('PMAS(25,1)=125.00D0    !mass of Higgs', 
            'MSTP(1) = 4', 
            'MSEL=7               ! User defined processes', 
            'MWID(7)=2            !use width of bprime as defined by PMAS', 
            'MSTJ(1)=1            ! Fragmentation/hadronization on or off', 
            'MSTP(61)=1           ! Parton showering on or off', 
            'PMAS(5,1)=4.8        ! b quark mass', 
            'PMAS(6,1)=172.5      ! t quark mass', 
            'PMAS(7,1) = 575.0D0     ! bprime quark mass', 
            'PMAS(7,2) = 5.750D0                       ', 
            'PMAS(7,3) = 57.50D0                ', 
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
            'MDME(58,1)=0     ! Z0 b', 
            'MDME(59,1)=0     ! W u', 
            'MDME(60,1)=0     ! W c', 
            'MDME(61,1)=1     ! W t', 
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
            'MDME(210,1)=1     !Higgs decay into dd', 
            'MDME(211,1)=1     !Higgs decay into uu', 
            'MDME(212,1)=1     !Higgs decay into ss', 
            'MDME(213,1)=1     !Higgs decay into cc', 
            'MDME(214,1)=1     !Higgs decay into bb', 
            'MDME(215,1)=1     !Higgs decay into tt', 
            'MDME(216,1)=1     !Higgs decay into', 
            'MDME(217,1)=1     !Higgs decay into Higgs decay', 
            'MDME(218,1)=1     !Higgs decay into e nu e', 
            'MDME(219,1)=1     !Higgs decay into mu nu mu', 
            'MDME(220,1)=1     !Higgs decay into tau nu tau', 
            'MDME(221,1)=1     !Higgs decay into Higgs decay', 
            'MDME(222,1)=1     !Higgs decay into g g', 
            'MDME(223,1)=1     !Higgs decay into gam gam', 
            'MDME(224,1)=1     !Higgs decay into gam Z', 
            'MDME(225,1)=1     !Higgs decay into Z Z', 
            'MDME(226,1)=1     !Higgs decay into W W'),
        parameterSets = cms.vstring('pythiaUESettings', 
            'processParameters')
    )
)


process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 

