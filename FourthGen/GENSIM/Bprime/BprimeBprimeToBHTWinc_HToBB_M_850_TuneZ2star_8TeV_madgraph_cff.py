import FWCore.ParameterSet.Config as cms

#from Configuration.Generator.PythiaUEZ2Settings_cfi import *
from Configuration.Generator.PythiaUEZ2starSettings_cfi import *
generator = cms.EDFilter("Pythia6HadronizerFilter",
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    comEnergy = cms.double(8000.0),
    PythiaParameters = cms.PSet(
        pythiaUESettingsBlock,
        processParameters = cms.vstring(
	    'PMAS(25,1)=125.00D0    !mass of Higgs', 
            'MSTP(1) = 4',
            'MSEL=7         ! User defined processes',
	    'MWID(7)=2',
            'MSTJ(1)=1       ! Fragmentation/hadronization on or off',
            'MSTP(61)=1      ! Parton showering on or off',
            'PMAS(5,1)=4.8   ! b quark mass', #from Spring11 4000040
            'PMAS(6,1)=172.5 ! t quark mass', #from Spring11 4000040
            'PMAS(7,1) = 850.0D0     ! bprime quarks mass',
            'PMAS(7,2) = 8.500D0    ! bprime quark width',
            'PMAS(7,3) = 85.00D0 ! Max value above which the BW shape is truncated', 
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
            'BRAT(61)  = 0.5D0',
            'BRAT(62)  = 0.0D0',
            'BRAT(63)  = 0.5D0',
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
            'MDME(226,1)=0     !Higgs decay into W W',
            'MDME(190,1)=1     !W decay into u dbar',
            'MDME(191,1)=1     !W decay into c dbar',
            'MDME(192,1)=1     !W decay into t dbar',
            'MDME(193,1)=-1    !W decay into t4 dbar',
            'MDME(194,1)=1     !W decay into u sbar',
            'MDME(195,1)=1     !W decay into c sbar',
            'MDME(196,1)=1     !W decay into t sbar',
            'MDME(197,1)=-1    !W decay into t4 sbar',
            'MDME(198,1)=1     !W decay into u bbar',
            'MDME(199,1)=1     !W decay into c bbar',
            'MDME(200,1)=1     !W decay into t bbar',
            'MDME(201,1)=-1    !W decay into t4 bbar',
            'MDME(202,1)=-1    !W decay into u b4bar',
            'MDME(203,1)=-1    !W decay into c b4bar',
            'MDME(204,1)=-1    !W decay into t b4bar',
            'MDME(205,1)=-1    !W decay into t4 b4bar',
            'MDME(206,1)=1     !W decay into e- nu_e',
            'MDME(207,1)=1     !W decay into mu nu_mu',
            'MDME(208,1)=1     !W decay into tau nu_tau',
            'MDME(209,1)=-1    !W decay into tau4 nu_tau4'),
        # This is a vector of ParameterSet names to be read, in this order
        parameterSets = cms.vstring('pythiaUESettings',
            'processParameters')
    ),
    jetMatching = cms.untracked.PSet(
       scheme = cms.string("Madgraph"),
       mode = cms.string("auto"),       # soup, or "inclusive" / "exclusive"
       MEMAIN_etaclmax = cms.double(5.0),
       MEMAIN_qcut = cms.double(-1),
       MEMAIN_nqmatch = cms.int32(-1),
       MEMAIN_minjets = cms.int32(-1),
       MEMAIN_maxjets = cms.int32(-1),
       MEMAIN_showerkt = cms.double(0),
       MEMAIN_excres = cms.string(''),
       outTree_flag = cms.int32(0)
    )
)

ProductionFilterSequence = cms.Sequence(generator)

