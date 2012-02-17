//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct 27 06:05:14 2011 by ROOT version 5.27/06b
// from TTree root/root
// found on file: results_9_1_sMl.root
//////////////////////////////////////////////////////////

#ifndef bPrimeNtupleMC_h
#define bPrimeNtupleMC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class bPrimeNtupleMC {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           EvtInfo_RunNo;
   Long64_t        EvtInfo_EvtNo;
   Int_t           EvtInfo_BxNo;
   Int_t           EvtInfo_LumiNo;
   Int_t           EvtInfo_Orbit;
   Int_t           EvtInfo_McFlag;
   Int_t           EvtInfo_McSigTag;
   Int_t           EvtInfo_McbprimeMode[2];
   Int_t           EvtInfo_MctprimeMode[2];
   Int_t           EvtInfo_McWMode[4];
   Int_t           EvtInfo_McZMode[2];
   Float_t         EvtInfo_McbprimeMass[2];
   Float_t         EvtInfo_MctprimeMass[2];
   Float_t         EvtInfo_MctopMass[2];
   Float_t         EvtInfo_McWMass[4];
   Float_t         EvtInfo_McZMass[2];
   Float_t         EvtInfo_McDauPt[14];
   Float_t         EvtInfo_McDauEta[14];
   Float_t         EvtInfo_McDauPhi[14];
   Int_t           EvtInfo_McDauPdgID[14];
   Int_t           EvtInfo_PDFid1;
   Int_t           EvtInfo_PDFid2;
   Float_t         EvtInfo_PDFx1;
   Float_t         EvtInfo_RhoPU;
   Float_t         EvtInfo_SigmaPU;
   Float_t         EvtInfo_PDFx2;
   Float_t         EvtInfo_PDFscale;
   Float_t         EvtInfo_PDFv1;
   Float_t         EvtInfo_PDFv2;
   Float_t         EvtInfo_MET;
   Float_t         EvtInfo_METPhi;
   Float_t         EvtInfo_RawMET;
   Float_t         EvtInfo_RawMETPhi;
   Float_t         EvtInfo_SumEt;
   Float_t         EvtInfo_METSig;
   Float_t         EvtInfo_eLong;
   Float_t         EvtInfo_MaxHadTower;
   Float_t         EvtInfo_MaxEmTower;
   Float_t         EvtInfo_FracHad;
   Float_t         EvtInfo_FracEm;
   Float_t         EvtInfo_GenMET;
   Float_t         EvtInfo_GenMETPhi;
   Float_t         EvtInfo_PFMET;
   Float_t         EvtInfo_PFMETPhi;
   Float_t         EvtInfo_PFRawMET;
   Float_t         EvtInfo_PFRawMETPhi;
   Float_t         EvtInfo_PFSumEt;
   Float_t         EvtInfo_PFMETSig;
   Float_t         EvtInfo_PFGenMET;
   Float_t         EvtInfo_PFGenMETPhi;
   Float_t         EvtInfo_PFMETx;
   Float_t         EvtInfo_PFMETy;
   Int_t           EvtInfo_TrgCount;
   Int_t           EvtInfo_nTrgBook;
   Char_t          EvtInfo_TrgBook[3020];   //[EvtInfo.nTrgBook]
   Int_t           EvtInfo_L1[128];
   Int_t           EvtInfo_TT[64];
   Float_t         EvtInfo_HighPurityFraction;
   Int_t           EvtInfo_NofTracks;
   Int_t           EvtInfo_nHLT;
   Bool_t          EvtInfo_HLTbits[239];   //[EvtInfo.nHLT]
   Int_t           EvtInfo_nBX;
   Int_t           EvtInfo_nPU[3];   //[EvtInfo.nBX]
   Int_t           EvtInfo_BXPU[3];   //[EvtInfo.nBX]
   Int_t           GenInfo_Size;
   Float_t         GenInfo_Pt[13];   //[GenInfo.Size]
   Float_t         GenInfo_Eta[13];   //[GenInfo.Size]
   Float_t         GenInfo_Phi[13];   //[GenInfo.Size]
   Float_t         GenInfo_Mass[13];   //[GenInfo.Size]
   Int_t           GenInfo_PdgID[13];   //[GenInfo.Size]
   Int_t           GenInfo_Status[13];   //[GenInfo.Size]
   Int_t           GenInfo_nMo[13];   //[GenInfo.Size]
   Int_t           GenInfo_nDa[13];   //[GenInfo.Size]
   Int_t           GenInfo_Mo1[13];   //[GenInfo.Size]
   Int_t           GenInfo_Mo2[13];   //[GenInfo.Size]
   Int_t           GenInfo_Da1[13];   //[GenInfo.Size]
   Int_t           GenInfo_Da2[13];   //[GenInfo.Size]
   Int_t           PFLepInfo_Size;
   Int_t           PFLepInfo_Index[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_isEcalDriven[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_isTrackerDriven[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_LeptonType[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_Charge[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_ChargeGsf[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_ChargeCtf[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_ChargeScPix[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_Pt[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_Et[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_Eta[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_caloEta[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_Phi[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_TrackIso[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_EcalIso[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_HcalIso[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_HcalDepth1Iso[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_HcalDepth2Iso[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_ChargedHadronIso[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_NeutralHadronIso[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_PhotonIso[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_CaloEnergy[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_e1x5[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_e2x5Max[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_e5x5[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_Px[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_Py[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_Pz[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_Energy[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_vertexZ[5];   //[PFLepInfo.Size]
   Bool_t          PFLepInfo_MuIDGlobalMuonPromptTight[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_MuInnerPtError[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_MuGlobalPtError[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_MuInnerTrackDz[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_MuInnerTrackD0[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_MuInnerTrackDxy_BS[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_MuInnerTrackDxy_PV[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_MuInnerTrackDxy_PVBS[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_MuInnerTrackNHits[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_MuNTrackerHits[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_MuCaloCompat[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_MuNChambers[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_MuNChambersMatchesSegment[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_MuNPixelLayers[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_MuNPixelLayersWMeasurement[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_MuNLostInnerHits[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_MuNLostOuterHits[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_MuNMuonhits[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_MuType[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_MuGlobalNormalizedChi2[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_ElTrackNLostHits[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_ElTrackDz[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_ElTrackD0[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_ElTrackDxy_BS[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_ElTrackDxy_PV[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_ElTrackDxy_PVBS[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_simpleEleId95relIso[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_simpleEleId90relIso[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_simpleEleId85relIso[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_simpleEleId80relIso[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_simpleEleId70relIso[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_simpleEleId60relIso[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_simpleEleId95cIso[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_simpleEleId90cIso[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_simpleEleId85cIso[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_simpleEleId80cIso[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_simpleEleId70cIso[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_simpleEleId60cIso[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_eidVeryLoose[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_eidLoose[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_eidMedium[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_eidTight[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_eidSuperTight[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_eidHyperTight1[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_eidHyperTight2[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_eidHyperTight3[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_eidHyperTight4[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_eidVeryLooseMC[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_eidLooseMC[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_eidMediumMC[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_eidTightMC[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_eidSuperTightMC[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_eidHyperTight1MC[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_eidHyperTight2MC[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_eidHyperTight3MC[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_eidHyperTight4MC[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_ElEoverP[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_EldeltaEta[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_EldeltaPhi[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_ElHadoverEm[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_ElsigmaIetaIeta[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_ElscSigmaIetaIeta[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_ElEnergyErr[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_ElMomentumErr[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_ElTrackNHits[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_ElSharedHitsFraction[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_dR_gsf_ctfTrack[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_dPt_gsf_ctfTrack[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_ElNClusters[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_ElClassification[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_ElFBrem[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_ElNumberOfBrems[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_NumberOfExpectedInnerHits[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_Eldist[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_Eldcot[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_Elconvradius[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_ElConvPoint_x[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_ElConvPoint_y[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_ElConvPoint_z[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_dcotdist[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_ElseedEoverP[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_ElEcalIso04[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_ElHcalIso04[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_GenPt[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_GenEta[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_GenPhi[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_GenPdgID[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_GenMCTag[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_TrgPt[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_TrgEta[5];   //[PFLepInfo.Size]
   Float_t         PFLepInfo_TrgPhi[5];   //[PFLepInfo.Size]
   Int_t           PFLepInfo_TrgID[5];   //[PFLepInfo.Size]
   Int_t           LepInfo_Size;
   Int_t           LepInfo_Index[22];   //[LepInfo.Size]
   Int_t           LepInfo_isEcalDriven[22];   //[LepInfo.Size]
   Int_t           LepInfo_isTrackerDriven[22];   //[LepInfo.Size]
   Int_t           LepInfo_LeptonType[22];   //[LepInfo.Size]
   Int_t           LepInfo_Charge[22];   //[LepInfo.Size]
   Int_t           LepInfo_ChargeGsf[22];   //[LepInfo.Size]
   Int_t           LepInfo_ChargeCtf[22];   //[LepInfo.Size]
   Int_t           LepInfo_ChargeScPix[22];   //[LepInfo.Size]
   Float_t         LepInfo_Pt[22];   //[LepInfo.Size]
   Float_t         LepInfo_Et[22];   //[LepInfo.Size]
   Float_t         LepInfo_Eta[22];   //[LepInfo.Size]
   Float_t         LepInfo_caloEta[22];   //[LepInfo.Size]
   Float_t         LepInfo_Phi[22];   //[LepInfo.Size]
   Float_t         LepInfo_TrackIso[22];   //[LepInfo.Size]
   Float_t         LepInfo_EcalIso[22];   //[LepInfo.Size]
   Float_t         LepInfo_HcalIso[22];   //[LepInfo.Size]
   Float_t         LepInfo_HcalDepth1Iso[22];   //[LepInfo.Size]
   Float_t         LepInfo_HcalDepth2Iso[22];   //[LepInfo.Size]
   Float_t         LepInfo_ChargedHadronIso[22];   //[LepInfo.Size]
   Float_t         LepInfo_NeutralHadronIso[22];   //[LepInfo.Size]
   Float_t         LepInfo_PhotonIso[22];   //[LepInfo.Size]
   Float_t         LepInfo_CaloEnergy[22];   //[LepInfo.Size]
   Float_t         LepInfo_e1x5[22];   //[LepInfo.Size]
   Float_t         LepInfo_e2x5Max[22];   //[LepInfo.Size]
   Float_t         LepInfo_e5x5[22];   //[LepInfo.Size]
   Float_t         LepInfo_Px[22];   //[LepInfo.Size]
   Float_t         LepInfo_Py[22];   //[LepInfo.Size]
   Float_t         LepInfo_Pz[22];   //[LepInfo.Size]
   Float_t         LepInfo_Energy[22];   //[LepInfo.Size]
   Float_t         LepInfo_vertexZ[22];   //[LepInfo.Size]
   Bool_t          LepInfo_MuIDGlobalMuonPromptTight[22];   //[LepInfo.Size]
   Float_t         LepInfo_MuInnerPtError[22];   //[LepInfo.Size]
   Float_t         LepInfo_MuGlobalPtError[22];   //[LepInfo.Size]
   Float_t         LepInfo_MuInnerTrackDz[22];   //[LepInfo.Size]
   Float_t         LepInfo_MuInnerTrackD0[22];   //[LepInfo.Size]
   Float_t         LepInfo_MuInnerTrackDxy_BS[22];   //[LepInfo.Size]
   Float_t         LepInfo_MuInnerTrackDxy_PV[22];   //[LepInfo.Size]
   Float_t         LepInfo_MuInnerTrackDxy_PVBS[22];   //[LepInfo.Size]
   Int_t           LepInfo_MuInnerTrackNHits[22];   //[LepInfo.Size]
   Int_t           LepInfo_MuNTrackerHits[22];   //[LepInfo.Size]
   Float_t         LepInfo_MuCaloCompat[22];   //[LepInfo.Size]
   Int_t           LepInfo_MuNChambers[22];   //[LepInfo.Size]
   Int_t           LepInfo_MuNChambersMatchesSegment[22];   //[LepInfo.Size]
   Int_t           LepInfo_MuNPixelLayers[22];   //[LepInfo.Size]
   Int_t           LepInfo_MuNPixelLayersWMeasurement[22];   //[LepInfo.Size]
   Int_t           LepInfo_MuNLostInnerHits[22];   //[LepInfo.Size]
   Int_t           LepInfo_MuNLostOuterHits[22];   //[LepInfo.Size]
   Int_t           LepInfo_MuNMuonhits[22];   //[LepInfo.Size]
   Int_t           LepInfo_MuType[22];   //[LepInfo.Size]
   Float_t         LepInfo_MuGlobalNormalizedChi2[22];   //[LepInfo.Size]
   Float_t         LepInfo_ElTrackNLostHits[22];   //[LepInfo.Size]
   Float_t         LepInfo_ElTrackDz[22];   //[LepInfo.Size]
   Float_t         LepInfo_ElTrackD0[22];   //[LepInfo.Size]
   Float_t         LepInfo_ElTrackDxy_BS[22];   //[LepInfo.Size]
   Float_t         LepInfo_ElTrackDxy_PV[22];   //[LepInfo.Size]
   Float_t         LepInfo_ElTrackDxy_PVBS[22];   //[LepInfo.Size]
   Float_t         LepInfo_simpleEleId95relIso[22];   //[LepInfo.Size]
   Float_t         LepInfo_simpleEleId90relIso[22];   //[LepInfo.Size]
   Float_t         LepInfo_simpleEleId85relIso[22];   //[LepInfo.Size]
   Float_t         LepInfo_simpleEleId80relIso[22];   //[LepInfo.Size]
   Float_t         LepInfo_simpleEleId70relIso[22];   //[LepInfo.Size]
   Float_t         LepInfo_simpleEleId60relIso[22];   //[LepInfo.Size]
   Float_t         LepInfo_simpleEleId95cIso[22];   //[LepInfo.Size]
   Float_t         LepInfo_simpleEleId90cIso[22];   //[LepInfo.Size]
   Float_t         LepInfo_simpleEleId85cIso[22];   //[LepInfo.Size]
   Float_t         LepInfo_simpleEleId80cIso[22];   //[LepInfo.Size]
   Float_t         LepInfo_simpleEleId70cIso[22];   //[LepInfo.Size]
   Float_t         LepInfo_simpleEleId60cIso[22];   //[LepInfo.Size]
   Float_t         LepInfo_eidVeryLoose[22];   //[LepInfo.Size]
   Float_t         LepInfo_eidLoose[22];   //[LepInfo.Size]
   Float_t         LepInfo_eidMedium[22];   //[LepInfo.Size]
   Float_t         LepInfo_eidTight[22];   //[LepInfo.Size]
   Float_t         LepInfo_eidSuperTight[22];   //[LepInfo.Size]
   Float_t         LepInfo_eidHyperTight1[22];   //[LepInfo.Size]
   Float_t         LepInfo_eidHyperTight2[22];   //[LepInfo.Size]
   Float_t         LepInfo_eidHyperTight3[22];   //[LepInfo.Size]
   Float_t         LepInfo_eidHyperTight4[22];   //[LepInfo.Size]
   Float_t         LepInfo_eidVeryLooseMC[22];   //[LepInfo.Size]
   Float_t         LepInfo_eidLooseMC[22];   //[LepInfo.Size]
   Float_t         LepInfo_eidMediumMC[22];   //[LepInfo.Size]
   Float_t         LepInfo_eidTightMC[22];   //[LepInfo.Size]
   Float_t         LepInfo_eidSuperTightMC[22];   //[LepInfo.Size]
   Float_t         LepInfo_eidHyperTight1MC[22];   //[LepInfo.Size]
   Float_t         LepInfo_eidHyperTight2MC[22];   //[LepInfo.Size]
   Float_t         LepInfo_eidHyperTight3MC[22];   //[LepInfo.Size]
   Float_t         LepInfo_eidHyperTight4MC[22];   //[LepInfo.Size]
   Float_t         LepInfo_ElEoverP[22];   //[LepInfo.Size]
   Float_t         LepInfo_EldeltaEta[22];   //[LepInfo.Size]
   Float_t         LepInfo_EldeltaPhi[22];   //[LepInfo.Size]
   Float_t         LepInfo_ElHadoverEm[22];   //[LepInfo.Size]
   Float_t         LepInfo_ElsigmaIetaIeta[22];   //[LepInfo.Size]
   Float_t         LepInfo_ElscSigmaIetaIeta[22];   //[LepInfo.Size]
   Float_t         LepInfo_ElEnergyErr[22];   //[LepInfo.Size]
   Float_t         LepInfo_ElMomentumErr[22];   //[LepInfo.Size]
   Int_t           LepInfo_ElTrackNHits[22];   //[LepInfo.Size]
   Float_t         LepInfo_ElSharedHitsFraction[22];   //[LepInfo.Size]
   Float_t         LepInfo_dR_gsf_ctfTrack[22];   //[LepInfo.Size]
   Float_t         LepInfo_dPt_gsf_ctfTrack[22];   //[LepInfo.Size]
   Int_t           LepInfo_ElNClusters[22];   //[LepInfo.Size]
   Int_t           LepInfo_ElClassification[22];   //[LepInfo.Size]
   Float_t         LepInfo_ElFBrem[22];   //[LepInfo.Size]
   Int_t           LepInfo_ElNumberOfBrems[22];   //[LepInfo.Size]
   Int_t           LepInfo_NumberOfExpectedInnerHits[22];   //[LepInfo.Size]
   Float_t         LepInfo_Eldist[22];   //[LepInfo.Size]
   Float_t         LepInfo_Eldcot[22];   //[LepInfo.Size]
   Float_t         LepInfo_Elconvradius[22];   //[LepInfo.Size]
   Float_t         LepInfo_ElConvPoint_x[22];   //[LepInfo.Size]
   Float_t         LepInfo_ElConvPoint_y[22];   //[LepInfo.Size]
   Float_t         LepInfo_ElConvPoint_z[22];   //[LepInfo.Size]
   Float_t         LepInfo_dcotdist[22];   //[LepInfo.Size]
   Float_t         LepInfo_ElseedEoverP[22];   //[LepInfo.Size]
   Float_t         LepInfo_ElEcalIso04[22];   //[LepInfo.Size]
   Float_t         LepInfo_ElHcalIso04[22];   //[LepInfo.Size]
   Float_t         LepInfo_GenPt[22];   //[LepInfo.Size]
   Float_t         LepInfo_GenEta[22];   //[LepInfo.Size]
   Float_t         LepInfo_GenPhi[22];   //[LepInfo.Size]
   Int_t           LepInfo_GenPdgID[22];   //[LepInfo.Size]
   Int_t           LepInfo_GenMCTag[22];   //[LepInfo.Size]
   Float_t         LepInfo_TrgPt[22];   //[LepInfo.Size]
   Float_t         LepInfo_TrgEta[22];   //[LepInfo.Size]
   Float_t         LepInfo_TrgPhi[22];   //[LepInfo.Size]
   Int_t           LepInfo_TrgID[22];   //[LepInfo.Size]
   Int_t           VertexInfo_Size;
   Int_t           VertexInfo_isValid[53];   //[VertexInfo.Size]
   Bool_t          VertexInfo_isFake[53];   //[VertexInfo.Size]
   Int_t           VertexInfo_Type[53];   //[VertexInfo.Size]
   Float_t         VertexInfo_Ndof[53];   //[VertexInfo.Size]
   Float_t         VertexInfo_NormalizedChi2[53];   //[VertexInfo.Size]
   Float_t         VertexInfo_Pt_Sum[53];   //[VertexInfo.Size]
   Float_t         VertexInfo_x[53];   //[VertexInfo.Size]
   Float_t         VertexInfo_y[53];   //[VertexInfo.Size]
   Float_t         VertexInfo_z[53];   //[VertexInfo.Size]
   Float_t         VertexInfo_Rho[53];   //[VertexInfo.Size]
   Int_t           PFJetInfo_Size;
   Int_t           PFJetInfo_Index[12];   //[PFJetInfo.Size]
   Int_t           PFJetInfo_NTracks[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_Et[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_Pt[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_Unc[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_Eta[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_Phi[12];   //[PFJetInfo.Size]
   Int_t           PFJetInfo_JetIDLOOSE[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_JetCharge[12];   //[PFJetInfo.Size]
   Int_t           PFJetInfo_NConstituents[12];   //[PFJetInfo.Size]
   Int_t           PFJetInfo_NCH[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_CEF[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_NHF[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_NEF[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_CHF[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_JVAlpha[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_JVBeta[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_PtCorrRaw[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_PtCorrL2[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_PtCorrL3[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_PtCorrL7g[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_PtCorrL7uds[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_PtCorrL7c[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_PtCorrL7b[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_JetBProbBJetTags[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_JetProbBJetTags[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_TrackCountHiPurBJetTags[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_TrackCountHiEffBJetTags[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_SimpleSVBJetTags[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_SimpleSVHEBJetTags[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_SimpleSVHPBJetTags[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_CombinedSVBJetTags[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_CombinedSVMVABJetTags[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_SoftElecByIP3dBJetTags[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_SoftElecByPtBJetTags[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_SoftMuonBJetTags[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_SoftMuonByIP3dBJetTags[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_SoftMuonByPtBJetTags[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_GenJetPt[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_GenJetEta[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_GenJetPhi[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_GenPt[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_GenEta[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_GenPhi[12];   //[PFJetInfo.Size]
   Int_t           PFJetInfo_GenPdgID[12];   //[PFJetInfo.Size]
   Int_t           PFJetInfo_GenFlavor[12];   //[PFJetInfo.Size]
   Int_t           PFJetInfo_GenMCTag[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_Px[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_Py[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_Pz[12];   //[PFJetInfo.Size]
   Float_t         PFJetInfo_Energy[12];   //[PFJetInfo.Size]
   Int_t           JetInfo_Size;
   Int_t           JetInfo_Index[12];   //[JetInfo.Size]
   Int_t           JetInfo_NTracks[12];   //[JetInfo.Size]
   Float_t         JetInfo_Et[12];   //[JetInfo.Size]
   Float_t         JetInfo_Pt[12];   //[JetInfo.Size]
   Float_t         JetInfo_Unc[12];   //[JetInfo.Size]
   Float_t         JetInfo_Eta[12];   //[JetInfo.Size]
   Float_t         JetInfo_Phi[12];   //[JetInfo.Size]
   Int_t           JetInfo_JetIDLOOSE[12];   //[JetInfo.Size]
   Float_t         JetInfo_JetCharge[12];   //[JetInfo.Size]
   Int_t           JetInfo_NConstituents[12];   //[JetInfo.Size]
   Int_t           JetInfo_NCH[12];   //[JetInfo.Size]
   Float_t         JetInfo_CEF[12];   //[JetInfo.Size]
   Float_t         JetInfo_NHF[12];   //[JetInfo.Size]
   Float_t         JetInfo_NEF[12];   //[JetInfo.Size]
   Float_t         JetInfo_CHF[12];   //[JetInfo.Size]
   Float_t         JetInfo_JVAlpha[12];   //[JetInfo.Size]
   Float_t         JetInfo_JVBeta[12];   //[JetInfo.Size]
   Float_t         JetInfo_PtCorrRaw[12];   //[JetInfo.Size]
   Float_t         JetInfo_PtCorrL2[12];   //[JetInfo.Size]
   Float_t         JetInfo_PtCorrL3[12];   //[JetInfo.Size]
   Float_t         JetInfo_PtCorrL7g[12];   //[JetInfo.Size]
   Float_t         JetInfo_PtCorrL7uds[12];   //[JetInfo.Size]
   Float_t         JetInfo_PtCorrL7c[12];   //[JetInfo.Size]
   Float_t         JetInfo_PtCorrL7b[12];   //[JetInfo.Size]
   Float_t         JetInfo_JetBProbBJetTags[12];   //[JetInfo.Size]
   Float_t         JetInfo_JetProbBJetTags[12];   //[JetInfo.Size]
   Float_t         JetInfo_TrackCountHiPurBJetTags[12];   //[JetInfo.Size]
   Float_t         JetInfo_TrackCountHiEffBJetTags[12];   //[JetInfo.Size]
   Float_t         JetInfo_SimpleSVBJetTags[12];   //[JetInfo.Size]
   Float_t         JetInfo_SimpleSVHEBJetTags[12];   //[JetInfo.Size]
   Float_t         JetInfo_SimpleSVHPBJetTags[12];   //[JetInfo.Size]
   Float_t         JetInfo_CombinedSVBJetTags[12];   //[JetInfo.Size]
   Float_t         JetInfo_CombinedSVMVABJetTags[12];   //[JetInfo.Size]
   Float_t         JetInfo_SoftElecByIP3dBJetTags[12];   //[JetInfo.Size]
   Float_t         JetInfo_SoftElecByPtBJetTags[12];   //[JetInfo.Size]
   Float_t         JetInfo_SoftMuonBJetTags[12];   //[JetInfo.Size]
   Float_t         JetInfo_SoftMuonByIP3dBJetTags[12];   //[JetInfo.Size]
   Float_t         JetInfo_SoftMuonByPtBJetTags[12];   //[JetInfo.Size]
   Float_t         JetInfo_GenJetPt[12];   //[JetInfo.Size]
   Float_t         JetInfo_GenJetEta[12];   //[JetInfo.Size]
   Float_t         JetInfo_GenJetPhi[12];   //[JetInfo.Size]
   Float_t         JetInfo_GenPt[12];   //[JetInfo.Size]
   Float_t         JetInfo_GenEta[12];   //[JetInfo.Size]
   Float_t         JetInfo_GenPhi[12];   //[JetInfo.Size]
   Int_t           JetInfo_GenPdgID[12];   //[JetInfo.Size]
   Int_t           JetInfo_GenFlavor[12];   //[JetInfo.Size]
   Int_t           JetInfo_GenMCTag[12];   //[JetInfo.Size]
   Float_t         JetInfo_Px[12];   //[JetInfo.Size]
   Float_t         JetInfo_Py[12];   //[JetInfo.Size]
   Float_t         JetInfo_Pz[12];   //[JetInfo.Size]
   Float_t         JetInfo_Energy[12];   //[JetInfo.Size]
   Int_t           PairInfo_Size;
   Int_t           PairInfo_Index[231];   //[PairInfo.Size]
   Int_t           PairInfo_Type[231];   //[PairInfo.Size]
   Int_t           PairInfo_Obj1Index[231];   //[PairInfo.Size]
   Int_t           PairInfo_Obj2Index[231];   //[PairInfo.Size]
   Float_t         PairInfo_Mass[231];   //[PairInfo.Size]
   Float_t         PairInfo_Pt[231];   //[PairInfo.Size]
   Float_t         PairInfo_Eta[231];   //[PairInfo.Size]
   Float_t         PairInfo_Phi[231];   //[PairInfo.Size]
   Float_t         PairInfo_GenMass[231];   //[PairInfo.Size]
   Float_t         PairInfo_GenPt[231];   //[PairInfo.Size]
   Float_t         PairInfo_GenEta[231];   //[PairInfo.Size]
   Float_t         PairInfo_GenPhi[231];   //[PairInfo.Size]
   Int_t           PairInfo_GenPdgID[231];   //[PairInfo.Size]

   // List of branches
   TBranch        *b_EvtInfo_RunNo;   //!
   TBranch        *b_EvtInfo_EvtNo;   //!
   TBranch        *b_EvtInfo_BxNo;   //!
   TBranch        *b_EvtInfo_LumiNo;   //!
   TBranch        *b_EvtInfo_Orbit;   //!
   TBranch        *b_EvtInfo_McFlag;   //!
   TBranch        *b_EvtInfo_McSigTag;   //!
   TBranch        *b_EvtInfo_McbprimeMode;   //!
   TBranch        *b_EvtInfo_MctprimeMode;   //!
   TBranch        *b_EvtInfo_McWMode;   //!
   TBranch        *b_EvtInfo_McZMode;   //!
   TBranch        *b_EvtInfo_McbprimeMass;   //!
   TBranch        *b_EvtInfo_MctprimeMass;   //!
   TBranch        *b_EvtInfo_MctopMass;   //!
   TBranch        *b_EvtInfo_McWMass;   //!
   TBranch        *b_EvtInfo_McZMass;   //!
   TBranch        *b_EvtInfo_McDauPt;   //!
   TBranch        *b_EvtInfo_McDauEta;   //!
   TBranch        *b_EvtInfo_McDauPhi;   //!
   TBranch        *b_EvtInfo_McDauPdgID;   //!
   TBranch        *b_EvtInfo_PDFid1;   //!
   TBranch        *b_EvtInfo_PDFid2;   //!
   TBranch        *b_EvtInfo_PDFx1;   //!
   TBranch        *b_EvtInfo_RhoPU;   //!
   TBranch        *b_EvtInfo_SigmaPU;   //!
   TBranch        *b_EvtInfo_PDFx2;   //!
   TBranch        *b_EvtInfo_PDFscale;   //!
   TBranch        *b_EvtInfo_PDFv1;   //!
   TBranch        *b_EvtInfo_PDFv2;   //!
   TBranch        *b_EvtInfo_MET;   //!
   TBranch        *b_EvtInfo_METPhi;   //!
   TBranch        *b_EvtInfo_RawMET;   //!
   TBranch        *b_EvtInfo_RawMETPhi;   //!
   TBranch        *b_EvtInfo_SumEt;   //!
   TBranch        *b_EvtInfo_METSig;   //!
   TBranch        *b_EvtInfo_eLong;   //!
   TBranch        *b_EvtInfo_MaxHadTower;   //!
   TBranch        *b_EvtInfo_MaxEmTower;   //!
   TBranch        *b_EvtInfo_FracHad;   //!
   TBranch        *b_EvtInfo_FracEm;   //!
   TBranch        *b_EvtInfo_GenMET;   //!
   TBranch        *b_EvtInfo_GenMETPhi;   //!
   TBranch        *b_EvtInfo_PFMET;   //!
   TBranch        *b_EvtInfo_PFMETPhi;   //!
   TBranch        *b_EvtInfo_PFRawMET;   //!
   TBranch        *b_EvtInfo_PFRawMETPhi;   //!
   TBranch        *b_EvtInfo_PFSumEt;   //!
   TBranch        *b_EvtInfo_PFMETSig;   //!
   TBranch        *b_EvtInfo_PFGenMET;   //!
   TBranch        *b_EvtInfo_PFGenMETPhi;   //!
   TBranch        *b_EvtInfo_PFMETx;   //!
   TBranch        *b_EvtInfo_PFMETy;   //!
   TBranch        *b_EvtInfo_TrgCount;   //!
   TBranch        *b_EvtInfo_nTrgBook;   //!
   TBranch        *b_EvtInfo_TrgBook;   //!
   TBranch        *b_EvtInfo_L1;   //!
   TBranch        *b_EvtInfo_TT;   //!
   TBranch        *b_EvtInfo_HighPurityFraction;   //!
   TBranch        *b_EvtInfo_NofTracks;   //!
   TBranch        *b_EvtInfo_nHLT;   //!
   TBranch        *b_EvtInfo_HLTbits;   //!
   TBranch        *b_EvtInfo_nBX;   //!
   TBranch        *b_EvtInfo_nPU;   //!
   TBranch        *b_EvtInfo_BXPU;   //!
   TBranch        *b_GenInfo_Size;   //!
   TBranch        *b_GenInfo_Pt;   //!
   TBranch        *b_GenInfo_Eta;   //!
   TBranch        *b_GenInfo_Phi;   //!
   TBranch        *b_GenInfo_Mass;   //!
   TBranch        *b_GenInfo_PdgID;   //!
   TBranch        *b_GenInfo_Status;   //!
   TBranch        *b_GenInfo_nMo;   //!
   TBranch        *b_GenInfo_nDa;   //!
   TBranch        *b_GenInfo_Mo1;   //!
   TBranch        *b_GenInfo_Mo2;   //!
   TBranch        *b_GenInfo_Da1;   //!
   TBranch        *b_GenInfo_Da2;   //!
   TBranch        *b_PFLepInfo_Size;   //!
   TBranch        *b_PFLepInfo_Index;   //!
   TBranch        *b_PFLepInfo_isEcalDriven;   //!
   TBranch        *b_PFLepInfo_isTrackerDriven;   //!
   TBranch        *b_PFLepInfo_LeptonType;   //!
   TBranch        *b_PFLepInfo_Charge;   //!
   TBranch        *b_PFLepInfo_ChargeGsf;   //!
   TBranch        *b_PFLepInfo_ChargeCtf;   //!
   TBranch        *b_PFLepInfo_ChargeScPix;   //!
   TBranch        *b_PFLepInfo_Pt;   //!
   TBranch        *b_PFLepInfo_Et;   //!
   TBranch        *b_PFLepInfo_Eta;   //!
   TBranch        *b_PFLepInfo_caloEta;   //!
   TBranch        *b_PFLepInfo_Phi;   //!
   TBranch        *b_PFLepInfo_TrackIso;   //!
   TBranch        *b_PFLepInfo_EcalIso;   //!
   TBranch        *b_PFLepInfo_HcalIso;   //!
   TBranch        *b_PFLepInfo_HcalDepth1Iso;   //!
   TBranch        *b_PFLepInfo_HcalDepth2Iso;   //!
   TBranch        *b_PFLepInfo_ChargedHadronIso;   //!
   TBranch        *b_PFLepInfo_NeutralHadronIso;   //!
   TBranch        *b_PFLepInfo_PhotonIso;   //!
   TBranch        *b_PFLepInfo_CaloEnergy;   //!
   TBranch        *b_PFLepInfo_e1x5;   //!
   TBranch        *b_PFLepInfo_e2x5Max;   //!
   TBranch        *b_PFLepInfo_e5x5;   //!
   TBranch        *b_PFLepInfo_Px;   //!
   TBranch        *b_PFLepInfo_Py;   //!
   TBranch        *b_PFLepInfo_Pz;   //!
   TBranch        *b_PFLepInfo_Energy;   //!
   TBranch        *b_PFLepInfo_vertexZ;   //!
   TBranch        *b_PFLepInfo_MuIDGlobalMuonPromptTight;   //!
   TBranch        *b_PFLepInfo_MuInnerPtError;   //!
   TBranch        *b_PFLepInfo_MuGlobalPtError;   //!
   TBranch        *b_PFLepInfo_MuInnerTrackDz;   //!
   TBranch        *b_PFLepInfo_MuInnerTrackD0;   //!
   TBranch        *b_PFLepInfo_MuInnerTrackDxy_BS;   //!
   TBranch        *b_PFLepInfo_MuInnerTrackDxy_PV;   //!
   TBranch        *b_PFLepInfo_MuInnerTrackDxy_PVBS;   //!
   TBranch        *b_PFLepInfo_MuInnerTrackNHits;   //!
   TBranch        *b_PFLepInfo_MuNTrackerHits;   //!
   TBranch        *b_PFLepInfo_MuCaloCompat;   //!
   TBranch        *b_PFLepInfo_MuNChambers;   //!
   TBranch        *b_PFLepInfo_MuNChambersMatchesSegment;   //!
   TBranch        *b_PFLepInfo_MuNPixelLayers;   //!
   TBranch        *b_PFLepInfo_MuNPixelLayersWMeasurement;   //!
   TBranch        *b_PFLepInfo_MuNLostInnerHits;   //!
   TBranch        *b_PFLepInfo_MuNLostOuterHits;   //!
   TBranch        *b_PFLepInfo_MuNMuonhits;   //!
   TBranch        *b_PFLepInfo_MuType;   //!
   TBranch        *b_PFLepInfo_MuGlobalNormalizedChi2;   //!
   TBranch        *b_PFLepInfo_ElTrackNLostHits;   //!
   TBranch        *b_PFLepInfo_ElTrackDz;   //!
   TBranch        *b_PFLepInfo_ElTrackD0;   //!
   TBranch        *b_PFLepInfo_ElTrackDxy_BS;   //!
   TBranch        *b_PFLepInfo_ElTrackDxy_PV;   //!
   TBranch        *b_PFLepInfo_ElTrackDxy_PVBS;   //!
   TBranch        *b_PFLepInfo_simpleEleId95relIso;   //!
   TBranch        *b_PFLepInfo_simpleEleId90relIso;   //!
   TBranch        *b_PFLepInfo_simpleEleId85relIso;   //!
   TBranch        *b_PFLepInfo_simpleEleId80relIso;   //!
   TBranch        *b_PFLepInfo_simpleEleId70relIso;   //!
   TBranch        *b_PFLepInfo_simpleEleId60relIso;   //!
   TBranch        *b_PFLepInfo_simpleEleId95cIso;   //!
   TBranch        *b_PFLepInfo_simpleEleId90cIso;   //!
   TBranch        *b_PFLepInfo_simpleEleId85cIso;   //!
   TBranch        *b_PFLepInfo_simpleEleId80cIso;   //!
   TBranch        *b_PFLepInfo_simpleEleId70cIso;   //!
   TBranch        *b_PFLepInfo_simpleEleId60cIso;   //!
   TBranch        *b_PFLepInfo_eidVeryLoose;   //!
   TBranch        *b_PFLepInfo_eidLoose;   //!
   TBranch        *b_PFLepInfo_eidMedium;   //!
   TBranch        *b_PFLepInfo_eidTight;   //!
   TBranch        *b_PFLepInfo_eidSuperTight;   //!
   TBranch        *b_PFLepInfo_eidHyperTight1;   //!
   TBranch        *b_PFLepInfo_eidHyperTight2;   //!
   TBranch        *b_PFLepInfo_eidHyperTight3;   //!
   TBranch        *b_PFLepInfo_eidHyperTight4;   //!
   TBranch        *b_PFLepInfo_eidVeryLooseMC;   //!
   TBranch        *b_PFLepInfo_eidLooseMC;   //!
   TBranch        *b_PFLepInfo_eidMediumMC;   //!
   TBranch        *b_PFLepInfo_eidTightMC;   //!
   TBranch        *b_PFLepInfo_eidSuperTightMC;   //!
   TBranch        *b_PFLepInfo_eidHyperTight1MC;   //!
   TBranch        *b_PFLepInfo_eidHyperTight2MC;   //!
   TBranch        *b_PFLepInfo_eidHyperTight3MC;   //!
   TBranch        *b_PFLepInfo_eidHyperTight4MC;   //!
   TBranch        *b_PFLepInfo_ElEoverP;   //!
   TBranch        *b_PFLepInfo_EldeltaEta;   //!
   TBranch        *b_PFLepInfo_EldeltaPhi;   //!
   TBranch        *b_PFLepInfo_ElHadoverEm;   //!
   TBranch        *b_PFLepInfo_ElsigmaIetaIeta;   //!
   TBranch        *b_PFLepInfo_ElscSigmaIetaIeta;   //!
   TBranch        *b_PFLepInfo_ElEnergyErr;   //!
   TBranch        *b_PFLepInfo_ElMomentumErr;   //!
   TBranch        *b_PFLepInfo_ElTrackNHits;   //!
   TBranch        *b_PFLepInfo_ElSharedHitsFraction;   //!
   TBranch        *b_PFLepInfo_dR_gsf_ctfTrack;   //!
   TBranch        *b_PFLepInfo_dPt_gsf_ctfTrack;   //!
   TBranch        *b_PFLepInfo_ElNClusters;   //!
   TBranch        *b_PFLepInfo_ElClassification;   //!
   TBranch        *b_PFLepInfo_ElFBrem;   //!
   TBranch        *b_PFLepInfo_ElNumberOfBrems;   //!
   TBranch        *b_PFLepInfo_NumberOfExpectedInnerHits;   //!
   TBranch        *b_PFLepInfo_Eldist;   //!
   TBranch        *b_PFLepInfo_Eldcot;   //!
   TBranch        *b_PFLepInfo_Elconvradius;   //!
   TBranch        *b_PFLepInfo_ElConvPoint_x;   //!
   TBranch        *b_PFLepInfo_ElConvPoint_y;   //!
   TBranch        *b_PFLepInfo_ElConvPoint_z;   //!
   TBranch        *b_PFLepInfo_dcotdist;   //!
   TBranch        *b_PFLepInfo_ElseedEoverP;   //!
   TBranch        *b_PFLepInfo_ElEcalIso04;   //!
   TBranch        *b_PFLepInfo_ElHcalIso04;   //!
   TBranch        *b_PFLepInfo_GenPt;   //!
   TBranch        *b_PFLepInfo_GenEta;   //!
   TBranch        *b_PFLepInfo_GenPhi;   //!
   TBranch        *b_PFLepInfo_GenPdgID;   //!
   TBranch        *b_PFLepInfo_GenMCTag;   //!
   TBranch        *b_PFLepInfo_TrgPt;   //!
   TBranch        *b_PFLepInfo_TrgEta;   //!
   TBranch        *b_PFLepInfo_TrgPhi;   //!
   TBranch        *b_PFLepInfo_TrgID;   //!
   TBranch        *b_LepInfo_Size;   //!
   TBranch        *b_LepInfo_Index;   //!
   TBranch        *b_LepInfo_isEcalDriven;   //!
   TBranch        *b_LepInfo_isTrackerDriven;   //!
   TBranch        *b_LepInfo_LeptonType;   //!
   TBranch        *b_LepInfo_Charge;   //!
   TBranch        *b_LepInfo_ChargeGsf;   //!
   TBranch        *b_LepInfo_ChargeCtf;   //!
   TBranch        *b_LepInfo_ChargeScPix;   //!
   TBranch        *b_LepInfo_Pt;   //!
   TBranch        *b_LepInfo_Et;   //!
   TBranch        *b_LepInfo_Eta;   //!
   TBranch        *b_LepInfo_caloEta;   //!
   TBranch        *b_LepInfo_Phi;   //!
   TBranch        *b_LepInfo_TrackIso;   //!
   TBranch        *b_LepInfo_EcalIso;   //!
   TBranch        *b_LepInfo_HcalIso;   //!
   TBranch        *b_LepInfo_HcalDepth1Iso;   //!
   TBranch        *b_LepInfo_HcalDepth2Iso;   //!
   TBranch        *b_LepInfo_ChargedHadronIso;   //!
   TBranch        *b_LepInfo_NeutralHadronIso;   //!
   TBranch        *b_LepInfo_PhotonIso;   //!
   TBranch        *b_LepInfo_CaloEnergy;   //!
   TBranch        *b_LepInfo_e1x5;   //!
   TBranch        *b_LepInfo_e2x5Max;   //!
   TBranch        *b_LepInfo_e5x5;   //!
   TBranch        *b_LepInfo_Px;   //!
   TBranch        *b_LepInfo_Py;   //!
   TBranch        *b_LepInfo_Pz;   //!
   TBranch        *b_LepInfo_Energy;   //!
   TBranch        *b_LepInfo_vertexZ;   //!
   TBranch        *b_LepInfo_MuIDGlobalMuonPromptTight;   //!
   TBranch        *b_LepInfo_MuInnerPtError;   //!
   TBranch        *b_LepInfo_MuGlobalPtError;   //!
   TBranch        *b_LepInfo_MuInnerTrackDz;   //!
   TBranch        *b_LepInfo_MuInnerTrackD0;   //!
   TBranch        *b_LepInfo_MuInnerTrackDxy_BS;   //!
   TBranch        *b_LepInfo_MuInnerTrackDxy_PV;   //!
   TBranch        *b_LepInfo_MuInnerTrackDxy_PVBS;   //!
   TBranch        *b_LepInfo_MuInnerTrackNHits;   //!
   TBranch        *b_LepInfo_MuNTrackerHits;   //!
   TBranch        *b_LepInfo_MuCaloCompat;   //!
   TBranch        *b_LepInfo_MuNChambers;   //!
   TBranch        *b_LepInfo_MuNChambersMatchesSegment;   //!
   TBranch        *b_LepInfo_MuNPixelLayers;   //!
   TBranch        *b_LepInfo_MuNPixelLayersWMeasurement;   //!
   TBranch        *b_LepInfo_MuNLostInnerHits;   //!
   TBranch        *b_LepInfo_MuNLostOuterHits;   //!
   TBranch        *b_LepInfo_MuNMuonhits;   //!
   TBranch        *b_LepInfo_MuType;   //!
   TBranch        *b_LepInfo_MuGlobalNormalizedChi2;   //!
   TBranch        *b_LepInfo_ElTrackNLostHits;   //!
   TBranch        *b_LepInfo_ElTrackDz;   //!
   TBranch        *b_LepInfo_ElTrackD0;   //!
   TBranch        *b_LepInfo_ElTrackDxy_BS;   //!
   TBranch        *b_LepInfo_ElTrackDxy_PV;   //!
   TBranch        *b_LepInfo_ElTrackDxy_PVBS;   //!
   TBranch        *b_LepInfo_simpleEleId95relIso;   //!
   TBranch        *b_LepInfo_simpleEleId90relIso;   //!
   TBranch        *b_LepInfo_simpleEleId85relIso;   //!
   TBranch        *b_LepInfo_simpleEleId80relIso;   //!
   TBranch        *b_LepInfo_simpleEleId70relIso;   //!
   TBranch        *b_LepInfo_simpleEleId60relIso;   //!
   TBranch        *b_LepInfo_simpleEleId95cIso;   //!
   TBranch        *b_LepInfo_simpleEleId90cIso;   //!
   TBranch        *b_LepInfo_simpleEleId85cIso;   //!
   TBranch        *b_LepInfo_simpleEleId80cIso;   //!
   TBranch        *b_LepInfo_simpleEleId70cIso;   //!
   TBranch        *b_LepInfo_simpleEleId60cIso;   //!
   TBranch        *b_LepInfo_eidVeryLoose;   //!
   TBranch        *b_LepInfo_eidLoose;   //!
   TBranch        *b_LepInfo_eidMedium;   //!
   TBranch        *b_LepInfo_eidTight;   //!
   TBranch        *b_LepInfo_eidSuperTight;   //!
   TBranch        *b_LepInfo_eidHyperTight1;   //!
   TBranch        *b_LepInfo_eidHyperTight2;   //!
   TBranch        *b_LepInfo_eidHyperTight3;   //!
   TBranch        *b_LepInfo_eidHyperTight4;   //!
   TBranch        *b_LepInfo_eidVeryLooseMC;   //!
   TBranch        *b_LepInfo_eidLooseMC;   //!
   TBranch        *b_LepInfo_eidMediumMC;   //!
   TBranch        *b_LepInfo_eidTightMC;   //!
   TBranch        *b_LepInfo_eidSuperTightMC;   //!
   TBranch        *b_LepInfo_eidHyperTight1MC;   //!
   TBranch        *b_LepInfo_eidHyperTight2MC;   //!
   TBranch        *b_LepInfo_eidHyperTight3MC;   //!
   TBranch        *b_LepInfo_eidHyperTight4MC;   //!
   TBranch        *b_LepInfo_ElEoverP;   //!
   TBranch        *b_LepInfo_EldeltaEta;   //!
   TBranch        *b_LepInfo_EldeltaPhi;   //!
   TBranch        *b_LepInfo_ElHadoverEm;   //!
   TBranch        *b_LepInfo_ElsigmaIetaIeta;   //!
   TBranch        *b_LepInfo_ElscSigmaIetaIeta;   //!
   TBranch        *b_LepInfo_ElEnergyErr;   //!
   TBranch        *b_LepInfo_ElMomentumErr;   //!
   TBranch        *b_LepInfo_ElTrackNHits;   //!
   TBranch        *b_LepInfo_ElSharedHitsFraction;   //!
   TBranch        *b_LepInfo_dR_gsf_ctfTrack;   //!
   TBranch        *b_LepInfo_dPt_gsf_ctfTrack;   //!
   TBranch        *b_LepInfo_ElNClusters;   //!
   TBranch        *b_LepInfo_ElClassification;   //!
   TBranch        *b_LepInfo_ElFBrem;   //!
   TBranch        *b_LepInfo_ElNumberOfBrems;   //!
   TBranch        *b_LepInfo_NumberOfExpectedInnerHits;   //!
   TBranch        *b_LepInfo_Eldist;   //!
   TBranch        *b_LepInfo_Eldcot;   //!
   TBranch        *b_LepInfo_Elconvradius;   //!
   TBranch        *b_LepInfo_ElConvPoint_x;   //!
   TBranch        *b_LepInfo_ElConvPoint_y;   //!
   TBranch        *b_LepInfo_ElConvPoint_z;   //!
   TBranch        *b_LepInfo_dcotdist;   //!
   TBranch        *b_LepInfo_ElseedEoverP;   //!
   TBranch        *b_LepInfo_ElEcalIso04;   //!
   TBranch        *b_LepInfo_ElHcalIso04;   //!
   TBranch        *b_LepInfo_GenPt;   //!
   TBranch        *b_LepInfo_GenEta;   //!
   TBranch        *b_LepInfo_GenPhi;   //!
   TBranch        *b_LepInfo_GenPdgID;   //!
   TBranch        *b_LepInfo_GenMCTag;   //!
   TBranch        *b_LepInfo_TrgPt;   //!
   TBranch        *b_LepInfo_TrgEta;   //!
   TBranch        *b_LepInfo_TrgPhi;   //!
   TBranch        *b_LepInfo_TrgID;   //!
   TBranch        *b_VertexInfo_Size;   //!
   TBranch        *b_VertexInfo_isValid;   //!
   TBranch        *b_VertexInfo_isFake;   //!
   TBranch        *b_VertexInfo_Type;   //!
   TBranch        *b_VertexInfo_Ndof;   //!
   TBranch        *b_VertexInfo_NormalizedChi2;   //!
   TBranch        *b_VertexInfo_Pt_Sum;   //!
   TBranch        *b_VertexInfo_x;   //!
   TBranch        *b_VertexInfo_y;   //!
   TBranch        *b_VertexInfo_z;   //!
   TBranch        *b_VertexInfo_Rho;   //!
   TBranch        *b_PFJetInfo_Size;   //!
   TBranch        *b_PFJetInfo_Index;   //!
   TBranch        *b_PFJetInfo_NTracks;   //!
   TBranch        *b_PFJetInfo_Et;   //!
   TBranch        *b_PFJetInfo_Pt;   //!
   TBranch        *b_PFJetInfo_Unc;   //!
   TBranch        *b_PFJetInfo_Eta;   //!
   TBranch        *b_PFJetInfo_Phi;   //!
   TBranch        *b_PFJetInfo_JetIDLOOSE;   //!
   TBranch        *b_PFJetInfo_JetCharge;   //!
   TBranch        *b_PFJetInfo_NConstituents;   //!
   TBranch        *b_PFJetInfo_NCH;   //!
   TBranch        *b_PFJetInfo_CEF;   //!
   TBranch        *b_PFJetInfo_NHF;   //!
   TBranch        *b_PFJetInfo_NEF;   //!
   TBranch        *b_PFJetInfo_CHF;   //!
   TBranch        *b_PFJetInfo_JVAlpha;   //!
   TBranch        *b_PFJetInfo_JVBeta;   //!
   TBranch        *b_PFJetInfo_PtCorrRaw;   //!
   TBranch        *b_PFJetInfo_PtCorrL2;   //!
   TBranch        *b_PFJetInfo_PtCorrL3;   //!
   TBranch        *b_PFJetInfo_PtCorrL7g;   //!
   TBranch        *b_PFJetInfo_PtCorrL7uds;   //!
   TBranch        *b_PFJetInfo_PtCorrL7c;   //!
   TBranch        *b_PFJetInfo_PtCorrL7b;   //!
   TBranch        *b_PFJetInfo_JetBProbBJetTags;   //!
   TBranch        *b_PFJetInfo_JetProbBJetTags;   //!
   TBranch        *b_PFJetInfo_TrackCountHiPurBJetTags;   //!
   TBranch        *b_PFJetInfo_TrackCountHiEffBJetTags;   //!
   TBranch        *b_PFJetInfo_SimpleSVBJetTags;   //!
   TBranch        *b_PFJetInfo_SimpleSVHEBJetTags;   //!
   TBranch        *b_PFJetInfo_SimpleSVHPBJetTags;   //!
   TBranch        *b_PFJetInfo_CombinedSVBJetTags;   //!
   TBranch        *b_PFJetInfo_CombinedSVMVABJetTags;   //!
   TBranch        *b_PFJetInfo_SoftElecByIP3dBJetTags;   //!
   TBranch        *b_PFJetInfo_SoftElecByPtBJetTags;   //!
   TBranch        *b_PFJetInfo_SoftMuonBJetTags;   //!
   TBranch        *b_PFJetInfo_SoftMuonByIP3dBJetTags;   //!
   TBranch        *b_PFJetInfo_SoftMuonByPtBJetTags;   //!
   TBranch        *b_PFJetInfo_GenJetPt;   //!
   TBranch        *b_PFJetInfo_GenJetEta;   //!
   TBranch        *b_PFJetInfo_GenJetPhi;   //!
   TBranch        *b_PFJetInfo_GenPt;   //!
   TBranch        *b_PFJetInfo_GenEta;   //!
   TBranch        *b_PFJetInfo_GenPhi;   //!
   TBranch        *b_PFJetInfo_GenPdgID;   //!
   TBranch        *b_PFJetInfo_GenFlavor;   //!
   TBranch        *b_PFJetInfo_GenMCTag;   //!
   TBranch        *b_PFJetInfo_Px;   //!
   TBranch        *b_PFJetInfo_Py;   //!
   TBranch        *b_PFJetInfo_Pz;   //!
   TBranch        *b_PFJetInfo_Energy;   //!
   TBranch        *b_JetInfo_Size;   //!
   TBranch        *b_JetInfo_Index;   //!
   TBranch        *b_JetInfo_NTracks;   //!
   TBranch        *b_JetInfo_Et;   //!
   TBranch        *b_JetInfo_Pt;   //!
   TBranch        *b_JetInfo_Unc;   //!
   TBranch        *b_JetInfo_Eta;   //!
   TBranch        *b_JetInfo_Phi;   //!
   TBranch        *b_JetInfo_JetIDLOOSE;   //!
   TBranch        *b_JetInfo_JetCharge;   //!
   TBranch        *b_JetInfo_NConstituents;   //!
   TBranch        *b_JetInfo_NCH;   //!
   TBranch        *b_JetInfo_CEF;   //!
   TBranch        *b_JetInfo_NHF;   //!
   TBranch        *b_JetInfo_NEF;   //!
   TBranch        *b_JetInfo_CHF;   //!
   TBranch        *b_JetInfo_JVAlpha;   //!
   TBranch        *b_JetInfo_JVBeta;   //!
   TBranch        *b_JetInfo_PtCorrRaw;   //!
   TBranch        *b_JetInfo_PtCorrL2;   //!
   TBranch        *b_JetInfo_PtCorrL3;   //!
   TBranch        *b_JetInfo_PtCorrL7g;   //!
   TBranch        *b_JetInfo_PtCorrL7uds;   //!
   TBranch        *b_JetInfo_PtCorrL7c;   //!
   TBranch        *b_JetInfo_PtCorrL7b;   //!
   TBranch        *b_JetInfo_JetBProbBJetTags;   //!
   TBranch        *b_JetInfo_JetProbBJetTags;   //!
   TBranch        *b_JetInfo_TrackCountHiPurBJetTags;   //!
   TBranch        *b_JetInfo_TrackCountHiEffBJetTags;   //!
   TBranch        *b_JetInfo_SimpleSVBJetTags;   //!
   TBranch        *b_JetInfo_SimpleSVHEBJetTags;   //!
   TBranch        *b_JetInfo_SimpleSVHPBJetTags;   //!
   TBranch        *b_JetInfo_CombinedSVBJetTags;   //!
   TBranch        *b_JetInfo_CombinedSVMVABJetTags;   //!
   TBranch        *b_JetInfo_SoftElecByIP3dBJetTags;   //!
   TBranch        *b_JetInfo_SoftElecByPtBJetTags;   //!
   TBranch        *b_JetInfo_SoftMuonBJetTags;   //!
   TBranch        *b_JetInfo_SoftMuonByIP3dBJetTags;   //!
   TBranch        *b_JetInfo_SoftMuonByPtBJetTags;   //!
   TBranch        *b_JetInfo_GenJetPt;   //!
   TBranch        *b_JetInfo_GenJetEta;   //!
   TBranch        *b_JetInfo_GenJetPhi;   //!
   TBranch        *b_JetInfo_GenPt;   //!
   TBranch        *b_JetInfo_GenEta;   //!
   TBranch        *b_JetInfo_GenPhi;   //!
   TBranch        *b_JetInfo_GenPdgID;   //!
   TBranch        *b_JetInfo_GenFlavor;   //!
   TBranch        *b_JetInfo_GenMCTag;   //!
   TBranch        *b_JetInfo_Px;   //!
   TBranch        *b_JetInfo_Py;   //!
   TBranch        *b_JetInfo_Pz;   //!
   TBranch        *b_JetInfo_Energy;   //!
   TBranch        *b_PairInfo_Size;   //!
   TBranch        *b_PairInfo_Index;   //!
   TBranch        *b_PairInfo_Type;   //!
   TBranch        *b_PairInfo_Obj1Index;   //!
   TBranch        *b_PairInfo_Obj2Index;   //!
   TBranch        *b_PairInfo_Mass;   //!
   TBranch        *b_PairInfo_Pt;   //!
   TBranch        *b_PairInfo_Eta;   //!
   TBranch        *b_PairInfo_Phi;   //!
   TBranch        *b_PairInfo_GenMass;   //!
   TBranch        *b_PairInfo_GenPt;   //!
   TBranch        *b_PairInfo_GenEta;   //!
   TBranch        *b_PairInfo_GenPhi;   //!
   TBranch        *b_PairInfo_GenPdgID;   //!

   bPrimeNtupleMC(TTree *tree=0);
   virtual ~bPrimeNtupleMC();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef bPrimeNtupleMC_cxx
bPrimeNtupleMC::bPrimeNtupleMC(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("results_9_1_sMl.root");
      if (!f) {
         f = new TFile("results_9_1_sMl.root");
         f->cd("results_9_1_sMl.root:/bprimeKit");
      }
      tree = (TTree*)gDirectory->Get("root");

   }
   Init(tree);
}

bPrimeNtupleMC::~bPrimeNtupleMC()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t bPrimeNtupleMC::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t bPrimeNtupleMC::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void bPrimeNtupleMC::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EvtInfo.RunNo", &EvtInfo_RunNo, &b_EvtInfo_RunNo);
   fChain->SetBranchAddress("EvtInfo.EvtNo", &EvtInfo_EvtNo, &b_EvtInfo_EvtNo);
   fChain->SetBranchAddress("EvtInfo.BxNo", &EvtInfo_BxNo, &b_EvtInfo_BxNo);
   fChain->SetBranchAddress("EvtInfo.LumiNo", &EvtInfo_LumiNo, &b_EvtInfo_LumiNo);
   fChain->SetBranchAddress("EvtInfo.Orbit", &EvtInfo_Orbit, &b_EvtInfo_Orbit);
   fChain->SetBranchAddress("EvtInfo.McFlag", &EvtInfo_McFlag, &b_EvtInfo_McFlag);
   fChain->SetBranchAddress("EvtInfo.McSigTag", &EvtInfo_McSigTag, &b_EvtInfo_McSigTag);
   fChain->SetBranchAddress("EvtInfo.McbprimeMode", EvtInfo_McbprimeMode, &b_EvtInfo_McbprimeMode);
   fChain->SetBranchAddress("EvtInfo.MctprimeMode", EvtInfo_MctprimeMode, &b_EvtInfo_MctprimeMode);
   fChain->SetBranchAddress("EvtInfo.McWMode", EvtInfo_McWMode, &b_EvtInfo_McWMode);
   fChain->SetBranchAddress("EvtInfo.McZMode", EvtInfo_McZMode, &b_EvtInfo_McZMode);
   fChain->SetBranchAddress("EvtInfo.McbprimeMass", EvtInfo_McbprimeMass, &b_EvtInfo_McbprimeMass);
   fChain->SetBranchAddress("EvtInfo.MctprimeMass", EvtInfo_MctprimeMass, &b_EvtInfo_MctprimeMass);
   fChain->SetBranchAddress("EvtInfo.MctopMass", EvtInfo_MctopMass, &b_EvtInfo_MctopMass);
   fChain->SetBranchAddress("EvtInfo.McWMass", EvtInfo_McWMass, &b_EvtInfo_McWMass);
   fChain->SetBranchAddress("EvtInfo.McZMass", EvtInfo_McZMass, &b_EvtInfo_McZMass);
   fChain->SetBranchAddress("EvtInfo.McDauPt", EvtInfo_McDauPt, &b_EvtInfo_McDauPt);
   fChain->SetBranchAddress("EvtInfo.McDauEta", EvtInfo_McDauEta, &b_EvtInfo_McDauEta);
   fChain->SetBranchAddress("EvtInfo.McDauPhi", EvtInfo_McDauPhi, &b_EvtInfo_McDauPhi);
   fChain->SetBranchAddress("EvtInfo.McDauPdgID", EvtInfo_McDauPdgID, &b_EvtInfo_McDauPdgID);
   fChain->SetBranchAddress("EvtInfo.PDFid1", &EvtInfo_PDFid1, &b_EvtInfo_PDFid1);
   fChain->SetBranchAddress("EvtInfo.PDFid2", &EvtInfo_PDFid2, &b_EvtInfo_PDFid2);
   fChain->SetBranchAddress("EvtInfo.PDFx1", &EvtInfo_PDFx1, &b_EvtInfo_PDFx1);
   fChain->SetBranchAddress("EvtInfo.RhoPU", &EvtInfo_RhoPU, &b_EvtInfo_RhoPU);
   fChain->SetBranchAddress("EvtInfo.SigmaPU", &EvtInfo_SigmaPU, &b_EvtInfo_SigmaPU);
   fChain->SetBranchAddress("EvtInfo.PDFx2", &EvtInfo_PDFx2, &b_EvtInfo_PDFx2);
   fChain->SetBranchAddress("EvtInfo.PDFscale", &EvtInfo_PDFscale, &b_EvtInfo_PDFscale);
   fChain->SetBranchAddress("EvtInfo.PDFv1", &EvtInfo_PDFv1, &b_EvtInfo_PDFv1);
   fChain->SetBranchAddress("EvtInfo.PDFv2", &EvtInfo_PDFv2, &b_EvtInfo_PDFv2);
   fChain->SetBranchAddress("EvtInfo.MET", &EvtInfo_MET, &b_EvtInfo_MET);
   fChain->SetBranchAddress("EvtInfo.METPhi", &EvtInfo_METPhi, &b_EvtInfo_METPhi);
   fChain->SetBranchAddress("EvtInfo.RawMET", &EvtInfo_RawMET, &b_EvtInfo_RawMET);
   fChain->SetBranchAddress("EvtInfo.RawMETPhi", &EvtInfo_RawMETPhi, &b_EvtInfo_RawMETPhi);
   fChain->SetBranchAddress("EvtInfo.SumEt", &EvtInfo_SumEt, &b_EvtInfo_SumEt);
   fChain->SetBranchAddress("EvtInfo.METSig", &EvtInfo_METSig, &b_EvtInfo_METSig);
   fChain->SetBranchAddress("EvtInfo.eLong", &EvtInfo_eLong, &b_EvtInfo_eLong);
   fChain->SetBranchAddress("EvtInfo.MaxHadTower", &EvtInfo_MaxHadTower, &b_EvtInfo_MaxHadTower);
   fChain->SetBranchAddress("EvtInfo.MaxEmTower", &EvtInfo_MaxEmTower, &b_EvtInfo_MaxEmTower);
   fChain->SetBranchAddress("EvtInfo.FracHad", &EvtInfo_FracHad, &b_EvtInfo_FracHad);
   fChain->SetBranchAddress("EvtInfo.FracEm", &EvtInfo_FracEm, &b_EvtInfo_FracEm);
   fChain->SetBranchAddress("EvtInfo.GenMET", &EvtInfo_GenMET, &b_EvtInfo_GenMET);
   fChain->SetBranchAddress("EvtInfo.GenMETPhi", &EvtInfo_GenMETPhi, &b_EvtInfo_GenMETPhi);
   fChain->SetBranchAddress("EvtInfo.PFMET", &EvtInfo_PFMET, &b_EvtInfo_PFMET);
   fChain->SetBranchAddress("EvtInfo.PFMETPhi", &EvtInfo_PFMETPhi, &b_EvtInfo_PFMETPhi);
   fChain->SetBranchAddress("EvtInfo.PFRawMET", &EvtInfo_PFRawMET, &b_EvtInfo_PFRawMET);
   fChain->SetBranchAddress("EvtInfo.PFRawMETPhi", &EvtInfo_PFRawMETPhi, &b_EvtInfo_PFRawMETPhi);
   fChain->SetBranchAddress("EvtInfo.PFSumEt", &EvtInfo_PFSumEt, &b_EvtInfo_PFSumEt);
   fChain->SetBranchAddress("EvtInfo.PFMETSig", &EvtInfo_PFMETSig, &b_EvtInfo_PFMETSig);
   fChain->SetBranchAddress("EvtInfo.PFGenMET", &EvtInfo_PFGenMET, &b_EvtInfo_PFGenMET);
   fChain->SetBranchAddress("EvtInfo.PFGenMETPhi", &EvtInfo_PFGenMETPhi, &b_EvtInfo_PFGenMETPhi);
   fChain->SetBranchAddress("EvtInfo.PFMETx", &EvtInfo_PFMETx, &b_EvtInfo_PFMETx);
   fChain->SetBranchAddress("EvtInfo.PFMETy", &EvtInfo_PFMETy, &b_EvtInfo_PFMETy);
   fChain->SetBranchAddress("EvtInfo.TrgCount", &EvtInfo_TrgCount, &b_EvtInfo_TrgCount);
   fChain->SetBranchAddress("EvtInfo.nTrgBook", &EvtInfo_nTrgBook, &b_EvtInfo_nTrgBook);
   fChain->SetBranchAddress("EvtInfo.TrgBook", EvtInfo_TrgBook, &b_EvtInfo_TrgBook);
   fChain->SetBranchAddress("EvtInfo.L1", EvtInfo_L1, &b_EvtInfo_L1);
   fChain->SetBranchAddress("EvtInfo.TT", EvtInfo_TT, &b_EvtInfo_TT);
   fChain->SetBranchAddress("EvtInfo.HighPurityFraction", &EvtInfo_HighPurityFraction, &b_EvtInfo_HighPurityFraction);
   fChain->SetBranchAddress("EvtInfo.NofTracks", &EvtInfo_NofTracks, &b_EvtInfo_NofTracks);
   fChain->SetBranchAddress("EvtInfo.nHLT", &EvtInfo_nHLT, &b_EvtInfo_nHLT);
   fChain->SetBranchAddress("EvtInfo.HLTbits", EvtInfo_HLTbits, &b_EvtInfo_HLTbits);
   fChain->SetBranchAddress("EvtInfo.nBX", &EvtInfo_nBX, &b_EvtInfo_nBX);
   fChain->SetBranchAddress("EvtInfo.nPU", EvtInfo_nPU, &b_EvtInfo_nPU);
   fChain->SetBranchAddress("EvtInfo.BXPU", EvtInfo_BXPU, &b_EvtInfo_BXPU);
   fChain->SetBranchAddress("GenInfo.Size", &GenInfo_Size, &b_GenInfo_Size);
   fChain->SetBranchAddress("GenInfo.Pt", GenInfo_Pt, &b_GenInfo_Pt);
   fChain->SetBranchAddress("GenInfo.Eta", GenInfo_Eta, &b_GenInfo_Eta);
   fChain->SetBranchAddress("GenInfo.Phi", GenInfo_Phi, &b_GenInfo_Phi);
   fChain->SetBranchAddress("GenInfo.Mass", GenInfo_Mass, &b_GenInfo_Mass);
   fChain->SetBranchAddress("GenInfo.PdgID", GenInfo_PdgID, &b_GenInfo_PdgID);
   fChain->SetBranchAddress("GenInfo.Status", GenInfo_Status, &b_GenInfo_Status);
   fChain->SetBranchAddress("GenInfo.nMo", GenInfo_nMo, &b_GenInfo_nMo);
   fChain->SetBranchAddress("GenInfo.nDa", GenInfo_nDa, &b_GenInfo_nDa);
   fChain->SetBranchAddress("GenInfo.Mo1", GenInfo_Mo1, &b_GenInfo_Mo1);
   fChain->SetBranchAddress("GenInfo.Mo2", GenInfo_Mo2, &b_GenInfo_Mo2);
   fChain->SetBranchAddress("GenInfo.Da1", GenInfo_Da1, &b_GenInfo_Da1);
   fChain->SetBranchAddress("GenInfo.Da2", GenInfo_Da2, &b_GenInfo_Da2);
   fChain->SetBranchAddress("PFLepInfo.Size", &PFLepInfo_Size, &b_PFLepInfo_Size);
   fChain->SetBranchAddress("PFLepInfo.Index", PFLepInfo_Index, &b_PFLepInfo_Index);
   fChain->SetBranchAddress("PFLepInfo.isEcalDriven", PFLepInfo_isEcalDriven, &b_PFLepInfo_isEcalDriven);
   fChain->SetBranchAddress("PFLepInfo.isTrackerDriven", PFLepInfo_isTrackerDriven, &b_PFLepInfo_isTrackerDriven);
   fChain->SetBranchAddress("PFLepInfo.LeptonType", PFLepInfo_LeptonType, &b_PFLepInfo_LeptonType);
   fChain->SetBranchAddress("PFLepInfo.Charge", PFLepInfo_Charge, &b_PFLepInfo_Charge);
   fChain->SetBranchAddress("PFLepInfo.ChargeGsf", PFLepInfo_ChargeGsf, &b_PFLepInfo_ChargeGsf);
   fChain->SetBranchAddress("PFLepInfo.ChargeCtf", PFLepInfo_ChargeCtf, &b_PFLepInfo_ChargeCtf);
   fChain->SetBranchAddress("PFLepInfo.ChargeScPix", PFLepInfo_ChargeScPix, &b_PFLepInfo_ChargeScPix);
   fChain->SetBranchAddress("PFLepInfo.Pt", PFLepInfo_Pt, &b_PFLepInfo_Pt);
   fChain->SetBranchAddress("PFLepInfo.Et", PFLepInfo_Et, &b_PFLepInfo_Et);
   fChain->SetBranchAddress("PFLepInfo.Eta", PFLepInfo_Eta, &b_PFLepInfo_Eta);
   fChain->SetBranchAddress("PFLepInfo.caloEta", PFLepInfo_caloEta, &b_PFLepInfo_caloEta);
   fChain->SetBranchAddress("PFLepInfo.Phi", PFLepInfo_Phi, &b_PFLepInfo_Phi);
   fChain->SetBranchAddress("PFLepInfo.TrackIso", PFLepInfo_TrackIso, &b_PFLepInfo_TrackIso);
   fChain->SetBranchAddress("PFLepInfo.EcalIso", PFLepInfo_EcalIso, &b_PFLepInfo_EcalIso);
   fChain->SetBranchAddress("PFLepInfo.HcalIso", PFLepInfo_HcalIso, &b_PFLepInfo_HcalIso);
   fChain->SetBranchAddress("PFLepInfo.HcalDepth1Iso", PFLepInfo_HcalDepth1Iso, &b_PFLepInfo_HcalDepth1Iso);
   fChain->SetBranchAddress("PFLepInfo.HcalDepth2Iso", PFLepInfo_HcalDepth2Iso, &b_PFLepInfo_HcalDepth2Iso);
   fChain->SetBranchAddress("PFLepInfo.ChargedHadronIso", PFLepInfo_ChargedHadronIso, &b_PFLepInfo_ChargedHadronIso);
   fChain->SetBranchAddress("PFLepInfo.NeutralHadronIso", PFLepInfo_NeutralHadronIso, &b_PFLepInfo_NeutralHadronIso);
   fChain->SetBranchAddress("PFLepInfo.PhotonIso", PFLepInfo_PhotonIso, &b_PFLepInfo_PhotonIso);
   fChain->SetBranchAddress("PFLepInfo.CaloEnergy", PFLepInfo_CaloEnergy, &b_PFLepInfo_CaloEnergy);
   fChain->SetBranchAddress("PFLepInfo.e1x5", PFLepInfo_e1x5, &b_PFLepInfo_e1x5);
   fChain->SetBranchAddress("PFLepInfo.e2x5Max", PFLepInfo_e2x5Max, &b_PFLepInfo_e2x5Max);
   fChain->SetBranchAddress("PFLepInfo.e5x5", PFLepInfo_e5x5, &b_PFLepInfo_e5x5);
   fChain->SetBranchAddress("PFLepInfo.Px", PFLepInfo_Px, &b_PFLepInfo_Px);
   fChain->SetBranchAddress("PFLepInfo.Py", PFLepInfo_Py, &b_PFLepInfo_Py);
   fChain->SetBranchAddress("PFLepInfo.Pz", PFLepInfo_Pz, &b_PFLepInfo_Pz);
   fChain->SetBranchAddress("PFLepInfo.Energy", PFLepInfo_Energy, &b_PFLepInfo_Energy);
   fChain->SetBranchAddress("PFLepInfo.vertexZ", PFLepInfo_vertexZ, &b_PFLepInfo_vertexZ);
   fChain->SetBranchAddress("PFLepInfo.MuIDGlobalMuonPromptTight", PFLepInfo_MuIDGlobalMuonPromptTight, &b_PFLepInfo_MuIDGlobalMuonPromptTight);
   fChain->SetBranchAddress("PFLepInfo.MuInnerPtError", PFLepInfo_MuInnerPtError, &b_PFLepInfo_MuInnerPtError);
   fChain->SetBranchAddress("PFLepInfo.MuGlobalPtError", PFLepInfo_MuGlobalPtError, &b_PFLepInfo_MuGlobalPtError);
   fChain->SetBranchAddress("PFLepInfo.MuInnerTrackDz", PFLepInfo_MuInnerTrackDz, &b_PFLepInfo_MuInnerTrackDz);
   fChain->SetBranchAddress("PFLepInfo.MuInnerTrackD0", PFLepInfo_MuInnerTrackD0, &b_PFLepInfo_MuInnerTrackD0);
   fChain->SetBranchAddress("PFLepInfo.MuInnerTrackDxy_BS", PFLepInfo_MuInnerTrackDxy_BS, &b_PFLepInfo_MuInnerTrackDxy_BS);
   fChain->SetBranchAddress("PFLepInfo.MuInnerTrackDxy_PV", PFLepInfo_MuInnerTrackDxy_PV, &b_PFLepInfo_MuInnerTrackDxy_PV);
   fChain->SetBranchAddress("PFLepInfo.MuInnerTrackDxy_PVBS", PFLepInfo_MuInnerTrackDxy_PVBS, &b_PFLepInfo_MuInnerTrackDxy_PVBS);
   fChain->SetBranchAddress("PFLepInfo.MuInnerTrackNHits", PFLepInfo_MuInnerTrackNHits, &b_PFLepInfo_MuInnerTrackNHits);
   fChain->SetBranchAddress("PFLepInfo.MuNTrackerHits", PFLepInfo_MuNTrackerHits, &b_PFLepInfo_MuNTrackerHits);
   fChain->SetBranchAddress("PFLepInfo.MuCaloCompat", PFLepInfo_MuCaloCompat, &b_PFLepInfo_MuCaloCompat);
   fChain->SetBranchAddress("PFLepInfo.MuNChambers", PFLepInfo_MuNChambers, &b_PFLepInfo_MuNChambers);
   fChain->SetBranchAddress("PFLepInfo.MuNChambersMatchesSegment", PFLepInfo_MuNChambersMatchesSegment, &b_PFLepInfo_MuNChambersMatchesSegment);
   fChain->SetBranchAddress("PFLepInfo.MuNPixelLayers", PFLepInfo_MuNPixelLayers, &b_PFLepInfo_MuNPixelLayers);
   fChain->SetBranchAddress("PFLepInfo.MuNPixelLayersWMeasurement", PFLepInfo_MuNPixelLayersWMeasurement, &b_PFLepInfo_MuNPixelLayersWMeasurement);
   fChain->SetBranchAddress("PFLepInfo.MuNLostInnerHits", PFLepInfo_MuNLostInnerHits, &b_PFLepInfo_MuNLostInnerHits);
   fChain->SetBranchAddress("PFLepInfo.MuNLostOuterHits", PFLepInfo_MuNLostOuterHits, &b_PFLepInfo_MuNLostOuterHits);
   fChain->SetBranchAddress("PFLepInfo.MuNMuonhits", PFLepInfo_MuNMuonhits, &b_PFLepInfo_MuNMuonhits);
   fChain->SetBranchAddress("PFLepInfo.MuType", PFLepInfo_MuType, &b_PFLepInfo_MuType);
   fChain->SetBranchAddress("PFLepInfo.MuGlobalNormalizedChi2", PFLepInfo_MuGlobalNormalizedChi2, &b_PFLepInfo_MuGlobalNormalizedChi2);
   fChain->SetBranchAddress("PFLepInfo.ElTrackNLostHits", PFLepInfo_ElTrackNLostHits, &b_PFLepInfo_ElTrackNLostHits);
   fChain->SetBranchAddress("PFLepInfo.ElTrackDz", PFLepInfo_ElTrackDz, &b_PFLepInfo_ElTrackDz);
   fChain->SetBranchAddress("PFLepInfo.ElTrackD0", PFLepInfo_ElTrackD0, &b_PFLepInfo_ElTrackD0);
   fChain->SetBranchAddress("PFLepInfo.ElTrackDxy_BS", PFLepInfo_ElTrackDxy_BS, &b_PFLepInfo_ElTrackDxy_BS);
   fChain->SetBranchAddress("PFLepInfo.ElTrackDxy_PV", PFLepInfo_ElTrackDxy_PV, &b_PFLepInfo_ElTrackDxy_PV);
   fChain->SetBranchAddress("PFLepInfo.ElTrackDxy_PVBS", PFLepInfo_ElTrackDxy_PVBS, &b_PFLepInfo_ElTrackDxy_PVBS);
   fChain->SetBranchAddress("PFLepInfo.simpleEleId95relIso", PFLepInfo_simpleEleId95relIso, &b_PFLepInfo_simpleEleId95relIso);
   fChain->SetBranchAddress("PFLepInfo.simpleEleId90relIso", PFLepInfo_simpleEleId90relIso, &b_PFLepInfo_simpleEleId90relIso);
   fChain->SetBranchAddress("PFLepInfo.simpleEleId85relIso", PFLepInfo_simpleEleId85relIso, &b_PFLepInfo_simpleEleId85relIso);
   fChain->SetBranchAddress("PFLepInfo.simpleEleId80relIso", PFLepInfo_simpleEleId80relIso, &b_PFLepInfo_simpleEleId80relIso);
   fChain->SetBranchAddress("PFLepInfo.simpleEleId70relIso", PFLepInfo_simpleEleId70relIso, &b_PFLepInfo_simpleEleId70relIso);
   fChain->SetBranchAddress("PFLepInfo.simpleEleId60relIso", PFLepInfo_simpleEleId60relIso, &b_PFLepInfo_simpleEleId60relIso);
   fChain->SetBranchAddress("PFLepInfo.simpleEleId95cIso", PFLepInfo_simpleEleId95cIso, &b_PFLepInfo_simpleEleId95cIso);
   fChain->SetBranchAddress("PFLepInfo.simpleEleId90cIso", PFLepInfo_simpleEleId90cIso, &b_PFLepInfo_simpleEleId90cIso);
   fChain->SetBranchAddress("PFLepInfo.simpleEleId85cIso", PFLepInfo_simpleEleId85cIso, &b_PFLepInfo_simpleEleId85cIso);
   fChain->SetBranchAddress("PFLepInfo.simpleEleId80cIso", PFLepInfo_simpleEleId80cIso, &b_PFLepInfo_simpleEleId80cIso);
   fChain->SetBranchAddress("PFLepInfo.simpleEleId70cIso", PFLepInfo_simpleEleId70cIso, &b_PFLepInfo_simpleEleId70cIso);
   fChain->SetBranchAddress("PFLepInfo.simpleEleId60cIso", PFLepInfo_simpleEleId60cIso, &b_PFLepInfo_simpleEleId60cIso);
   fChain->SetBranchAddress("PFLepInfo.eidVeryLoose", PFLepInfo_eidVeryLoose, &b_PFLepInfo_eidVeryLoose);
   fChain->SetBranchAddress("PFLepInfo.eidLoose", PFLepInfo_eidLoose, &b_PFLepInfo_eidLoose);
   fChain->SetBranchAddress("PFLepInfo.eidMedium", PFLepInfo_eidMedium, &b_PFLepInfo_eidMedium);
   fChain->SetBranchAddress("PFLepInfo.eidTight", PFLepInfo_eidTight, &b_PFLepInfo_eidTight);
   fChain->SetBranchAddress("PFLepInfo.eidSuperTight", PFLepInfo_eidSuperTight, &b_PFLepInfo_eidSuperTight);
   fChain->SetBranchAddress("PFLepInfo.eidHyperTight1", PFLepInfo_eidHyperTight1, &b_PFLepInfo_eidHyperTight1);
   fChain->SetBranchAddress("PFLepInfo.eidHyperTight2", PFLepInfo_eidHyperTight2, &b_PFLepInfo_eidHyperTight2);
   fChain->SetBranchAddress("PFLepInfo.eidHyperTight3", PFLepInfo_eidHyperTight3, &b_PFLepInfo_eidHyperTight3);
   fChain->SetBranchAddress("PFLepInfo.eidHyperTight4", PFLepInfo_eidHyperTight4, &b_PFLepInfo_eidHyperTight4);
   fChain->SetBranchAddress("PFLepInfo.eidVeryLooseMC", PFLepInfo_eidVeryLooseMC, &b_PFLepInfo_eidVeryLooseMC);
   fChain->SetBranchAddress("PFLepInfo.eidLooseMC", PFLepInfo_eidLooseMC, &b_PFLepInfo_eidLooseMC);
   fChain->SetBranchAddress("PFLepInfo.eidMediumMC", PFLepInfo_eidMediumMC, &b_PFLepInfo_eidMediumMC);
   fChain->SetBranchAddress("PFLepInfo.eidTightMC", PFLepInfo_eidTightMC, &b_PFLepInfo_eidTightMC);
   fChain->SetBranchAddress("PFLepInfo.eidSuperTightMC", PFLepInfo_eidSuperTightMC, &b_PFLepInfo_eidSuperTightMC);
   fChain->SetBranchAddress("PFLepInfo.eidHyperTight1MC", PFLepInfo_eidHyperTight1MC, &b_PFLepInfo_eidHyperTight1MC);
   fChain->SetBranchAddress("PFLepInfo.eidHyperTight2MC", PFLepInfo_eidHyperTight2MC, &b_PFLepInfo_eidHyperTight2MC);
   fChain->SetBranchAddress("PFLepInfo.eidHyperTight3MC", PFLepInfo_eidHyperTight3MC, &b_PFLepInfo_eidHyperTight3MC);
   fChain->SetBranchAddress("PFLepInfo.eidHyperTight4MC", PFLepInfo_eidHyperTight4MC, &b_PFLepInfo_eidHyperTight4MC);
   fChain->SetBranchAddress("PFLepInfo.ElEoverP", PFLepInfo_ElEoverP, &b_PFLepInfo_ElEoverP);
   fChain->SetBranchAddress("PFLepInfo.EldeltaEta", PFLepInfo_EldeltaEta, &b_PFLepInfo_EldeltaEta);
   fChain->SetBranchAddress("PFLepInfo.EldeltaPhi", PFLepInfo_EldeltaPhi, &b_PFLepInfo_EldeltaPhi);
   fChain->SetBranchAddress("PFLepInfo.ElHadoverEm", PFLepInfo_ElHadoverEm, &b_PFLepInfo_ElHadoverEm);
   fChain->SetBranchAddress("PFLepInfo.ElsigmaIetaIeta", PFLepInfo_ElsigmaIetaIeta, &b_PFLepInfo_ElsigmaIetaIeta);
   fChain->SetBranchAddress("PFLepInfo.ElscSigmaIetaIeta", PFLepInfo_ElscSigmaIetaIeta, &b_PFLepInfo_ElscSigmaIetaIeta);
   fChain->SetBranchAddress("PFLepInfo.ElEnergyErr", PFLepInfo_ElEnergyErr, &b_PFLepInfo_ElEnergyErr);
   fChain->SetBranchAddress("PFLepInfo.ElMomentumErr", PFLepInfo_ElMomentumErr, &b_PFLepInfo_ElMomentumErr);
   fChain->SetBranchAddress("PFLepInfo.ElTrackNHits", PFLepInfo_ElTrackNHits, &b_PFLepInfo_ElTrackNHits);
   fChain->SetBranchAddress("PFLepInfo.ElSharedHitsFraction", PFLepInfo_ElSharedHitsFraction, &b_PFLepInfo_ElSharedHitsFraction);
   fChain->SetBranchAddress("PFLepInfo.dR_gsf_ctfTrack", PFLepInfo_dR_gsf_ctfTrack, &b_PFLepInfo_dR_gsf_ctfTrack);
   fChain->SetBranchAddress("PFLepInfo.dPt_gsf_ctfTrack", PFLepInfo_dPt_gsf_ctfTrack, &b_PFLepInfo_dPt_gsf_ctfTrack);
   fChain->SetBranchAddress("PFLepInfo.ElNClusters", PFLepInfo_ElNClusters, &b_PFLepInfo_ElNClusters);
   fChain->SetBranchAddress("PFLepInfo.ElClassification", PFLepInfo_ElClassification, &b_PFLepInfo_ElClassification);
   fChain->SetBranchAddress("PFLepInfo.ElFBrem", PFLepInfo_ElFBrem, &b_PFLepInfo_ElFBrem);
   fChain->SetBranchAddress("PFLepInfo.ElNumberOfBrems", PFLepInfo_ElNumberOfBrems, &b_PFLepInfo_ElNumberOfBrems);
   fChain->SetBranchAddress("PFLepInfo.NumberOfExpectedInnerHits", PFLepInfo_NumberOfExpectedInnerHits, &b_PFLepInfo_NumberOfExpectedInnerHits);
   fChain->SetBranchAddress("PFLepInfo.Eldist", PFLepInfo_Eldist, &b_PFLepInfo_Eldist);
   fChain->SetBranchAddress("PFLepInfo.Eldcot", PFLepInfo_Eldcot, &b_PFLepInfo_Eldcot);
   fChain->SetBranchAddress("PFLepInfo.Elconvradius", PFLepInfo_Elconvradius, &b_PFLepInfo_Elconvradius);
   fChain->SetBranchAddress("PFLepInfo.ElConvPoint_x", PFLepInfo_ElConvPoint_x, &b_PFLepInfo_ElConvPoint_x);
   fChain->SetBranchAddress("PFLepInfo.ElConvPoint_y", PFLepInfo_ElConvPoint_y, &b_PFLepInfo_ElConvPoint_y);
   fChain->SetBranchAddress("PFLepInfo.ElConvPoint_z", PFLepInfo_ElConvPoint_z, &b_PFLepInfo_ElConvPoint_z);
   fChain->SetBranchAddress("PFLepInfo.dcotdist", PFLepInfo_dcotdist, &b_PFLepInfo_dcotdist);
   fChain->SetBranchAddress("PFLepInfo.ElseedEoverP", PFLepInfo_ElseedEoverP, &b_PFLepInfo_ElseedEoverP);
   fChain->SetBranchAddress("PFLepInfo.ElEcalIso04", PFLepInfo_ElEcalIso04, &b_PFLepInfo_ElEcalIso04);
   fChain->SetBranchAddress("PFLepInfo.ElHcalIso04", PFLepInfo_ElHcalIso04, &b_PFLepInfo_ElHcalIso04);
   fChain->SetBranchAddress("PFLepInfo.GenPt", PFLepInfo_GenPt, &b_PFLepInfo_GenPt);
   fChain->SetBranchAddress("PFLepInfo.GenEta", PFLepInfo_GenEta, &b_PFLepInfo_GenEta);
   fChain->SetBranchAddress("PFLepInfo.GenPhi", PFLepInfo_GenPhi, &b_PFLepInfo_GenPhi);
   fChain->SetBranchAddress("PFLepInfo.GenPdgID", PFLepInfo_GenPdgID, &b_PFLepInfo_GenPdgID);
   fChain->SetBranchAddress("PFLepInfo.GenMCTag", PFLepInfo_GenMCTag, &b_PFLepInfo_GenMCTag);
   fChain->SetBranchAddress("PFLepInfo.TrgPt", PFLepInfo_TrgPt, &b_PFLepInfo_TrgPt);
   fChain->SetBranchAddress("PFLepInfo.TrgEta", PFLepInfo_TrgEta, &b_PFLepInfo_TrgEta);
   fChain->SetBranchAddress("PFLepInfo.TrgPhi", PFLepInfo_TrgPhi, &b_PFLepInfo_TrgPhi);
   fChain->SetBranchAddress("PFLepInfo.TrgID", PFLepInfo_TrgID, &b_PFLepInfo_TrgID);
   fChain->SetBranchAddress("LepInfo.Size", &LepInfo_Size, &b_LepInfo_Size);
   fChain->SetBranchAddress("LepInfo.Index", LepInfo_Index, &b_LepInfo_Index);
   fChain->SetBranchAddress("LepInfo.isEcalDriven", LepInfo_isEcalDriven, &b_LepInfo_isEcalDriven);
   fChain->SetBranchAddress("LepInfo.isTrackerDriven", LepInfo_isTrackerDriven, &b_LepInfo_isTrackerDriven);
   fChain->SetBranchAddress("LepInfo.LeptonType", LepInfo_LeptonType, &b_LepInfo_LeptonType);
   fChain->SetBranchAddress("LepInfo.Charge", LepInfo_Charge, &b_LepInfo_Charge);
   fChain->SetBranchAddress("LepInfo.ChargeGsf", LepInfo_ChargeGsf, &b_LepInfo_ChargeGsf);
   fChain->SetBranchAddress("LepInfo.ChargeCtf", LepInfo_ChargeCtf, &b_LepInfo_ChargeCtf);
   fChain->SetBranchAddress("LepInfo.ChargeScPix", LepInfo_ChargeScPix, &b_LepInfo_ChargeScPix);
   fChain->SetBranchAddress("LepInfo.Pt", LepInfo_Pt, &b_LepInfo_Pt);
   fChain->SetBranchAddress("LepInfo.Et", LepInfo_Et, &b_LepInfo_Et);
   fChain->SetBranchAddress("LepInfo.Eta", LepInfo_Eta, &b_LepInfo_Eta);
   fChain->SetBranchAddress("LepInfo.caloEta", LepInfo_caloEta, &b_LepInfo_caloEta);
   fChain->SetBranchAddress("LepInfo.Phi", LepInfo_Phi, &b_LepInfo_Phi);
   fChain->SetBranchAddress("LepInfo.TrackIso", LepInfo_TrackIso, &b_LepInfo_TrackIso);
   fChain->SetBranchAddress("LepInfo.EcalIso", LepInfo_EcalIso, &b_LepInfo_EcalIso);
   fChain->SetBranchAddress("LepInfo.HcalIso", LepInfo_HcalIso, &b_LepInfo_HcalIso);
   fChain->SetBranchAddress("LepInfo.HcalDepth1Iso", LepInfo_HcalDepth1Iso, &b_LepInfo_HcalDepth1Iso);
   fChain->SetBranchAddress("LepInfo.HcalDepth2Iso", LepInfo_HcalDepth2Iso, &b_LepInfo_HcalDepth2Iso);
   fChain->SetBranchAddress("LepInfo.ChargedHadronIso", LepInfo_ChargedHadronIso, &b_LepInfo_ChargedHadronIso);
   fChain->SetBranchAddress("LepInfo.NeutralHadronIso", LepInfo_NeutralHadronIso, &b_LepInfo_NeutralHadronIso);
   fChain->SetBranchAddress("LepInfo.PhotonIso", LepInfo_PhotonIso, &b_LepInfo_PhotonIso);
   fChain->SetBranchAddress("LepInfo.CaloEnergy", LepInfo_CaloEnergy, &b_LepInfo_CaloEnergy);
   fChain->SetBranchAddress("LepInfo.e1x5", LepInfo_e1x5, &b_LepInfo_e1x5);
   fChain->SetBranchAddress("LepInfo.e2x5Max", LepInfo_e2x5Max, &b_LepInfo_e2x5Max);
   fChain->SetBranchAddress("LepInfo.e5x5", LepInfo_e5x5, &b_LepInfo_e5x5);
   fChain->SetBranchAddress("LepInfo.Px", LepInfo_Px, &b_LepInfo_Px);
   fChain->SetBranchAddress("LepInfo.Py", LepInfo_Py, &b_LepInfo_Py);
   fChain->SetBranchAddress("LepInfo.Pz", LepInfo_Pz, &b_LepInfo_Pz);
   fChain->SetBranchAddress("LepInfo.Energy", LepInfo_Energy, &b_LepInfo_Energy);
   fChain->SetBranchAddress("LepInfo.vertexZ", LepInfo_vertexZ, &b_LepInfo_vertexZ);
   fChain->SetBranchAddress("LepInfo.MuIDGlobalMuonPromptTight", LepInfo_MuIDGlobalMuonPromptTight, &b_LepInfo_MuIDGlobalMuonPromptTight);
   fChain->SetBranchAddress("LepInfo.MuInnerPtError", LepInfo_MuInnerPtError, &b_LepInfo_MuInnerPtError);
   fChain->SetBranchAddress("LepInfo.MuGlobalPtError", LepInfo_MuGlobalPtError, &b_LepInfo_MuGlobalPtError);
   fChain->SetBranchAddress("LepInfo.MuInnerTrackDz", LepInfo_MuInnerTrackDz, &b_LepInfo_MuInnerTrackDz);
   fChain->SetBranchAddress("LepInfo.MuInnerTrackD0", LepInfo_MuInnerTrackD0, &b_LepInfo_MuInnerTrackD0);
   fChain->SetBranchAddress("LepInfo.MuInnerTrackDxy_BS", LepInfo_MuInnerTrackDxy_BS, &b_LepInfo_MuInnerTrackDxy_BS);
   fChain->SetBranchAddress("LepInfo.MuInnerTrackDxy_PV", LepInfo_MuInnerTrackDxy_PV, &b_LepInfo_MuInnerTrackDxy_PV);
   fChain->SetBranchAddress("LepInfo.MuInnerTrackDxy_PVBS", LepInfo_MuInnerTrackDxy_PVBS, &b_LepInfo_MuInnerTrackDxy_PVBS);
   fChain->SetBranchAddress("LepInfo.MuInnerTrackNHits", LepInfo_MuInnerTrackNHits, &b_LepInfo_MuInnerTrackNHits);
   fChain->SetBranchAddress("LepInfo.MuNTrackerHits", LepInfo_MuNTrackerHits, &b_LepInfo_MuNTrackerHits);
   fChain->SetBranchAddress("LepInfo.MuCaloCompat", LepInfo_MuCaloCompat, &b_LepInfo_MuCaloCompat);
   fChain->SetBranchAddress("LepInfo.MuNChambers", LepInfo_MuNChambers, &b_LepInfo_MuNChambers);
   fChain->SetBranchAddress("LepInfo.MuNChambersMatchesSegment", LepInfo_MuNChambersMatchesSegment, &b_LepInfo_MuNChambersMatchesSegment);
   fChain->SetBranchAddress("LepInfo.MuNPixelLayers", LepInfo_MuNPixelLayers, &b_LepInfo_MuNPixelLayers);
   fChain->SetBranchAddress("LepInfo.MuNPixelLayersWMeasurement", LepInfo_MuNPixelLayersWMeasurement, &b_LepInfo_MuNPixelLayersWMeasurement);
   fChain->SetBranchAddress("LepInfo.MuNLostInnerHits", LepInfo_MuNLostInnerHits, &b_LepInfo_MuNLostInnerHits);
   fChain->SetBranchAddress("LepInfo.MuNLostOuterHits", LepInfo_MuNLostOuterHits, &b_LepInfo_MuNLostOuterHits);
   fChain->SetBranchAddress("LepInfo.MuNMuonhits", LepInfo_MuNMuonhits, &b_LepInfo_MuNMuonhits);
   fChain->SetBranchAddress("LepInfo.MuType", LepInfo_MuType, &b_LepInfo_MuType);
   fChain->SetBranchAddress("LepInfo.MuGlobalNormalizedChi2", LepInfo_MuGlobalNormalizedChi2, &b_LepInfo_MuGlobalNormalizedChi2);
   fChain->SetBranchAddress("LepInfo.ElTrackNLostHits", LepInfo_ElTrackNLostHits, &b_LepInfo_ElTrackNLostHits);
   fChain->SetBranchAddress("LepInfo.ElTrackDz", LepInfo_ElTrackDz, &b_LepInfo_ElTrackDz);
   fChain->SetBranchAddress("LepInfo.ElTrackD0", LepInfo_ElTrackD0, &b_LepInfo_ElTrackD0);
   fChain->SetBranchAddress("LepInfo.ElTrackDxy_BS", LepInfo_ElTrackDxy_BS, &b_LepInfo_ElTrackDxy_BS);
   fChain->SetBranchAddress("LepInfo.ElTrackDxy_PV", LepInfo_ElTrackDxy_PV, &b_LepInfo_ElTrackDxy_PV);
   fChain->SetBranchAddress("LepInfo.ElTrackDxy_PVBS", LepInfo_ElTrackDxy_PVBS, &b_LepInfo_ElTrackDxy_PVBS);
   fChain->SetBranchAddress("LepInfo.simpleEleId95relIso", LepInfo_simpleEleId95relIso, &b_LepInfo_simpleEleId95relIso);
   fChain->SetBranchAddress("LepInfo.simpleEleId90relIso", LepInfo_simpleEleId90relIso, &b_LepInfo_simpleEleId90relIso);
   fChain->SetBranchAddress("LepInfo.simpleEleId85relIso", LepInfo_simpleEleId85relIso, &b_LepInfo_simpleEleId85relIso);
   fChain->SetBranchAddress("LepInfo.simpleEleId80relIso", LepInfo_simpleEleId80relIso, &b_LepInfo_simpleEleId80relIso);
   fChain->SetBranchAddress("LepInfo.simpleEleId70relIso", LepInfo_simpleEleId70relIso, &b_LepInfo_simpleEleId70relIso);
   fChain->SetBranchAddress("LepInfo.simpleEleId60relIso", LepInfo_simpleEleId60relIso, &b_LepInfo_simpleEleId60relIso);
   fChain->SetBranchAddress("LepInfo.simpleEleId95cIso", LepInfo_simpleEleId95cIso, &b_LepInfo_simpleEleId95cIso);
   fChain->SetBranchAddress("LepInfo.simpleEleId90cIso", LepInfo_simpleEleId90cIso, &b_LepInfo_simpleEleId90cIso);
   fChain->SetBranchAddress("LepInfo.simpleEleId85cIso", LepInfo_simpleEleId85cIso, &b_LepInfo_simpleEleId85cIso);
   fChain->SetBranchAddress("LepInfo.simpleEleId80cIso", LepInfo_simpleEleId80cIso, &b_LepInfo_simpleEleId80cIso);
   fChain->SetBranchAddress("LepInfo.simpleEleId70cIso", LepInfo_simpleEleId70cIso, &b_LepInfo_simpleEleId70cIso);
   fChain->SetBranchAddress("LepInfo.simpleEleId60cIso", LepInfo_simpleEleId60cIso, &b_LepInfo_simpleEleId60cIso);
   fChain->SetBranchAddress("LepInfo.eidVeryLoose", LepInfo_eidVeryLoose, &b_LepInfo_eidVeryLoose);
   fChain->SetBranchAddress("LepInfo.eidLoose", LepInfo_eidLoose, &b_LepInfo_eidLoose);
   fChain->SetBranchAddress("LepInfo.eidMedium", LepInfo_eidMedium, &b_LepInfo_eidMedium);
   fChain->SetBranchAddress("LepInfo.eidTight", LepInfo_eidTight, &b_LepInfo_eidTight);
   fChain->SetBranchAddress("LepInfo.eidSuperTight", LepInfo_eidSuperTight, &b_LepInfo_eidSuperTight);
   fChain->SetBranchAddress("LepInfo.eidHyperTight1", LepInfo_eidHyperTight1, &b_LepInfo_eidHyperTight1);
   fChain->SetBranchAddress("LepInfo.eidHyperTight2", LepInfo_eidHyperTight2, &b_LepInfo_eidHyperTight2);
   fChain->SetBranchAddress("LepInfo.eidHyperTight3", LepInfo_eidHyperTight3, &b_LepInfo_eidHyperTight3);
   fChain->SetBranchAddress("LepInfo.eidHyperTight4", LepInfo_eidHyperTight4, &b_LepInfo_eidHyperTight4);
   fChain->SetBranchAddress("LepInfo.eidVeryLooseMC", LepInfo_eidVeryLooseMC, &b_LepInfo_eidVeryLooseMC);
   fChain->SetBranchAddress("LepInfo.eidLooseMC", LepInfo_eidLooseMC, &b_LepInfo_eidLooseMC);
   fChain->SetBranchAddress("LepInfo.eidMediumMC", LepInfo_eidMediumMC, &b_LepInfo_eidMediumMC);
   fChain->SetBranchAddress("LepInfo.eidTightMC", LepInfo_eidTightMC, &b_LepInfo_eidTightMC);
   fChain->SetBranchAddress("LepInfo.eidSuperTightMC", LepInfo_eidSuperTightMC, &b_LepInfo_eidSuperTightMC);
   fChain->SetBranchAddress("LepInfo.eidHyperTight1MC", LepInfo_eidHyperTight1MC, &b_LepInfo_eidHyperTight1MC);
   fChain->SetBranchAddress("LepInfo.eidHyperTight2MC", LepInfo_eidHyperTight2MC, &b_LepInfo_eidHyperTight2MC);
   fChain->SetBranchAddress("LepInfo.eidHyperTight3MC", LepInfo_eidHyperTight3MC, &b_LepInfo_eidHyperTight3MC);
   fChain->SetBranchAddress("LepInfo.eidHyperTight4MC", LepInfo_eidHyperTight4MC, &b_LepInfo_eidHyperTight4MC);
   fChain->SetBranchAddress("LepInfo.ElEoverP", LepInfo_ElEoverP, &b_LepInfo_ElEoverP);
   fChain->SetBranchAddress("LepInfo.EldeltaEta", LepInfo_EldeltaEta, &b_LepInfo_EldeltaEta);
   fChain->SetBranchAddress("LepInfo.EldeltaPhi", LepInfo_EldeltaPhi, &b_LepInfo_EldeltaPhi);
   fChain->SetBranchAddress("LepInfo.ElHadoverEm", LepInfo_ElHadoverEm, &b_LepInfo_ElHadoverEm);
   fChain->SetBranchAddress("LepInfo.ElsigmaIetaIeta", LepInfo_ElsigmaIetaIeta, &b_LepInfo_ElsigmaIetaIeta);
   fChain->SetBranchAddress("LepInfo.ElscSigmaIetaIeta", LepInfo_ElscSigmaIetaIeta, &b_LepInfo_ElscSigmaIetaIeta);
   fChain->SetBranchAddress("LepInfo.ElEnergyErr", LepInfo_ElEnergyErr, &b_LepInfo_ElEnergyErr);
   fChain->SetBranchAddress("LepInfo.ElMomentumErr", LepInfo_ElMomentumErr, &b_LepInfo_ElMomentumErr);
   fChain->SetBranchAddress("LepInfo.ElTrackNHits", LepInfo_ElTrackNHits, &b_LepInfo_ElTrackNHits);
   fChain->SetBranchAddress("LepInfo.ElSharedHitsFraction", LepInfo_ElSharedHitsFraction, &b_LepInfo_ElSharedHitsFraction);
   fChain->SetBranchAddress("LepInfo.dR_gsf_ctfTrack", LepInfo_dR_gsf_ctfTrack, &b_LepInfo_dR_gsf_ctfTrack);
   fChain->SetBranchAddress("LepInfo.dPt_gsf_ctfTrack", LepInfo_dPt_gsf_ctfTrack, &b_LepInfo_dPt_gsf_ctfTrack);
   fChain->SetBranchAddress("LepInfo.ElNClusters", LepInfo_ElNClusters, &b_LepInfo_ElNClusters);
   fChain->SetBranchAddress("LepInfo.ElClassification", LepInfo_ElClassification, &b_LepInfo_ElClassification);
   fChain->SetBranchAddress("LepInfo.ElFBrem", LepInfo_ElFBrem, &b_LepInfo_ElFBrem);
   fChain->SetBranchAddress("LepInfo.ElNumberOfBrems", LepInfo_ElNumberOfBrems, &b_LepInfo_ElNumberOfBrems);
   fChain->SetBranchAddress("LepInfo.NumberOfExpectedInnerHits", LepInfo_NumberOfExpectedInnerHits, &b_LepInfo_NumberOfExpectedInnerHits);
   fChain->SetBranchAddress("LepInfo.Eldist", LepInfo_Eldist, &b_LepInfo_Eldist);
   fChain->SetBranchAddress("LepInfo.Eldcot", LepInfo_Eldcot, &b_LepInfo_Eldcot);
   fChain->SetBranchAddress("LepInfo.Elconvradius", LepInfo_Elconvradius, &b_LepInfo_Elconvradius);
   fChain->SetBranchAddress("LepInfo.ElConvPoint_x", LepInfo_ElConvPoint_x, &b_LepInfo_ElConvPoint_x);
   fChain->SetBranchAddress("LepInfo.ElConvPoint_y", LepInfo_ElConvPoint_y, &b_LepInfo_ElConvPoint_y);
   fChain->SetBranchAddress("LepInfo.ElConvPoint_z", LepInfo_ElConvPoint_z, &b_LepInfo_ElConvPoint_z);
   fChain->SetBranchAddress("LepInfo.dcotdist", LepInfo_dcotdist, &b_LepInfo_dcotdist);
   fChain->SetBranchAddress("LepInfo.ElseedEoverP", LepInfo_ElseedEoverP, &b_LepInfo_ElseedEoverP);
   fChain->SetBranchAddress("LepInfo.ElEcalIso04", LepInfo_ElEcalIso04, &b_LepInfo_ElEcalIso04);
   fChain->SetBranchAddress("LepInfo.ElHcalIso04", LepInfo_ElHcalIso04, &b_LepInfo_ElHcalIso04);
   fChain->SetBranchAddress("LepInfo.GenPt", LepInfo_GenPt, &b_LepInfo_GenPt);
   fChain->SetBranchAddress("LepInfo.GenEta", LepInfo_GenEta, &b_LepInfo_GenEta);
   fChain->SetBranchAddress("LepInfo.GenPhi", LepInfo_GenPhi, &b_LepInfo_GenPhi);
   fChain->SetBranchAddress("LepInfo.GenPdgID", LepInfo_GenPdgID, &b_LepInfo_GenPdgID);
   fChain->SetBranchAddress("LepInfo.GenMCTag", LepInfo_GenMCTag, &b_LepInfo_GenMCTag);
   fChain->SetBranchAddress("LepInfo.TrgPt", LepInfo_TrgPt, &b_LepInfo_TrgPt);
   fChain->SetBranchAddress("LepInfo.TrgEta", LepInfo_TrgEta, &b_LepInfo_TrgEta);
   fChain->SetBranchAddress("LepInfo.TrgPhi", LepInfo_TrgPhi, &b_LepInfo_TrgPhi);
   fChain->SetBranchAddress("LepInfo.TrgID", LepInfo_TrgID, &b_LepInfo_TrgID);
   fChain->SetBranchAddress("VertexInfo.Size", &VertexInfo_Size, &b_VertexInfo_Size);
   fChain->SetBranchAddress("VertexInfo.isValid", VertexInfo_isValid, &b_VertexInfo_isValid);
   fChain->SetBranchAddress("VertexInfo.isFake", VertexInfo_isFake, &b_VertexInfo_isFake);
   fChain->SetBranchAddress("VertexInfo.Type", VertexInfo_Type, &b_VertexInfo_Type);
   fChain->SetBranchAddress("VertexInfo.Ndof", VertexInfo_Ndof, &b_VertexInfo_Ndof);
   fChain->SetBranchAddress("VertexInfo.NormalizedChi2", VertexInfo_NormalizedChi2, &b_VertexInfo_NormalizedChi2);
   fChain->SetBranchAddress("VertexInfo.Pt_Sum", VertexInfo_Pt_Sum, &b_VertexInfo_Pt_Sum);
   fChain->SetBranchAddress("VertexInfo.x", VertexInfo_x, &b_VertexInfo_x);
   fChain->SetBranchAddress("VertexInfo.y", VertexInfo_y, &b_VertexInfo_y);
   fChain->SetBranchAddress("VertexInfo.z", VertexInfo_z, &b_VertexInfo_z);
   fChain->SetBranchAddress("VertexInfo.Rho", VertexInfo_Rho, &b_VertexInfo_Rho);
   fChain->SetBranchAddress("PFJetInfo.Size", &PFJetInfo_Size, &b_PFJetInfo_Size);
   fChain->SetBranchAddress("PFJetInfo.Index", PFJetInfo_Index, &b_PFJetInfo_Index);
   fChain->SetBranchAddress("PFJetInfo.NTracks", PFJetInfo_NTracks, &b_PFJetInfo_NTracks);
   fChain->SetBranchAddress("PFJetInfo.Et", PFJetInfo_Et, &b_PFJetInfo_Et);
   fChain->SetBranchAddress("PFJetInfo.Pt", PFJetInfo_Pt, &b_PFJetInfo_Pt);
   fChain->SetBranchAddress("PFJetInfo.Unc", PFJetInfo_Unc, &b_PFJetInfo_Unc);
   fChain->SetBranchAddress("PFJetInfo.Eta", PFJetInfo_Eta, &b_PFJetInfo_Eta);
   fChain->SetBranchAddress("PFJetInfo.Phi", PFJetInfo_Phi, &b_PFJetInfo_Phi);
   fChain->SetBranchAddress("PFJetInfo.JetIDLOOSE", PFJetInfo_JetIDLOOSE, &b_PFJetInfo_JetIDLOOSE);
   fChain->SetBranchAddress("PFJetInfo.JetCharge", PFJetInfo_JetCharge, &b_PFJetInfo_JetCharge);
   fChain->SetBranchAddress("PFJetInfo.NConstituents", PFJetInfo_NConstituents, &b_PFJetInfo_NConstituents);
   fChain->SetBranchAddress("PFJetInfo.NCH", PFJetInfo_NCH, &b_PFJetInfo_NCH);
   fChain->SetBranchAddress("PFJetInfo.CEF", PFJetInfo_CEF, &b_PFJetInfo_CEF);
   fChain->SetBranchAddress("PFJetInfo.NHF", PFJetInfo_NHF, &b_PFJetInfo_NHF);
   fChain->SetBranchAddress("PFJetInfo.NEF", PFJetInfo_NEF, &b_PFJetInfo_NEF);
   fChain->SetBranchAddress("PFJetInfo.CHF", PFJetInfo_CHF, &b_PFJetInfo_CHF);
   fChain->SetBranchAddress("PFJetInfo.JVAlpha", PFJetInfo_JVAlpha, &b_PFJetInfo_JVAlpha);
   fChain->SetBranchAddress("PFJetInfo.JVBeta", PFJetInfo_JVBeta, &b_PFJetInfo_JVBeta);
   fChain->SetBranchAddress("PFJetInfo.PtCorrRaw", PFJetInfo_PtCorrRaw, &b_PFJetInfo_PtCorrRaw);
   fChain->SetBranchAddress("PFJetInfo.PtCorrL2", PFJetInfo_PtCorrL2, &b_PFJetInfo_PtCorrL2);
   fChain->SetBranchAddress("PFJetInfo.PtCorrL3", PFJetInfo_PtCorrL3, &b_PFJetInfo_PtCorrL3);
   fChain->SetBranchAddress("PFJetInfo.PtCorrL7g", PFJetInfo_PtCorrL7g, &b_PFJetInfo_PtCorrL7g);
   fChain->SetBranchAddress("PFJetInfo.PtCorrL7uds", PFJetInfo_PtCorrL7uds, &b_PFJetInfo_PtCorrL7uds);
   fChain->SetBranchAddress("PFJetInfo.PtCorrL7c", PFJetInfo_PtCorrL7c, &b_PFJetInfo_PtCorrL7c);
   fChain->SetBranchAddress("PFJetInfo.PtCorrL7b", PFJetInfo_PtCorrL7b, &b_PFJetInfo_PtCorrL7b);
   fChain->SetBranchAddress("PFJetInfo.JetBProbBJetTags", PFJetInfo_JetBProbBJetTags, &b_PFJetInfo_JetBProbBJetTags);
   fChain->SetBranchAddress("PFJetInfo.JetProbBJetTags", PFJetInfo_JetProbBJetTags, &b_PFJetInfo_JetProbBJetTags);
   fChain->SetBranchAddress("PFJetInfo.TrackCountHiPurBJetTags", PFJetInfo_TrackCountHiPurBJetTags, &b_PFJetInfo_TrackCountHiPurBJetTags);
   fChain->SetBranchAddress("PFJetInfo.TrackCountHiEffBJetTags", PFJetInfo_TrackCountHiEffBJetTags, &b_PFJetInfo_TrackCountHiEffBJetTags);
   fChain->SetBranchAddress("PFJetInfo.SimpleSVBJetTags", PFJetInfo_SimpleSVBJetTags, &b_PFJetInfo_SimpleSVBJetTags);
   fChain->SetBranchAddress("PFJetInfo.SimpleSVHEBJetTags", PFJetInfo_SimpleSVHEBJetTags, &b_PFJetInfo_SimpleSVHEBJetTags);
   fChain->SetBranchAddress("PFJetInfo.SimpleSVHPBJetTags", PFJetInfo_SimpleSVHPBJetTags, &b_PFJetInfo_SimpleSVHPBJetTags);
   fChain->SetBranchAddress("PFJetInfo.CombinedSVBJetTags", PFJetInfo_CombinedSVBJetTags, &b_PFJetInfo_CombinedSVBJetTags);
   fChain->SetBranchAddress("PFJetInfo.CombinedSVMVABJetTags", PFJetInfo_CombinedSVMVABJetTags, &b_PFJetInfo_CombinedSVMVABJetTags);
   fChain->SetBranchAddress("PFJetInfo.SoftElecByIP3dBJetTags", PFJetInfo_SoftElecByIP3dBJetTags, &b_PFJetInfo_SoftElecByIP3dBJetTags);
   fChain->SetBranchAddress("PFJetInfo.SoftElecByPtBJetTags", PFJetInfo_SoftElecByPtBJetTags, &b_PFJetInfo_SoftElecByPtBJetTags);
   fChain->SetBranchAddress("PFJetInfo.SoftMuonBJetTags", PFJetInfo_SoftMuonBJetTags, &b_PFJetInfo_SoftMuonBJetTags);
   fChain->SetBranchAddress("PFJetInfo.SoftMuonByIP3dBJetTags", PFJetInfo_SoftMuonByIP3dBJetTags, &b_PFJetInfo_SoftMuonByIP3dBJetTags);
   fChain->SetBranchAddress("PFJetInfo.SoftMuonByPtBJetTags", PFJetInfo_SoftMuonByPtBJetTags, &b_PFJetInfo_SoftMuonByPtBJetTags);
   fChain->SetBranchAddress("PFJetInfo.GenJetPt", PFJetInfo_GenJetPt, &b_PFJetInfo_GenJetPt);
   fChain->SetBranchAddress("PFJetInfo.GenJetEta", PFJetInfo_GenJetEta, &b_PFJetInfo_GenJetEta);
   fChain->SetBranchAddress("PFJetInfo.GenJetPhi", PFJetInfo_GenJetPhi, &b_PFJetInfo_GenJetPhi);
   fChain->SetBranchAddress("PFJetInfo.GenPt", PFJetInfo_GenPt, &b_PFJetInfo_GenPt);
   fChain->SetBranchAddress("PFJetInfo.GenEta", PFJetInfo_GenEta, &b_PFJetInfo_GenEta);
   fChain->SetBranchAddress("PFJetInfo.GenPhi", PFJetInfo_GenPhi, &b_PFJetInfo_GenPhi);
   fChain->SetBranchAddress("PFJetInfo.GenPdgID", PFJetInfo_GenPdgID, &b_PFJetInfo_GenPdgID);
   fChain->SetBranchAddress("PFJetInfo.GenFlavor", PFJetInfo_GenFlavor, &b_PFJetInfo_GenFlavor);
   fChain->SetBranchAddress("PFJetInfo.GenMCTag", PFJetInfo_GenMCTag, &b_PFJetInfo_GenMCTag);
   fChain->SetBranchAddress("PFJetInfo.Px", PFJetInfo_Px, &b_PFJetInfo_Px);
   fChain->SetBranchAddress("PFJetInfo.Py", PFJetInfo_Py, &b_PFJetInfo_Py);
   fChain->SetBranchAddress("PFJetInfo.Pz", PFJetInfo_Pz, &b_PFJetInfo_Pz);
   fChain->SetBranchAddress("PFJetInfo.Energy", PFJetInfo_Energy, &b_PFJetInfo_Energy);
   fChain->SetBranchAddress("JetInfo.Size", &JetInfo_Size, &b_JetInfo_Size);
   fChain->SetBranchAddress("JetInfo.Index", JetInfo_Index, &b_JetInfo_Index);
   fChain->SetBranchAddress("JetInfo.NTracks", JetInfo_NTracks, &b_JetInfo_NTracks);
   fChain->SetBranchAddress("JetInfo.Et", JetInfo_Et, &b_JetInfo_Et);
   fChain->SetBranchAddress("JetInfo.Pt", JetInfo_Pt, &b_JetInfo_Pt);
   fChain->SetBranchAddress("JetInfo.Unc", JetInfo_Unc, &b_JetInfo_Unc);
   fChain->SetBranchAddress("JetInfo.Eta", JetInfo_Eta, &b_JetInfo_Eta);
   fChain->SetBranchAddress("JetInfo.Phi", JetInfo_Phi, &b_JetInfo_Phi);
   fChain->SetBranchAddress("JetInfo.JetIDLOOSE", JetInfo_JetIDLOOSE, &b_JetInfo_JetIDLOOSE);
   fChain->SetBranchAddress("JetInfo.JetCharge", JetInfo_JetCharge, &b_JetInfo_JetCharge);
   fChain->SetBranchAddress("JetInfo.NConstituents", JetInfo_NConstituents, &b_JetInfo_NConstituents);
   fChain->SetBranchAddress("JetInfo.NCH", JetInfo_NCH, &b_JetInfo_NCH);
   fChain->SetBranchAddress("JetInfo.CEF", JetInfo_CEF, &b_JetInfo_CEF);
   fChain->SetBranchAddress("JetInfo.NHF", JetInfo_NHF, &b_JetInfo_NHF);
   fChain->SetBranchAddress("JetInfo.NEF", JetInfo_NEF, &b_JetInfo_NEF);
   fChain->SetBranchAddress("JetInfo.CHF", JetInfo_CHF, &b_JetInfo_CHF);
   fChain->SetBranchAddress("JetInfo.JVAlpha", JetInfo_JVAlpha, &b_JetInfo_JVAlpha);
   fChain->SetBranchAddress("JetInfo.JVBeta", JetInfo_JVBeta, &b_JetInfo_JVBeta);
   fChain->SetBranchAddress("JetInfo.PtCorrRaw", JetInfo_PtCorrRaw, &b_JetInfo_PtCorrRaw);
   fChain->SetBranchAddress("JetInfo.PtCorrL2", JetInfo_PtCorrL2, &b_JetInfo_PtCorrL2);
   fChain->SetBranchAddress("JetInfo.PtCorrL3", JetInfo_PtCorrL3, &b_JetInfo_PtCorrL3);
   fChain->SetBranchAddress("JetInfo.PtCorrL7g", JetInfo_PtCorrL7g, &b_JetInfo_PtCorrL7g);
   fChain->SetBranchAddress("JetInfo.PtCorrL7uds", JetInfo_PtCorrL7uds, &b_JetInfo_PtCorrL7uds);
   fChain->SetBranchAddress("JetInfo.PtCorrL7c", JetInfo_PtCorrL7c, &b_JetInfo_PtCorrL7c);
   fChain->SetBranchAddress("JetInfo.PtCorrL7b", JetInfo_PtCorrL7b, &b_JetInfo_PtCorrL7b);
   fChain->SetBranchAddress("JetInfo.JetBProbBJetTags", JetInfo_JetBProbBJetTags, &b_JetInfo_JetBProbBJetTags);
   fChain->SetBranchAddress("JetInfo.JetProbBJetTags", JetInfo_JetProbBJetTags, &b_JetInfo_JetProbBJetTags);
   fChain->SetBranchAddress("JetInfo.TrackCountHiPurBJetTags", JetInfo_TrackCountHiPurBJetTags, &b_JetInfo_TrackCountHiPurBJetTags);
   fChain->SetBranchAddress("JetInfo.TrackCountHiEffBJetTags", JetInfo_TrackCountHiEffBJetTags, &b_JetInfo_TrackCountHiEffBJetTags);
   fChain->SetBranchAddress("JetInfo.SimpleSVBJetTags", JetInfo_SimpleSVBJetTags, &b_JetInfo_SimpleSVBJetTags);
   fChain->SetBranchAddress("JetInfo.SimpleSVHEBJetTags", JetInfo_SimpleSVHEBJetTags, &b_JetInfo_SimpleSVHEBJetTags);
   fChain->SetBranchAddress("JetInfo.SimpleSVHPBJetTags", JetInfo_SimpleSVHPBJetTags, &b_JetInfo_SimpleSVHPBJetTags);
   fChain->SetBranchAddress("JetInfo.CombinedSVBJetTags", JetInfo_CombinedSVBJetTags, &b_JetInfo_CombinedSVBJetTags);
   fChain->SetBranchAddress("JetInfo.CombinedSVMVABJetTags", JetInfo_CombinedSVMVABJetTags, &b_JetInfo_CombinedSVMVABJetTags);
   fChain->SetBranchAddress("JetInfo.SoftElecByIP3dBJetTags", JetInfo_SoftElecByIP3dBJetTags, &b_JetInfo_SoftElecByIP3dBJetTags);
   fChain->SetBranchAddress("JetInfo.SoftElecByPtBJetTags", JetInfo_SoftElecByPtBJetTags, &b_JetInfo_SoftElecByPtBJetTags);
   fChain->SetBranchAddress("JetInfo.SoftMuonBJetTags", JetInfo_SoftMuonBJetTags, &b_JetInfo_SoftMuonBJetTags);
   fChain->SetBranchAddress("JetInfo.SoftMuonByIP3dBJetTags", JetInfo_SoftMuonByIP3dBJetTags, &b_JetInfo_SoftMuonByIP3dBJetTags);
   fChain->SetBranchAddress("JetInfo.SoftMuonByPtBJetTags", JetInfo_SoftMuonByPtBJetTags, &b_JetInfo_SoftMuonByPtBJetTags);
   fChain->SetBranchAddress("JetInfo.GenJetPt", JetInfo_GenJetPt, &b_JetInfo_GenJetPt);
   fChain->SetBranchAddress("JetInfo.GenJetEta", JetInfo_GenJetEta, &b_JetInfo_GenJetEta);
   fChain->SetBranchAddress("JetInfo.GenJetPhi", JetInfo_GenJetPhi, &b_JetInfo_GenJetPhi);
   fChain->SetBranchAddress("JetInfo.GenPt", JetInfo_GenPt, &b_JetInfo_GenPt);
   fChain->SetBranchAddress("JetInfo.GenEta", JetInfo_GenEta, &b_JetInfo_GenEta);
   fChain->SetBranchAddress("JetInfo.GenPhi", JetInfo_GenPhi, &b_JetInfo_GenPhi);
   fChain->SetBranchAddress("JetInfo.GenPdgID", JetInfo_GenPdgID, &b_JetInfo_GenPdgID);
   fChain->SetBranchAddress("JetInfo.GenFlavor", JetInfo_GenFlavor, &b_JetInfo_GenFlavor);
   fChain->SetBranchAddress("JetInfo.GenMCTag", JetInfo_GenMCTag, &b_JetInfo_GenMCTag);
   fChain->SetBranchAddress("JetInfo.Px", JetInfo_Px, &b_JetInfo_Px);
   fChain->SetBranchAddress("JetInfo.Py", JetInfo_Py, &b_JetInfo_Py);
   fChain->SetBranchAddress("JetInfo.Pz", JetInfo_Pz, &b_JetInfo_Pz);
   fChain->SetBranchAddress("JetInfo.Energy", JetInfo_Energy, &b_JetInfo_Energy);
   fChain->SetBranchAddress("PairInfo.Size", &PairInfo_Size, &b_PairInfo_Size);
   fChain->SetBranchAddress("PairInfo.Index", PairInfo_Index, &b_PairInfo_Index);
   fChain->SetBranchAddress("PairInfo.Type", PairInfo_Type, &b_PairInfo_Type);
   fChain->SetBranchAddress("PairInfo.Obj1Index", PairInfo_Obj1Index, &b_PairInfo_Obj1Index);
   fChain->SetBranchAddress("PairInfo.Obj2Index", PairInfo_Obj2Index, &b_PairInfo_Obj2Index);
   fChain->SetBranchAddress("PairInfo.Mass", PairInfo_Mass, &b_PairInfo_Mass);
   fChain->SetBranchAddress("PairInfo.Pt", PairInfo_Pt, &b_PairInfo_Pt);
   fChain->SetBranchAddress("PairInfo.Eta", PairInfo_Eta, &b_PairInfo_Eta);
   fChain->SetBranchAddress("PairInfo.Phi", PairInfo_Phi, &b_PairInfo_Phi);
   fChain->SetBranchAddress("PairInfo.GenMass", PairInfo_GenMass, &b_PairInfo_GenMass);
   fChain->SetBranchAddress("PairInfo.GenPt", PairInfo_GenPt, &b_PairInfo_GenPt);
   fChain->SetBranchAddress("PairInfo.GenEta", PairInfo_GenEta, &b_PairInfo_GenEta);
   fChain->SetBranchAddress("PairInfo.GenPhi", PairInfo_GenPhi, &b_PairInfo_GenPhi);
   fChain->SetBranchAddress("PairInfo.GenPdgID", PairInfo_GenPdgID, &b_PairInfo_GenPdgID);
   Notify();
}

Bool_t bPrimeNtupleMC::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void bPrimeNtupleMC::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t bPrimeNtupleMC::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef bPrimeNtupleMC_cxx
