#ifndef BPRIMECANDSELECTOR_H  
#define BPRIMECANDSELECTOR_H  

#include <TFile.h>
#include <TTree.h>
#include <string>
#include <TChain.h>
#include <TH1D.h> 
#include "MultileptonBprime/OS2LSelector/interface/OS2LSelector.h" 
#include <algorithm>
#include <string>

struct Histograms ; 
class BookHistograms ; 
class EvtInfoBranches ; 
class LepInfoBranches ; 

static const int nmax =  15 ; 
static const int NSTBins =  7 ; 

class BprimeCandSelector : public OS2LSelector {

  public: 

  BprimeCandSelector (const std::string& infiles=0,const bool& isData=0,const std::string&jsonfile=0,const std::string&outfile=0) ;	  
  ~BprimeCandSelector() ; 
  virtual void process () ;

  private:

  virtual void evtLoop () ;

  std::string inputfiles_ ;
  TString outfileName_ ; 
  TFile outfile_ ;
  bool isData_ ; 
  std::string jsonfile_ ; 
  TChain* thisChain_ ; 

  unsigned int evtProc_       ; 
  unsigned int evtAccept_     ; 
  unsigned int ngoodElectrons_; 
  unsigned int ngoodMuons_    ; 

  //TH1D* m_hcutFlow_ ; 
  //TH1D* m_hpte_ ; 
  //TH1D* m_hetae_ ;
  //TH1D* m_hptee_ ; 
  //TH1D* m_hetaee_ ; 
  //TH1D* m_hnjets_ ; 
  //TH1D* m_hnbtags_ ; 

//  Histograms* hist_ ; 
  BookHistograms* hist_ ; 
  BookHistograms* histSTLt500_ ; 
  BookHistograms* histSTGt500_ ; 

  BookHistograms* hist_DYee_DY_ ; 
  BookHistograms* hist_DYee_DY4Jets_ ; 
  BookHistograms* hist_DYee_DYSTLt500_ ; 
  BookHistograms* hist_DYee_DYSTGt500_ ; 
  BookHistograms* hist_DYee_DYSTGt500NoBTag_ ; 

  BookHistograms* hist_DYmumu_DY_ ; 
  BookHistograms* hist_DYmumu_DY4Jets_ ; 
  BookHistograms* hist_DYmumu_DYSTLt500_ ; 
  BookHistograms* hist_DYmumu_DYSTGt500_ ; 
  BookHistograms* hist_DYmumu_DYSTGt500NoBTag_ ; 

  BookHistograms* hist_EMu_DY_ ; 
  BookHistograms* hist_EMu_DY4Jets_ ; 
  BookHistograms* hist_EMu_DYSTLt500_ ; 
  BookHistograms* hist_EMu_DYSTGt500_ ; 
  BookHistograms* hist_EMu_DYSTGt500NoBTag_ ; 
  
  TH1D* m_hevtsel_      ; 
  TH1D* m_hevtsel_ee_   ; 
  TH1D* m_hevtsel_emu_  ; 
  TH1D* m_hevtsel_mumu_ ; 

  TH1D* m_hevtselSig_ee_  [100] ; 
  TH1D* m_hevtselSig_emu_ [100] ; 
  TH1D* m_hevtselSig_mumu_[100] ; 

  static double STBins[NSTBins] ; 

  TH1D* m_hst[3][3] ; 
  TH1D* m_hst_sig[3][3][66] ; 


};

#endif


