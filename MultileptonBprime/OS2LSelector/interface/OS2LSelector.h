#ifndef OS2LSELECTOR_H  
#define OS2LSELECTOR_H  

#include <TFile.h>
#include <TTree.h>
#include <string>
#include <TChain.h>
#include <TH1D.h> 

struct Histograms ; 
class EvtInfoBranches ; 
class LepInfoBranches ; 

static const int nmaxZ =  15 ; 

class OS2LSelector {

  public: 

  OS2LSelector (const std::string& infiles=0,const bool& isData=0,const std::string&jsonfile=0,const std::string&outfile=0) ;	  
  ~OS2LSelector() ; 
  virtual void process () ;

  protected: 

  virtual bool isGoodElectron (EvtInfoBranches*, LepInfoBranches*, int&, int&) ;
  virtual bool isGoodMuon     (LepInfoBranches*, int&, int&) ;
  virtual bool isGoodTau      (LepInfoBranches*, int&, int&) ; 

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

  TH1D* m_hcutFlow_ ; 
  TH1D* m_hnvtx_PUWt_ ; 
  TH1D* m_hnvtx_noPUWt_ ; 
  TH1D* m_hnvtx_PUWt_HLT_ ; 
  TH1D* m_hnvtx_noPUWt_HLT_ ; 
  TH1D* m_hnvtx_PUWt_ZWin_ ; 
  TH1D* m_hnvtx_noPUWt_ZWin_ ; 
  TH1D* m_hnvtx_PUWt_ZVeto_ ; 
  TH1D* m_hnvtx_noPUWt_ZVeto_ ; 
  TH1D* m_hpte_ ; 
  TH1D* m_hetae_ ;
  TH1D* m_hptee_ ; 
  TH1D* m_hetaee_ ; 
  TH1D* m_hmee_ ; 
  TH1D* m_hptmu_ ; 
  TH1D* m_hetamu_ ;
  TH1D* m_hptmumu_ ; 
  TH1D* m_hetamumu_ ; 
  TH1D* m_hmmumu_ ; 
  TH1D* m_hptemu_ ; 
  TH1D* m_hetaemu_ ; 
  TH1D* m_hmemu_ ; 
  TH1D* m_hnjets_ ; 
  TH1D* m_hnbtags_ ; 

  Histograms* hist_ ; 

};

#endif
